#include <math.h>
#include <itpp/itcomm.h>
#include "../../Library/Common/helper.hpp"
#include "../../Library/TrellisShaping/TrellisShaping.hpp"
#include "../../Library/TrellisShaping/TrellisCode.hpp"
#include "../../Library/TrellisShaping/QAM.hpp"
#include "../../Library/ClippingAndFiltering/ClippingAndFiltering.hpp"

using std::max;
using itpp::pow2;
using itpp::dec2bin;
using itpp::linspace;
using itpp::round_i;
using itpp::bin;
using itpp::vec;
using itpp::bvec;
using itpp::ivec;
using itpp::cvec;
using itpp::mat;
using helper::is_power_of_4;
using helper::demultiplex;
using helper::multiplex;
using helper::hard_decide_LLRs;
using helper::INFTY;
using helper::is_infinity;
using complex = std::complex<double>;


const static mat clippingRatioTable(
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0   0.1   0.2   0.3   0.6;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0   0.5   0.9   1.2   1.3;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0   0.9   1.3   1.7;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0"
);

TrellisShaping::TrellisShaping(const int m, const TSMetricType metric, const TSMappingType mappingType, const TSTerminationType terminationType)
{
    this->metric = metric;
    this->rate = (log2(m) - 1.0) / log2(m);

    this->set_constellation(m, mappingType);
    this->set_shaping_method(metric);
    this->set_termination_type(terminationType);
    this->set_convolutional_code();

    this->isValid = true;
}

TrellisShaping::~TrellisShaping()
{
}

double TrellisShaping::get_rate()
{
    return this->rate;
}

double TrellisShaping::get_optimum_clipping_ratio(const double scalingFactor)
{
    double clippingRatio = clippingRatioTable(this->QAMModulator.get_k(), round_i(10 * scalingFactor));

    it_assert_debug(clippingRatio > 0, "TrellisShaping::get_optimum_clipping_ratio(): Undefined parameters set is designated.");

    return clippingRatio;
}

void TrellisShaping::set_shaping_method(TSMetricType metricType)
{
    switch (metricType) {
        case TSMetricType::Autocorrelation:
            this->shapingMethod = TSShapingMethod::Sign;
            this->isSetUp = true;
            break;
        case TSMetricType::Clipping:
            this->shapingMethod = TSShapingMethod::Tail;
            break;
        default:
            it_error("Invalid shaping method is given.");
    };
}

void TrellisShaping::set_termination_type(const TSTerminationType terminationType)
{
    switch (terminationType) {
        case TSTerminationType::Unterminated:
            break;
        case TSTerminationType::Tailbite:
            break;
        default:
            it_error("TrellisCode::set_termination_type(): Invalid type of termination.");
    }
    this->terminationType = terminationType;
}

void TrellisShaping::set_convolutional_code()
{
    this->generator.set_type(TrellisCodeType::Generator);
    this->inverseSymdrome.set_type(TrellisCodeType::InverseSyndrome);
    this->parityChecker.set_type(TrellisCodeType::ParityChecker);
    this->compoundDecoder.set_type(TrellisCodeType::CompoundDecoder);

    this->numberOfShapingBits = 2;
}

void TrellisShaping::set_constellation(const int m, const TSMappingType mappingType)
{
    it_assert_debug(is_power_of_4(m) && m > 4, "TrellisShaping::set_constellation(): Invalid size of constellation.");

    this->mappingType = mappingType;
    this->QAMModulator.set_constellation(m, mappingType);

    if (this->metric == TSMetricType::Autocorrelation) {
        this->QAMModulator.make_correlation_table();
    }
}


bvec TrellisShaping::add_shaping_bits(const bvec& codeBits, const bvec& shapingBits)
{
    it_assert_debug(codeBits.size() == this->QAMModulator.get_k(), "TrellisShaping::add_shaping_bits(): Invalid code length of given vector.");
    it_assert_debug(shapingBits.size() == this->numberOfShapingBits, "TrellisShaping::add_shaping_bits(): Invalid code length of given vector.");

    bvec shapedBits(codeBits);
    if (this->shapingMethod == TSShapingMethod::Sign) {
        bvec MSBs = codeBits.left(this->numberOfShapingBits) + shapingBits;
        shapedBits.set_subvector(0, MSBs);
    } else if (this->shapingMethod == TSShapingMethod::Tail) {
        bvec LSBs = codeBits.right(this->numberOfShapingBits) + shapingBits;
        shapedBits.set_subvector(shapedBits.size() - LSBs.size(), LSBs);
    }
    return shapedBits;
}

double TrellisShaping::calculate_metric(const cvec& symbols, const cvec& autocorrelation, const int stage)
{
    double branchMetric = 0;

    switch (this->metric) {
        case TSMetricType::Autocorrelation:
            for (int m = 1; m < stage; m++) {
                branchMetric += 2.0 * real(conj(autocorrelation[m]) * this->QAMModulator.get_correlation(symbols[stage], symbols[stage - m]));
            }
            if (this->mappingType == TSMappingType::Type2) {
                for (int m = 1; m <= stage; m++) {
                    branchMetric += this->QAMModulator.get_squared_correlation(symbols[stage], symbols[stage - m]);
                }
            }
            break;
        case TSMetricType::Clipping:
            branchMetric += norm(symbols[stage] - this->referenceSymbols[stage]);
            break;
        default:
            it_error("Invalid metric type is given.");
    }

    return branchMetric;
}


double TrellisShaping::calculate_transit_probability(const int currentState, int nextState, const complex& symbol, const double N0)
{
    const int constellationSize = this->QAMModulator.get_M();
    cvec constellations(constellationSize / 2);
    bvec input = dec2bin(this->numberOfShapingBits, this->compoundDecoder.get_input(currentState, nextState));

    double probability = 1;
    for (int i = 0; i < this->numberOfShapingBits; i++) {
        int index;
        if (this->shapingMethod == TSShapingMethod::Sign) {
            index = i;
        } else {
            index = this->QAMModulator.get_k() - this->numberOfShapingBits + i;
        }

        constellations.set_subvector(0, this->QAMModulator.get_constellation_points(index, input[i]));

        double sum = 0;
        for (int j = 0; j < constellations.size(); j++) {
            sum += (1.0 / constellationSize) / sqrt(M_PI * N0) * exp(- norm(symbol - constellations[j]) / N0);
        }
        probability *= sum;
    }

    return log(probability);
}

void TrellisShaping::make_reference_symbols(const bvec& codes, double scalingFactor, double clippingRatio)
{
    const int constellationSize = this->QAMModulator.get_M() / 4;
    const int blockLength = this->QAMModulator.get_k();
    clippingRatio = this->get_optimum_clipping_ratio(scalingFactor);
    ClippingAndFiltering CAFProccesser(clippingRatio);
    bvec MSBs, LSBs;

    it_assert_debug(this->isValid, "TrellisShaping::make_reference_symbols(): Constellation and type are not set.");
    it_assert_debug(this->metric == TSMetricType::Clipping, "TrellisShaping::make_reference_symbols(): This is only for clipping metric.");
    it_assert_debug(codes.size() % blockLength == 0, "TrellisShaping::make_reference_symbols(): Invalid data size of given vecotr.");

    demultiplex(codes, MSBs, LSBs, blockLength - this->numberOfShapingBits, this->numberOfShapingBits);

    // Save original QAM constellation temporarily.
    this->QAMModulator.push_constellation();

    // Set a new QAM constellation for reference symbols.
    scalingFactor *= sqrt(pow2(this->numberOfShapingBits) * (constellationSize - 1) / (constellationSize * 4 - 1));
    this->QAMModulator.set_constellation(constellationSize, this->mappingType, scalingFactor);

    // Make the replica symbols sequence by clipping and set as reference symbols.
    this->referenceSymbols = CAFProccesser.clip_and_filter(this->QAMModulator.modulate_bits(MSBs));
    this->referenceSymbols = CAFProccesser.normalize_attenuation(this->referenceSymbols);

    this->QAMModulator.pop_constellation();
    this->isSetUp = true;
}

bvec TrellisShaping::form_syndrome(const bvec& outerCodes)
{
    int numberOfNonshapingBits = this->QAMModulator.get_k() - this->numberOfShapingBits;
    bvec MSBs, LSBs;
    bvec innerCodes;

    it_assert_debug(outerCodes.size() % (this->QAMModulator.get_k() - 1) == 0, "TrellisShaping::form_syndrome(): Invalid data size of given vecotr.");

    if (this->shapingMethod == TSShapingMethod::Sign) {
        demultiplex(outerCodes, MSBs, LSBs, 1, numberOfNonshapingBits);

        // Encode with inverse symdrome.
        MSBs = this->inverseSymdrome.encode(MSBs);

        multiplex(MSBs, LSBs, innerCodes, this->numberOfShapingBits, numberOfNonshapingBits);
    } else {
        demultiplex(outerCodes, MSBs, LSBs, numberOfNonshapingBits, 1);

        // Encode with inverse symdrome.
        LSBs = this->inverseSymdrome.encode(LSBs);

        multiplex(MSBs, LSBs, innerCodes, numberOfNonshapingBits, this->numberOfShapingBits);
    }

    it_assert_debug(innerCodes.size() % this->QAMModulator.get_k() == 0, "TrellisShaping::form_syndrome(): Invalid data size of given vecotr.");

    return innerCodes;
}

cvec TrellisShaping::trellis_coded_modulate(const bvec& codes)
{
    switch (this->terminationType) {
        case TSTerminationType::Unterminated:
            return encode_unterminated(codes);
        case TSTerminationType:: Tailbite:
            return encode_tailbite(codes);
        default:
            it_error("Invalid termination method is given.");
    }
}

cvec TrellisShaping::encode_unterminated(const bvec& codes)
{
    int blockLength = this->QAMModulator.get_k();
    int dataSize = codes.size() / blockLength;
    bvec tempOriginalBits(blockLength);
    bvec tempShapingBits(this->numberOfShapingBits);
    cvec tempSymbols;
    cvec tempAutocorrelation;

    it_assert_debug(this->isValid, "TrellisShaping::encode_unterminated(): Constellation and type are not set.");
    it_assert_debug(this->isSetUp, "TrellisShaping::encode_unterminated(): Setup has still not been completed.");
    it_assert_debug(dataSize == this->referenceSymbols.size(), "TrellisShaping::encode_unterminated(): Invalid size of given vector.");

    this->generator.init_trellis(dataSize);
    this->generator.nodeMetrics(0, 0) = 0;

    for (int node = 0; node < dataSize; node++) {
        tempSymbols.set_size(node + 1);
        tempAutocorrelation.set_size(node + 1);

        for (int currentState = 0; currentState < this->generator.get_number_of_states(); currentState++) {
            if (is_infinity(this->generator.nodeMetrics(node, currentState))) {
                continue;
            }

            // Get temporarily symbols sequence.
            tempSymbols.set_subvector(0, this->generator.get_symbols(node, currentState));
            // tempAutocorrelation.set_subvector(0, this->generator.get_autocorrelation(node, currentState));

            ivec nextStates = this->generator.get_next_states(currentState);

            for (int i = 0; i < nextStates.size(); i++) {
                // Make and append a new shped symbol.
                tempOriginalBits.set_subvector(0, codes.mid(node * blockLength, blockLength));
                tempShapingBits.set_subvector(0, dec2bin(this->numberOfShapingBits, this->generator.get_output(currentState, nextStates[i])));
                tempSymbols[node] = this->QAMModulator.modulate_bits(this->add_shaping_bits(tempOriginalBits, tempShapingBits)).get(0);

                // Calculate branch metric.
                double branchMetric = this->calculate_metric(tempSymbols, tempAutocorrelation, node);

                if (this->generator.nodeMetrics(node, currentState) + branchMetric <= this->generator.nodeMetrics(node + 1, nextStates[i])) {
                    // Update the metric.
                    this->generator.nodeMetrics(node + 1, nextStates[i]) = this->generator.nodeMetrics(node, currentState) + branchMetric;

                    // Reserve previous state and selected outputs.
                    this->generator.pathMemories(node + 1, nextStates[i]) = currentState;
                    this->generator.outputSymbols(node + 1, nextStates[i]) = tempSymbols[node];
                }
            }
        }
        //
        // if (this->metric == TSMetricType::Autocorrelation) {
        //     this->generator.autocorrelation(node + 1, nextStates[i]) = this->QAMModulator.get_correlation(tempSymbols[node - 1] , tempSymbols[node]);
        // }
    }

    // Find minimum metric.
    int state = this->generator.find_optimum_state(dataSize);

    // Trace back to calculate the output symbol.
    return this->generator.get_symbols(dataSize, state);
}

cvec TrellisShaping::encode_tailbite(const bvec& codes)
{
    int blockLength = this->QAMModulator.get_k();
    int dataSize = codes.size() / blockLength;
    int memorySize = this->generator.get_memory_size();
    bvec tempOriginalBits(blockLength);
    bvec tempShapingBits(this->numberOfShapingBits);
    cvec tempSymbols;
    cvec tempAutocorrelation;

    it_assert_debug(this->isValid, "TrellisShaping::encode_tailbite(): Constellation and type are not set.");
    it_assert_debug(this->isSetUp, "TrellisShaping::encode_tailbite(): Setup has still not been completed.");
    it_assert_debug(this->referenceSymbols.size() == 0 || dataSize == this->referenceSymbols.size(), "TrellisShaping::encode_tailbite(): Invalid size of given vector.");

    this->generator.init_trellis(dataSize);
    this->generator.nodeMetrics(0, 0) = 0;

    for (int node = 0; node < dataSize - memorySize; node++) {
        tempSymbols.set_size(node + 1);
        tempAutocorrelation.set_size(node + 1);

        for (int currentState = 0; currentState < this->generator.get_number_of_states(); currentState++) {
            if (is_infinity(this->generator.nodeMetrics(node, currentState))) {
                continue;
            }

            // Get temporarily symbols sequence.
            tempSymbols.set_subvector(0, this->generator.get_symbols(node, currentState));
            // tempAutocorrelation.set_subvector(0, this->generator.get_autocorrelation(node, currentState));

            ivec nextStates = this->generator.get_next_states(currentState);

            for (int i = 0; i < nextStates.size(); i++) {
                // Make and append a new shped symbol.
                tempOriginalBits.set_subvector(0, codes.mid(node * blockLength, blockLength));
                tempShapingBits.set_subvector(0, dec2bin(this->numberOfShapingBits, this->generator.get_output(currentState, nextStates[i])));
                tempSymbols[node] = this->QAMModulator.modulate_bits(this->add_shaping_bits(tempOriginalBits, tempShapingBits)).get(0);

                // Calculate branch metric.
                double branchMetric = this->calculate_metric(tempSymbols, tempAutocorrelation, node);

                if (this->generator.nodeMetrics(node, currentState) + branchMetric <= this->generator.nodeMetrics(node + 1, nextStates[i])) {
                    // Update the metric.
                    this->generator.nodeMetrics(node + 1, nextStates[i]) = this->generator.nodeMetrics(node, currentState) + branchMetric;

                    // Reserve previous state and selected outputs.
                    this->generator.pathMemories(node + 1, nextStates[i]) = currentState;
                    this->generator.outputSymbols(node + 1, nextStates[i]) = tempSymbols[node];
                }
            }
        }
    }

    // Find minimum metric.
    int state = this->generator.find_optimum_state(dataSize - memorySize);

    // Trace back to calculate the output symbol.
    cvec symbols = this->generator.get_symbols(dataSize - memorySize, state);

    // Append original
    for (int node = dataSize - memorySize; node < dataSize; node++) {
        tempOriginalBits.set_subvector(0, codes.mid(node * blockLength, blockLength));
        int nextState = this->generator.get_next_state(state, 0);
        tempShapingBits.set_subvector(0, dec2bin(this->numberOfShapingBits, this->generator.get_output(state, nextState)));
        symbols.ins(symbols.size(), this->QAMModulator.modulate_bits(this->add_shaping_bits(tempOriginalBits, tempShapingBits)).get(0));

        state = nextState;
    }

    return symbols;
}

bvec TrellisShaping::demodulate_bits(const cvec& symbols, const double N0)
{
    int numberOfNonshapingBits = this->QAMModulator.get_k() - this->numberOfShapingBits;
    bvec MSBs, LSBs;
    bvec outerCodes;

    bvec innerCodes = hard_decide_LLRs(this->QAMModulator.demodulate_soft_bits(symbols, N0));

    it_assert_debug(innerCodes.size() % this->QAMModulator.get_k() == 0, "TrellisShaping::form_syndrome(): Invalid data size of given vecotr.");

    if (this->shapingMethod == TSShapingMethod::Sign) {
        demultiplex(innerCodes, MSBs, LSBs, this->numberOfShapingBits, numberOfNonshapingBits);

        // Decode with inverse symdrome.
        MSBs = this->parityChecker.encode(MSBs);

        multiplex(MSBs, LSBs, outerCodes, 1, numberOfNonshapingBits);
    } else {
        demultiplex(innerCodes, MSBs, LSBs, numberOfNonshapingBits, this->numberOfShapingBits);

        // Decode with inverse symdrome.
        LSBs = this->parityChecker.encode(LSBs);

        multiplex(MSBs, LSBs, outerCodes, numberOfNonshapingBits, 1);
    }

    it_assert_debug(outerCodes.size() % (this->QAMModulator.get_k() - 1) == 0, "TrellisShaping::form_syndrome(): Invalid data size of given vecotr.");

    return outerCodes;
}

vec TrellisShaping::demodulate_soft_bits(const cvec& symbols, const double N0)
{
    switch (this->terminationType) {
        case TSTerminationType::Unterminated:
            return decode_by_SOVA(symbols, N0);
        case TSTerminationType:: Tailbite:
            return decode_by_BCJR(symbols, N0);
        default:
            it_error("Invalid termination method is given.");
    }
}

vec TrellisShaping::decode_by_SOVA(const cvec& symbols, const double N0)
{
    int dataSize = symbols.size();
    int blockLength = this->QAMModulator.get_k();
    int numberOfNonshapingBits = blockLength - this->numberOfShapingBits;
    vec LLRs, nonshapedLLRs, shapedLLRs;

    it_assert_debug(this->isValid, "TrellisShaping::decode_by_SOVA(): Constellation and type are not set.");
    it_assert_debug(this->isSetUp, "TrellisShaping::decode_by_SOVA(): Setup has still not been completed.");
    it_assert_debug(dataSize == this->referenceSymbols.size(), "TrellisShaping::decode_by_SOVA(): Invalid size of given vector.");

    // Caluclulate raw LLRs for non-shaped bits.
    LLRs = this->QAMModulator.demodulate_soft_bits(symbols, N0);

    if (this->shapingMethod == TSShapingMethod::Sign) {
        demultiplex(LLRs, shapedLLRs, nonshapedLLRs, this->numberOfShapingBits, numberOfNonshapingBits);
    } else {
        demultiplex(LLRs, nonshapedLLRs, shapedLLRs, numberOfNonshapingBits, this->numberOfShapingBits);
    }

    // Start max-log-MAP algorithm for shaped bits.
    this->compoundDecoder.init_trellis(dataSize);
    this->compoundDecoder.forwardMetrics(0, 0) = log(1);

    // First, calculate joint probabilities for forward recursion.
    for (int node = 1; node < dataSize; node++) {
        for (int currentState = 0; currentState < this->compoundDecoder.get_number_of_states(); currentState++) {
            ivec previousStates = this->compoundDecoder.get_previous_states(currentState);

            for (int i = 0; i < previousStates.size(); i++) {
                if (is_infinity(this->compoundDecoder.forwardMetrics(node - 1, previousStates[i]))) {
                    continue;
                }

                double metric = this->calculate_transit_probability(previousStates[i], currentState, symbols[node - 1], N0)
                    + this->compoundDecoder.forwardMetrics(node - 1, previousStates[i]);

                this->compoundDecoder.forwardMetrics(node, currentState) = max(this->compoundDecoder.forwardMetrics(node, currentState), metric);
            }
        }
    }

    // Next, obtain LLR of each shaped bit.
    shapedLLRs.set_size(dataSize * 1);

    for (int node = 0; node < dataSize; node++) {
        double metricFor0 = - INFTY;
        double metricFor1 = - INFTY;

        for (int currentState = 0; currentState < this->compoundDecoder.get_number_of_states(); currentState++) {
            ivec nextStates = this->compoundDecoder.get_next_states(currentState);

            for (int i = 0; i < nextStates.size(); i++) {
                int output = this->compoundDecoder.get_output(currentState, nextStates[i]);

                double metric = this->calculate_transit_probability(currentState, nextStates[i], symbols[node], N0)
                    + this->compoundDecoder.forwardMetrics(node, currentState);

                if (output == 0) {
                    metricFor0 = max(metricFor0, metric);
                } else if (output == 1) {
                    metricFor1 = max(metricFor1, metric);
                }
            }
        }
        shapedLLRs[node] = metricFor0 - metricFor1;
    }

    if (this->shapingMethod == TSShapingMethod::Sign) {
        multiplex(shapedLLRs, nonshapedLLRs, LLRs, 1, numberOfNonshapingBits);
    } else {
        multiplex(nonshapedLLRs, shapedLLRs, LLRs, numberOfNonshapingBits, 1);
    }

    return LLRs;
}

vec TrellisShaping::decode_by_BCJR(const cvec& symbols, const double N0)
{
    int dataSize = symbols.size();
    int blockLength = this->QAMModulator.get_k();
    int numberOfNonshapingBits = blockLength - this->numberOfShapingBits;
    vec LLRs, nonshapedLLRs, shapedLLRs;

    it_assert_debug(this->isValid, "TrellisShaping::decode_by_BCJR(): Constellation and type are not set.");
    it_assert_debug(this->isSetUp, "TrellisShaping::decode_by_BCJR(): Setup has still not been completed.");
    it_assert_debug(dataSize == this->referenceSymbols.size(), "TrellisShaping::decode_by_BCJR(): Invalid size of given vector.");

    // Caluclulate raw LLRs for non-shaped bits.
    LLRs = this->QAMModulator.demodulate_soft_bits(symbols, N0);

    if (this->shapingMethod == TSShapingMethod::Sign) {
        demultiplex(LLRs, shapedLLRs, nonshapedLLRs, this->numberOfShapingBits, numberOfNonshapingBits);
    } else {
        demultiplex(LLRs, nonshapedLLRs, shapedLLRs, numberOfNonshapingBits, this->numberOfShapingBits);
    }

    // Start max-log-MAP algorithm for shaped bits.
    this->compoundDecoder.init_trellis(dataSize);
    this->compoundDecoder.forwardMetrics(0, 0) = log(1);
    this->compoundDecoder.backwardMetrics(dataSize, 0) = log(1.0 / 2.0);
    this->compoundDecoder.backwardMetrics(dataSize, 4) = log(1.0 / 2.0);

    // First, calculate joint probabilities for forward recursion.
    for (int node = 1; node < dataSize; node++) {
        for (int currentState = 0; currentState < this->compoundDecoder.get_number_of_states(); currentState++) {
            ivec previousStates = this->compoundDecoder.get_previous_states(currentState);

            for (int i = 0; i < previousStates.size(); i++) {
                if (is_infinity(this->compoundDecoder.forwardMetrics(node - 1, previousStates[i]))) {
                    continue;
                }

                double metric = this->calculate_transit_probability(previousStates[i], currentState, symbols[node - 1], N0)
                    + this->compoundDecoder.forwardMetrics(node - 1, previousStates[i]);

                this->compoundDecoder.forwardMetrics(node, currentState) = max(this->compoundDecoder.forwardMetrics(node, currentState), metric);
            }
        }
    }

    // Second, calculate conditional probabilities for backward recursion.
    for (int node = dataSize - 1; node > 0; node--) {
        for (int currentState = 0; currentState < this->compoundDecoder.get_number_of_states(); currentState++) {
            ivec nextStates = this->compoundDecoder.get_next_states(currentState);

            for (int i = 0; i < nextStates.size(); i++) {
                if (is_infinity(this->compoundDecoder.backwardMetrics(node + 1, nextStates[i]))) {
                    continue;
                }

                double metric = this->calculate_transit_probability(currentState, nextStates[i], symbols[node], N0)
                    + this->compoundDecoder.backwardMetrics(node + 1, nextStates[i]);

                this->compoundDecoder.backwardMetrics(node, currentState) = max(this->compoundDecoder.backwardMetrics(node, currentState), metric);
            }
        }
    }

    // Finally, obtain LLR of each shaped bit.
    shapedLLRs.set_size(dataSize * 1);

    for (int node = 0; node < dataSize; node++) {
        double metricFor0 = - INFTY;
        double metricFor1 = - INFTY;

        for (int currentState = 0; currentState < this->compoundDecoder.get_number_of_states(); currentState++) {
            ivec nextStates = this->compoundDecoder.get_next_states(currentState);

            for (int i = 0; i < nextStates.size(); i++) {
                int output = this->compoundDecoder.get_output(currentState, nextStates[i]);

                double metric = this->calculate_transit_probability(currentState, nextStates[i], symbols[node], N0)
                    + this->compoundDecoder.forwardMetrics(node, currentState)
                    + this->compoundDecoder.backwardMetrics(node + 1, nextStates[i]);

                if (output == 0) {
                    metricFor0 = max(metricFor0, metric);
                } else if (output == 1) {
                    metricFor1 = max(metricFor1, metric);
                }
            }
        }
        shapedLLRs[node] = metricFor0 - metricFor1;
    }

    if (this->shapingMethod == TSShapingMethod::Sign) {
        multiplex(shapedLLRs, nonshapedLLRs, LLRs, 1, numberOfNonshapingBits);
    } else {
        multiplex(nonshapedLLRs, shapedLLRs, LLRs, numberOfNonshapingBits, 1);
    }

    return LLRs;
}
