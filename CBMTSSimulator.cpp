#include <iostream>
#include <iomanip>
#include <fstream>
#include <itpp/itcomm.h>
#include "config.hpp"
#include "./Library/Common/helper.hpp"
#include "./Library/Histogram/Histogram.hpp"
#include "./Library/Interleaver/RandomInterleaver.hpp"
#include "./Library/Interleaver/SRandomInterleaver.hpp"
#include "./Library/TrellisShaping/TrellisShaping.hpp"
#include "./Library/TurboCode/PuncturedTurboCode.hpp"


using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::ofstream;
using std::ios;
using itpp::GlobalRNG_randomize;
using itpp::vec;
using itpp::bvec;
using itpp::ivec;
using itpp::cvec;
using itpp::bmat;
using itpp::randb;
using itpp::linspace;
using itpp::pow2i;
using itpp::floor_i;
using itpp::dB;
using itpp::inv_dB;
using itpp::AWGN_Channel;
using itpp::BERC;
using itpp::BLERC;
using helper::join;
using helper::oversample;
using helper::downsample;
using helper::fft_normalized;
using helper::ifft_normalized;
using helper::assume_PAPR;
using helper::get_average_power;
using helper::hard_decide_LLRs;


void calculate_instantaneous_CCDF()
{
    const int constellationSize = pow2i(NUMBER_OF_BITS + 1);
    bvec transmittedBits;
    bvec transmittedCodes;
    cvec transmittedSymbols;
    cvec transmittedSignals;
    Histogram PSDHistogram(0, 13, 1001);
    TrellisShaping TSModulator(constellationSize, TSMetricType::Clipping, TSMappingType(TS_MAPPING_TYPE), TSTerminationType(TS_TERMINATION_TYPE));

    ofstream file(join("./Result/CBMTS/CCDF_Norm_", constellationSize, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_β-", SCALING_FACTOR, ".dat"), ios::out);

    for (int i = 0; i < NUMBER_OF_TRIALS; i++) {
        cout << "trial: " << i + 1 << "\r";

        transmittedBits = randb(NUMBER_OF_BITS * NUMBER_OF_SUBCARRIERS);

        transmittedCodes = TSModulator.form_syndrome(transmittedBits);

        TSModulator.make_reference_symbols(transmittedCodes, SCALING_FACTOR, CLIPPING_RATIO);

        transmittedSymbols = TSModulator.trellis_coded_modulate(transmittedCodes);

        transmittedSignals = ifft_normalized(oversample(transmittedSymbols, OVERSAMPLING_FACTOR));

        PSDHistogram.read_normalized_power(transmittedSignals);
    }

    cout << "\n";
    PSDHistogram.print_CCDF(file);
}

void calculate_average_power_reduction_capability()
{
    const int constellationSize = pow2i(NUMBER_OF_BITS + 1);
    bvec transmittedBits;
    bvec transmittedCodes;
    cvec transmittedSymbols;
    cvec transmittedSignals;
    double averagePowerOFDM, averagePowerTS;
    double capability = 0;
    itpp::QAM QAMModulator(constellationSize);
    TrellisShaping TSModulator(constellationSize, TSMetricType::Clipping, TSMappingType(TS_MAPPING_TYPE), TSTerminationType(TS_TERMINATION_TYPE));

    ofstream file(join("./Result/CBMTS/APR_Capa_", constellationSize, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_β-", SCALING_FACTOR, ".dat"), ios::out);

    for (int i = 0; i < NUMBER_OF_TRIALS; i++) {
        cout << "trial: " << i + 1 << "\r";

        // OFDM
        transmittedBits = randb((NUMBER_OF_BITS + 1) * NUMBER_OF_SUBCARRIERS);

        transmittedSymbols = QAMModulator.modulate_bits(transmittedBits);

        transmittedSignals = ifft_normalized(oversample(transmittedSymbols, OVERSAMPLING_FACTOR));

        averagePowerOFDM = get_average_power(transmittedSignals);

        // TS
        transmittedBits = randb(NUMBER_OF_BITS * NUMBER_OF_SUBCARRIERS);

        transmittedCodes = TSModulator.form_syndrome(transmittedBits);

        TSModulator.make_reference_symbols(transmittedCodes, SCALING_FACTOR, CLIPPING_RATIO);

        transmittedSymbols = TSModulator.trellis_coded_modulate(transmittedCodes);

        transmittedSignals = ifft_normalized(oversample(transmittedSymbols, OVERSAMPLING_FACTOR));

        averagePowerTS = get_average_power(transmittedSignals);

        capability += dB(averagePowerOFDM / averagePowerTS) / (double)(NUMBER_OF_TRIALS);
    }

    cout << capability << "[dB]" << "\n";
    file << capability << endl;
}

void calculate_clipping_ratio_character()
{
    const int constellationSize = pow2i(NUMBER_OF_BITS + 1);
    bvec transmittedBits;
    bvec transmittedCodes;
    cvec transmittedSymbols;
    cvec transmittedSignals;
    TrellisShaping TSModulator(constellationSize, TSMetricType::Clipping, TSMappingType(TS_MAPPING_TYPE), TSTerminationType(TS_TERMINATION_TYPE));
    double PAPR;
    ofstream file(join("./Result/CBMTS/Clipping_Chara_", constellationSize, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_β-", SCALING_FACTOR, ".dat"), ios::out);

    for (double  clippingRatio = 2.0; clippingRatio > 0; clippingRatio -= 0.1) {
        PAPR = 0;

        for (int i = 0; i < NUMBER_OF_TRIALS; i++) {
            cout << "clipping ratio: " << clippingRatio << " trial: " << i + 1 << " PAPR: " << PAPR / (i + 1) << "\r";

            transmittedBits = randb(NUMBER_OF_BITS * NUMBER_OF_SUBCARRIERS);

            transmittedCodes = TSModulator.form_syndrome(transmittedBits);

            TSModulator.make_reference_symbols(transmittedCodes, SCALING_FACTOR, clippingRatio);

            transmittedSymbols = TSModulator.trellis_coded_modulate(transmittedCodes);

            transmittedSignals = ifft_normalized(oversample(transmittedSymbols, OVERSAMPLING_FACTOR));

            PAPR += assume_PAPR(transmittedSignals);
        }

        cout << "\n";
        file << clippingRatio << " " << PAPR / NUMBER_OF_TRIALS << endl;
    }
}

void calculate_uncoded_BER_over_AWGN()
{
    const int constellationSize = pow2i(NUMBER_OF_BITS + 1);
    bvec transmittedBits;
    bvec transmittedCodes;
    cvec transmittedSymbols;
    cvec transmittedSignals;
    cvec receivedSignals;
    cvec receivedSymbols;
    vec receivedLLRs;
    bvec receivedBits;
    TrellisShaping TSModulator(constellationSize, TSMetricType::Clipping, TSMappingType(TS_MAPPING_TYPE), TSTerminationType(TS_TERMINATION_TYPE));
    AWGN_Channel AWGNCannel;
    BERC BERCounter;
    vec EbN0s = linspace(6, 22, 17);
    vec SNRs = TSModulator.get_rate() * log2(constellationSize) * inv_dB(EbN0s);
    ofstream file(join("./Result/CBMTS/BER_Uncod_", constellationSize, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_β-", SCALING_FACTOR, ".dat"), ios::out);

    for (int i = 0, size = EbN0s.size(); i < size; i++) {
        BERCounter.clear();

        for (int j = 0; j < NUMBER_OF_TRIALS; j++) {
            cerr << "trial: " << setw(5) << j + 1 << " Eb/N0(dB): " << setw(2) << EbN0s[i] << " BER: " << BERCounter.get_errorrate() << "\r";

            transmittedBits = randb(NUMBER_OF_BITS * NUMBER_OF_SUBCARRIERS);

            transmittedCodes = TSModulator.form_syndrome(transmittedBits);
            TSModulator.make_reference_symbols(transmittedCodes, SCALING_FACTOR, CLIPPING_RATIO);
            transmittedSymbols = TSModulator.trellis_coded_modulate(transmittedCodes);

            transmittedSignals = ifft_normalized(oversample(transmittedSymbols, OVERSAMPLING_FACTOR));

            AWGNCannel.set_noise(get_average_power(transmittedSymbols) / SNRs[i]);
            receivedSignals = AWGNCannel(transmittedSignals);

            receivedSymbols = downsample(fft_normalized(receivedSignals), OVERSAMPLING_FACTOR);

            receivedBits = TSModulator.demodulate_bits(receivedSymbols, AWGNCannel.get_noise());

            BERCounter.count(transmittedBits, receivedBits);
        }

        cout << "\n";
        file << EbN0s[i] << " " << BERCounter.get_errorrate() << endl;

        if (BERCounter.get_errorrate() < BER_LIMIT) break;
    }
}

void calculate_turbo_coded_BER_over_AWGN()
{
    const int constellationSize = pow2i(NUMBER_OF_BITS + 2);
    const int tailSize = 12;
    const int informationLength = NUMBER_OF_BITS * floor_i(NUMBER_OF_SUBCARRIERS - tailSize / static_cast<double>(NUMBER_OF_BITS + 1));
    const int codeLength = informationLength * (NUMBER_OF_BITS + 1) / NUMBER_OF_BITS + tailSize;
    bvec transmittedBits;
    bvec transmittedCodes;
    cvec transmittedSymbols;
    cvec transmittedSignals;
    cvec receivedSignals;
    cvec receivedSymbols;
    vec receivedLLRs;
    bvec receivedBits;
    TrellisShaping TSModulator(constellationSize, TSMetricType::Clipping, TSMappingType(TS_MAPPING_TYPE), TSTerminationType(TS_TERMINATION_TYPE));
    RandomInterleaver randomInterleaver((NUMBER_OF_BITS + 1) * NUMBER_OF_SUBCARRIERS);
    ivec generator("13 15");
    bmat punctureMatrix(3, 2 * NUMBER_OF_BITS);
    punctureMatrix.zeros();
    punctureMatrix.set_submatrix(0, 0, 0, 2 * NUMBER_OF_BITS - 1, 1);
    punctureMatrix(1, 1) = 1;
    punctureMatrix(2, 1) = 1;
    PuncturedTurboCode turboCode;
    turboCode.set_parameters(generator, generator, 4, SRandomInterleaver::make_interleaver(informationLength), NUMBER_OF_TURBO_ITERATIONS, "LOGMAP");
    turboCode.set_puncture_matrix(punctureMatrix);
    AWGN_Channel AWGNCannel;
    BERC BERCounter;
    BLERC FERCounter(informationLength);
    vec EbN0s = linspace(6, 23, 35);
    vec SNRs = TSModulator.get_rate() * turboCode.get_rate() * log2(constellationSize) * inv_dB(EbN0s);
    ofstream file1(join("./Result/CBMTS/BER_turbo_", constellationSize, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_β-", SCALING_FACTOR, "_Term-", TS_TERMINATION_TYPE, ".dat"), ios::out);
    ofstream file2(join("./Result/CBMTS/FER_turbo_", constellationSize, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_β-", SCALING_FACTOR, "_Term-", TS_TERMINATION_TYPE, ".dat"), ios::out);

    for (int i = 0, size = EbN0s.size(); i < size; i++) {
        BERCounter.clear();
        FERCounter.clear();

        for (int j = 0; j < NUMBER_OF_TRIALS; j++) {
            cerr << "trial: " << setw(5) << j + 1 << " Eb/N0(dB): " << setw(2) << EbN0s[i] << " BER: " << BERCounter.get_errorrate() << "\r";

            transmittedBits = randb(informationLength);

            // Encode by turbo code and fill random bits up to N * (B - 1).
            transmittedCodes = turboCode.encode(transmittedBits);
            transmittedCodes.ins(transmittedCodes.size(), randb((NUMBER_OF_BITS + 1) * NUMBER_OF_SUBCARRIERS - codeLength));

            // Random interleave the bit sequence to avoid burst errors.
            transmittedCodes = randomInterleaver.interleave(transmittedCodes);

            // Proccess trellis shaping to reduce PAPR.
            transmittedCodes = TSModulator.form_syndrome(transmittedCodes);
            TSModulator.make_reference_symbols(transmittedCodes, SCALING_FACTOR, CLIPPING_RATIO);
            transmittedSymbols = TSModulator.trellis_coded_modulate(transmittedCodes);

            transmittedSignals = ifft_normalized(oversample(transmittedSymbols, OVERSAMPLING_FACTOR));

            // Transmit signals over AWGN channel.
            AWGNCannel.set_noise(get_average_power(transmittedSymbols) / SNRs[i]);
            receivedSignals = AWGNCannel(transmittedSignals);

            receivedSymbols = downsample(fft_normalized(receivedSignals), OVERSAMPLING_FACTOR);

            // Demodulate and calculate LLRs with trellis.
            receivedLLRs = TSModulator.demodulate_soft_bits(receivedSymbols, AWGNCannel.get_noise());

            // Random interleave the bit sequence to avoid burst errors.
            receivedLLRs = randomInterleaver.deinterleave(receivedLLRs);

            // Decode information bits by turbo decoder.
            receivedBits = turboCode.decode(receivedLLRs.left(codeLength));

            BERCounter.count(transmittedBits, receivedBits);
            FERCounter.count(transmittedBits, receivedBits);
        }

        cout << "\n";
        file1 << EbN0s[i] << " " << BERCounter.get_errorrate() << endl;
        file2 << dB(SNRs[i]) << " " << FERCounter.get_errorrate() << endl;

        if (FERCounter.get_errorrate() < FER_LIMIT) break;
    }
}



int main(const int argc, const char* argv[])
{
    if (argc != 2) {
        cerr << "No arguments error." << endl;
        exit(1);
    }

    string type(argv[1]);
    cout << "Simulation runs now (type: " << type << ")!" << endl;

    // Set a random seed.
    GlobalRNG_randomize();

    if (type == "normal") {
        calculate_instantaneous_CCDF();
    } if (type == "average") {
        calculate_average_power_reduction_capability();
    } else if (type == "clipping") {
        calculate_clipping_ratio_character();
    } else if (type == "uncoded") {
        calculate_uncoded_BER_over_AWGN();
    } else if (type == "turbo") {
        calculate_turbo_coded_BER_over_AWGN();
    } else {
        cerr << "Invalid arguments error." << endl;
        exit(1);
    }

    cout << "Simulation finished successfully!" << endl;
    return 0;
}
