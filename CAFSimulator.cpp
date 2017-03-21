#include <iostream>
#include <iomanip>
#include <fstream>
#include <itpp/itcomm.h>
#include "config.hpp"
#include "./Library/Common/helper.hpp"
#include "./Library/Histogram/Histogram.hpp"
#include "./Library/ClippingAndFiltering/ClippingAndFiltering.hpp"
#include "./Library/Interleaver/RandomInterleaver.hpp"
#include "./Library/Interleaver/SRandomInterleaver.hpp"
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
using itpp::ivec;
using itpp::bvec;
using itpp::cvec;
using itpp::bmat;
using itpp::randb;
using itpp::linspace;
using itpp::pow2i;
using itpp::floor_i;
using itpp::dB;
using itpp::inv_dB;
using itpp::QAM;
using itpp::AWGN_Channel;
using itpp::BERC;
using itpp::BLERC;
using helper::join;
using helper::oversample;
using helper::downsample;
using helper::fft_normalized;
using helper::ifft_normalized;
using helper::get_average_power;
using helper::hard_decide_LLRs;


void calculate_instantaneous_CCDF()
{
    const int constellationSize = pow2i(NUMBER_OF_BITS);
    const int informationLength = NUMBER_OF_BITS * NUMBER_OF_SUBCARRIERS;
    bvec transmittedBits;
    cvec transmittedSymbols;
    cvec transmittedSignals;
    Histogram PSDHistogram(0, 13, 1001);
    QAM QAMModulator(constellationSize);
    ClippingAndFiltering CAFProccesser(CLIPPING_RATIO, OVERSAMPLING_FACTOR);
    ofstream file(join("./Result/CAF/CCDF_Norm_", constellationSize, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_γ-", CLIPPING_RATIO, ".dat"), ios::out);

    for (int i = 0; i < NUMBER_OF_TRIALS; i++) {
        cout << "trial: " << i + 1 << "\r";

        transmittedBits = randb(informationLength);

        transmittedSymbols = QAMModulator.modulate_bits(transmittedBits);

        for (int j = 0; j < NUMBER_OF_CLIPPING_ITERATIONS; j++) {
            transmittedSymbols = CAFProccesser.clip_and_filter(transmittedSymbols);
        }

        transmittedSignals = ifft_normalized(oversample(transmittedSymbols, OVERSAMPLING_FACTOR));

        PSDHistogram.read_normalized_power(transmittedSignals);
    }

    cout << "\n";
    PSDHistogram.print_CCDF(file);
}

void calculate_uncoded_BER_over_AWGN()
{
    const int constellationSize = pow2i(NUMBER_OF_BITS);
    const int informationLength = NUMBER_OF_BITS * NUMBER_OF_SUBCARRIERS;
    bvec transmittedBits;
    cvec transmittedSymbols;
    cvec transmittedSignals;
    cvec receivedSignals;
    cvec receivedSymbols;
    vec receivedLLRs;
    bvec receivedBits;
    cvec clippingAttenuation;
    QAM QAMModulator(constellationSize);
    ClippingAndFiltering CAFProccesser(CLIPPING_RATIO, OVERSAMPLING_FACTOR);
    clippingAttenuation.set_size(informationLength);
    clippingAttenuation.set_subvector(0, clippingAttenuation.size() - 1, CAFProccesser.calculate_attenuation_factor());
    AWGN_Channel AWGNCannel;
    BERC BERCounter;
    vec EbN0s = linspace(0, 22, 23);
    vec SNRs = log2(constellationSize) * inv_dB(EbN0s);
    ofstream file(join("./Result/CAF/BER_Uncod_", constellationSize, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_γ-", CLIPPING_RATIO, ".dat"), ios::out);

    for (size_t i = 0, size = EbN0s.size(); i < size; i++) {
        BERCounter.clear();

        for (int j = 0; j < NUMBER_OF_TRIALS; j++) {
            cerr << "trial: " << setw(5) << j + 1 << " Eb/N0(dB): " << setw(2) << EbN0s[i] << " BER: " << BERCounter.get_errorrate() << "\r";

            transmittedBits = randb(informationLength);

            transmittedSymbols = QAMModulator.modulate_bits(transmittedBits);

            for (int k = 0; k < NUMBER_OF_CLIPPING_ITERATIONS; k++) {
                transmittedSymbols = CAFProccesser.clip_and_filter(transmittedSymbols);
            }

            transmittedSignals = ifft_normalized(oversample(transmittedSymbols, OVERSAMPLING_FACTOR));

            AWGNCannel.set_noise(get_average_power(transmittedSymbols) / SNRs[i]);
            receivedSignals = AWGNCannel(transmittedSignals);

            receivedSymbols = downsample(fft_normalized(receivedSignals), OVERSAMPLING_FACTOR);

            receivedLLRs = QAMModulator.demodulate_soft_bits(receivedSymbols, clippingAttenuation, AWGNCannel.get_noise());

            receivedBits = hard_decide_LLRs(receivedLLRs);

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
    const int informationLength = NUMBER_OF_BITS * floor_i(NUMBER_OF_SUBCARRIERS - tailSize / static_cast<double>(NUMBER_OF_BITS + 2));
    const int codeLength = informationLength * (NUMBER_OF_BITS + 2) / NUMBER_OF_BITS + tailSize;
    bvec transmittedBits;
    bvec transmittedCodes;
    cvec transmittedSymbols;
    cvec transmittedSignals;
    cvec receivedSignals;
    cvec receivedSymbols;
    vec receivedLLRs;
    bvec receivedBits;
    cvec clippingAttenuation;
    QAM QAMModulator(constellationSize);
    RandomInterleaver interleaver((NUMBER_OF_BITS + 2) * NUMBER_OF_SUBCARRIERS);
    ivec generator("13 15");
    bmat punctureMatrix(3, NUMBER_OF_BITS);
    punctureMatrix.zeros();
    punctureMatrix.set_submatrix(0, 0, 0, NUMBER_OF_BITS - 1, 1);
    punctureMatrix(1, 1) = 1;
    punctureMatrix(2, NUMBER_OF_BITS * 3 / 2 - 6) = 1;
    PuncturedTurboCode turboCode;
    turboCode.set_parameters(generator, generator, 4, SRandomInterleaver::make_interleaver(informationLength), NUMBER_OF_TURBO_ITERATIONS, "LOGMAP");
    turboCode.set_puncture_matrix(punctureMatrix);
    ClippingAndFiltering CAFProccesser(CLIPPING_RATIO, OVERSAMPLING_FACTOR);
    clippingAttenuation.set_size((NUMBER_OF_BITS + 2) * NUMBER_OF_SUBCARRIERS);
    clippingAttenuation.set_subvector(0, clippingAttenuation.size() - 1, CAFProccesser.calculate_attenuation_factor());
    AWGN_Channel AWGNCannel;
    BERC BERCounter;
    BLERC FERCounter(informationLength);
    vec EbN0s = linspace(10, 20, 21);
    vec SNRs = turboCode.get_rate() * log2(constellationSize) * inv_dB(EbN0s);
    ofstream file1(join("./Result/CAF/BER_Turbo_", constellationSize, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_γ-", CLIPPING_RATIO, ".dat"), ios::out);
    ofstream file2(join("./Result/CAF/FER_Turbo_", constellationSize, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_γ-", CLIPPING_RATIO, ".dat"), ios::out);


    for (int i = 0, size = EbN0s.size(); i < size; i++) {
        BERCounter.clear();
        FERCounter.clear();

        for (int j = 0; j < NUMBER_OF_TRIALS; j++) {
            cerr << "trial: " << setw(5) << j + 1 << " Eb/N0(dB): " << setw(2) << EbN0s[i] << " BER: " << BERCounter.get_errorrate() << "\r";

            transmittedBits = randb(informationLength);

            // Encode by turbo code and fill random bits up to N * (B - 1).
            transmittedCodes = turboCode.encode(transmittedBits);
            transmittedCodes.ins(transmittedCodes.size(), randb((NUMBER_OF_BITS + 2) * NUMBER_OF_SUBCARRIERS - codeLength));

            // Random interleave the bit sequence to avoid burst errors.
            transmittedCodes = interleaver.interleave(transmittedCodes);

            transmittedSymbols = QAMModulator.modulate_bits(transmittedCodes);

            for (int k = 0; k < NUMBER_OF_CLIPPING_ITERATIONS; k++) {
                transmittedSymbols = CAFProccesser.clip_and_filter(transmittedSymbols);
            }

            transmittedSignals = ifft_normalized(oversample(transmittedSymbols, OVERSAMPLING_FACTOR));

            // Transmit signals over AWGN channel.
            AWGNCannel.set_noise(get_average_power(transmittedSymbols) / SNRs[i]);
            receivedSignals = AWGNCannel(transmittedSignals);

            receivedSymbols = downsample(fft_normalized(receivedSignals), OVERSAMPLING_FACTOR);
            receivedLLRs = QAMModulator.demodulate_soft_bits(receivedSymbols, clippingAttenuation, AWGNCannel.get_noise());

            // Deinterleave the bit sequence to avoid burst errors.
            receivedLLRs = interleaver.deinterleave(receivedLLRs);

            // Demodulate and calculate LLRs with trellis.
            receivedBits = turboCode.decode(receivedLLRs.left(codeLength));

            BERCounter.count(transmittedBits, receivedBits);
            FERCounter.count(transmittedBits, receivedBits);
        }

        cout << "\n";
        file1 << EbN0s[i] << " " << BERCounter.get_errorrate() << endl;
        file2 << dB(SNRs[i]) << " " << FERCounter.get_errorrate() << endl;

        if (BERCounter.get_errorrate() < BER_LIMIT || FERCounter.get_errorrate() < FER_LIMIT) break;
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
