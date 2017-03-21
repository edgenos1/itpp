#include <iostream>
#include <iomanip>
#include <fstream>
#include <itpp/itcomm.h>
#include "config.hpp"
#include "./Library/Common/helper.hpp"
#include "./Library/Histogram/Histogram.hpp"
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
using itpp::to_bvec;
using itpp::QAM;
using itpp::LDPC_Parity_Regular;
using itpp::LDPC_Generator_Systematic;
using itpp::LDPC_Code;
using itpp::Punctured_Turbo_Codec;
using itpp::AWGN_Channel;
using itpp::BERC;
using itpp::BLERC;
using helper::join;
using helper::oversample;
using helper::downsample;
using helper::fft_normalized;
using helper::ifft_normalized;
using helper::hard_decide_LLRs;


void calculate_instantaneous_CCDF()
{
    const int constellationSize = pow2i(NUMBER_OF_BITS);
    bvec transmittedBits;
    cvec transmittedSymbols;
    cvec transmittedSignals;
    Histogram PSDHistogram(0, 13, 1001);
    QAM QAMModulator(constellationSize);
    ofstream file(join("./Result/OFDM/CCDF_Norm_", constellationSize, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);

    for (int i = 0; i < NUMBER_OF_TRIALS; i++) {
        cout << "trial: " << i + 1 << "\r";

        transmittedBits = randb(NUMBER_OF_BITS * NUMBER_OF_SUBCARRIERS);

        transmittedSymbols = QAMModulator.modulate_bits(transmittedBits);

        transmittedSignals = ifft_normalized(oversample(transmittedSymbols, OVERSAMPLING_FACTOR));

        PSDHistogram.read_normalized_power(transmittedSignals);
    }
    cout << "\n";

    PSDHistogram.print_CCDF(file);
}

void calculate_uncoded_BER_over_AWGN()
{
    const int constellationSize = pow2i(NUMBER_OF_BITS);
    bvec transmittedBits;
    cvec transmittedSymbols;
    cvec transmittedSignals;
    cvec receivedSignals;
    cvec receivedSymbols;
    bvec receivedBits;
    QAM QAMModulator(constellationSize);
    AWGN_Channel AWGNCannel;
    BERC BERCounter;
    vec EbN0s = linspace(-3, 22, 23);
    vec SNRs = log2(constellationSize) * inv_dB(EbN0s);
    ofstream file(join("./Result/OFDM/BER_Uncod_", constellationSize, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);

    for (size_t i = 0, size = EbN0s.size(); i < size; i++) {
        AWGNCannel.set_noise(pow(SNRs[i], -1.0));
        BERCounter.clear();

        for (int j = 0; j < NUMBER_OF_TRIALS; j++) {
            cerr << "trial: " << setw(5) << j + 1 << " Eb/N0(dB): " << setw(2) << EbN0s[i] << " BER: " << BERCounter.get_errorrate() << "\r";

            transmittedBits = randb(NUMBER_OF_BITS * NUMBER_OF_SUBCARRIERS);

            transmittedSymbols = QAMModulator.modulate_bits(transmittedBits);

            transmittedSignals = ifft_normalized(oversample(transmittedSymbols, OVERSAMPLING_FACTOR));

            receivedSignals = AWGNCannel(transmittedSignals);

            receivedSymbols = downsample(fft_normalized(receivedSignals), OVERSAMPLING_FACTOR);

            receivedBits = hard_decide_LLRs(QAMModulator.demodulate_soft_bits(receivedSymbols, AWGNCannel.get_noise()));

            BERCounter.count(transmittedBits, receivedBits);
        }

        cout << "\n";
        file << EbN0s[i] << " " << BERCounter.get_errorrate() << endl;

        if (BERCounter.get_errorrate() < BER_LIMIT) break;
    }
}

void calculate_LDPC_coded_BER_over_AWGN()
{
    const int constellationSize = pow2i(NUMBER_OF_BITS);
    bvec transmittedBits;
    bvec transmittedCodes;
    cvec transmittedSymbols;
    cvec transmittedSignals;
    cvec receivedSignals;
    cvec receivedSymbols;
    bvec receivedCodes;
    bvec receivedBits;
    QAM QAMModulator(constellationSize);
    LDPC_Parity_Regular LDPCParity;
    LDPC_Generator_Systematic LDPCGenerator;
    LDPC_Code LDPCCode;
    AWGN_Channel AWGNCannel;
    BERC BERCounter;
    vec EbN0s = linspace(0, 22, 23);
    ofstream file(join("./Result/OFDM/BER_LDPC_", constellationSize, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);

    LDPCParity.generate(768, 3, 6);
    LDPCGenerator.construct(&LDPCParity);
    LDPCCode.set_code(&LDPCParity, &LDPCGenerator);
    LDPCCode.set_exit_conditions(50);
    vec SNRs = LDPCCode.get_rate() * log2(constellationSize) * inv_dB(EbN0s);

    for (size_t i = 0, size = EbN0s.size(); i < size; i++) {
        AWGNCannel.set_noise(pow(SNRs[i], -1.0) / 2.0);
        BERCounter.clear();

        for (int j = 0; j < NUMBER_OF_TRIALS; j++) {
            cerr << "trial: " << setw(5) << j + 1 << " Eb/N0(dB): " << setw(2) << EbN0s[i] << " BER: " << BERCounter.get_errorrate() << "\r";

            transmittedBits = randb(LDPCCode.get_rate() * NUMBER_OF_BITS * NUMBER_OF_SUBCARRIERS);

            transmittedCodes = LDPCCode.encode(transmittedBits);

            transmittedSymbols = QAMModulator.modulate_bits(transmittedCodes);

            transmittedSignals = ifft_normalized(oversample(transmittedSymbols, OVERSAMPLING_FACTOR));

            receivedSignals = AWGNCannel(transmittedSignals);

            receivedSymbols = downsample(fft_normalized(receivedSignals), OVERSAMPLING_FACTOR);

            receivedBits = LDPCCode.decode(QAMModulator.demodulate_soft_bits(receivedSymbols, AWGNCannel.get_noise()));

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
    AWGN_Channel AWGNCannel;
    BERC BERCounter;
    BLERC FERCounter(informationLength);
    vec EbN0s = linspace(6, 22, 17);
    vec SNRs = turboCode.get_rate() * log2(constellationSize) * inv_dB(EbN0s);
    ofstream file1(join("./Result/OFDM/BER_Turbo_", constellationSize, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);
    ofstream file2(join("./Result/OFDM/FER_Turbo_", constellationSize, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);


    for (int i = 0, size = EbN0s.size(); i < size; i++) {
        AWGNCannel.set_noise(pow(SNRs[i], -1.0));
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
            transmittedSignals = ifft_normalized(oversample(transmittedSymbols, OVERSAMPLING_FACTOR));

            // Transmit signals over AWGN channel.
            receivedSignals = AWGNCannel(transmittedSignals);

            receivedSymbols = downsample(fft_normalized(receivedSignals), OVERSAMPLING_FACTOR);
            receivedLLRs = QAMModulator.demodulate_soft_bits(receivedSymbols, AWGNCannel.get_noise());

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
    } else if (type == "ldpc") {
        calculate_LDPC_coded_BER_over_AWGN();
    } else if (type == "turbo") {
        calculate_turbo_coded_BER_over_AWGN();
    } else {
        cerr << "Invalid arguments error." << endl;
        exit(1);
    }

    cout << "Simulation finished successfully!" << endl;
    return 0;
}
