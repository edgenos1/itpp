#include <iostream>
#include <iomanip>
#include <fstream>
#include <itpp/itcomm.h>
#include "config.hpp"
#include "./Library/Common/helper.hpp"
#include "./Library/Histogram/Histogram.hpp"
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
using itpp::inv_dB;
using itpp::BPSK_c;
using itpp::AWGN_Channel;
using itpp::BERC;
using helper::join;


int main(const int argc, const char* argv[])
{
    // Set a random seed.
    GlobalRNG_randomize();

    const int informationLength = 10000;
    const int iterations = 15;
    bvec transmittedBits;
    bvec transmittedCodes;
    cvec transmittedSymbols;
    cvec receivedSymbols;
    vec receivedLLRs;
    bvec receivedBits;
    BPSK_c BPSKModulator;
    ivec generator("13 15");
    bmat punctureMatrix(
        "1, 1, 1, 1, 1, 1, 1, 1, 1, 1;"
        "0, 1, 0, 0, 0, 0, 0, 0, 0, 0;"
        "0, 1, 0, 0, 0, 0, 0, 0, 0, 0;"
    );
    PuncturedTurboCode turboCode;
    turboCode.set_parameters(generator, generator, 4, SRandomInterleaver::make_interleaver(informationLength), iterations, "LOGMAP");
    turboCode.set_puncture_matrix(punctureMatrix);
    AWGN_Channel AWGNCannel;
    BERC BERCounter;
    vec EbN0s = linspace(3.0, 3.5, 6);
    vec SNRs = inv_dB(EbN0s) * turboCode.get_rate();
    ofstream file(join("./Result/BER_Turbo_BPSK_", iterations, "-Itr", ".dat"), ios::out);

    for (size_t i = 0, size = EbN0s.size(); i < size; i++) {
        AWGNCannel.set_noise(pow(SNRs[i], -1.0));
        BERCounter.clear();

        for (int j = 0; j < NUMBER_OF_TRIALS; j++) {
            cerr << "trial: " << setw(5) << j + 1 << " Eb/N0(dB): " << setw(3) << EbN0s[i] << " BER: " << BERCounter.get_errorrate() << "\r";

            transmittedBits = randb(informationLength);

            // Encode by turbo code.
            transmittedCodes = turboCode.encode(transmittedBits);

            transmittedSymbols = BPSKModulator.modulate_bits(transmittedCodes);

            // Transmit signals over AWGN channel.
            receivedSymbols = AWGNCannel(transmittedSymbols);

            receivedLLRs = BPSKModulator.demodulate_soft_bits(receivedSymbols, AWGNCannel.get_noise());

            // Decode by turbo code.
            receivedBits = turboCode.decode(receivedLLRs);

            BERCounter.count(transmittedBits, receivedBits);
        }

        cout << "\n";
        file << EbN0s[i] << " " << BERCounter.get_errorrate() << endl;

        if (BERCounter.get_errorrate() < BER_LIMIT) break;
    }

    return 0;
}
