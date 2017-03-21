#include <iostream>
#include <iomanip>
#include <fstream>
#include <itpp/itcomm.h>
#include "config.hpp"
#include "./Library/Common/helper.hpp"
#include "./Library/Histogram/Histogram.hpp"
#include "./Library/SelectedMapping/SelectedMapping3.hpp"


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
using itpp::cvec;
using itpp::randb;
using itpp::linspace;
using itpp::pow2i;
using itpp::inv_dB;
using itpp::QAM;
using itpp::AWGN_Channel;
using itpp::BERC;
using helper::join;
using helper::oversample;
using helper::downsample;
using helper::fft_normalized;
using helper::ifft_normalized;


const int CONSTELLATION_SIZE = pow2i(NUMBER_OF_BITS);
const double CODE_RATE = (NUMBER_OF_BITS - 1) / NUMBER_OF_BITS;
const double BITS_PER_SYMBOL = log2(CONSTELLATION_SIZE);


void calculate_instantaneous_CCDF()
{
    bvec transmittedBits;
    cvec transmittedSymbols;
    cvec transmittedSignals;
    SelectedMapping3 SLM3Mapper(NUMBER_OF_CANDIDATES, NUMBER_OF_SUBCARRIERS, OVERSAMPLING_FACTOR);
    QAM QAMModulator(CONSTELLATION_SIZE);
    Histogram PSDHistogram(0, 13);
    ofstream file(join("./Result/SLM3/CCDFNorm_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub_", NUMBER_OF_CANDIDATES, "-Cand", ".dat"), ios::out);

    for (int i = 0; i < NUMBER_OF_TRIALS; i++) {
        cerr << "trial: " << i + 1 << "\r";

        transmittedBits = randb(NUMBER_OF_BITS * NUMBER_OF_SUBCARRIERS);

        transmittedSymbols = QAMModulator.modulate_bits(transmittedBits);

        transmittedSignals = SLM3Mapper.select(transmittedSymbols);

        PSDHistogram.read_normalized_power(transmittedSignals);
    }
    cout << "\n";

    PSDHistogram.print_CCDF(file);
}

void calculate_BER_over_AWGN()
{
    bvec transmittedBits;
    cvec transmittedSymbols;
    cvec transmittedSignals;
    cvec receivedSignals;
    cvec receivedSymbols;
    bvec receivedBits;
    QAM QAMModulator(CONSTELLATION_SIZE);
    AWGN_Channel AWGNCannel;
    BERC BERCounter;
    vec EbN0s = linspace(0, 22, 23);
    vec SNRs = CODE_RATE * BITS_PER_SYMBOL * inv_dB(EbN0s);
    ofstream file(join("./Result/SLM3/BER_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub_", NUMBER_OF_CANDIDATES, "-Cand", ".dat"), ios::out);

    for (size_t i = 0, size = EbN0s.size(); i < size; i++) {
        AWGNCannel.set_noise(1.0 / SNRs(i));
        BERCounter.clear();

        for (int j = 0; j < NUMBER_OF_TRIALS; j++) {
            cerr << "trial: " << setw(5) << j + 1 << " Eb/N0(dB): " << setw(2) << EbN0s[i] << " BER: " << BERCounter.get_errorrate() << "\r";

            transmittedBits = randb(NUMBER_OF_BITS * NUMBER_OF_SUBCARRIERS);

            transmittedSymbols = QAMModulator.modulate_bits(transmittedBits);

            transmittedSymbols = ifft_normalized(oversample(transmittedSymbols, OVERSAMPLING_FACTOR));

            receivedSymbols = AWGNCannel(transmittedSymbols);

            receivedSymbols = downsample(fft_normalized(receivedSymbols), OVERSAMPLING_FACTOR);

            receivedBits = QAMModulator.demodulate_bits(receivedSymbols);

            BERCounter.count(transmittedBits, receivedBits);
        }
        cout << "\n";

        file << EbN0s[i] << " " << BERCounter.get_errorrate() << endl;
        if (BERCounter.get_errorrate() < BER_LIMIT) break;
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
        calculate_BER_over_AWGN();
    } else {
        cerr << "Invalid arguments error." << endl;
        exit(1);
    }

    cout << "Simulation finished successfully!" << endl;
    return 0;
}
