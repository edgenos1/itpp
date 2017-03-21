#include <iostream>
#include <iomanip>
#include <fstream>
#include <itpp/itcomm.h>
#include "config.hpp"
#include "./Library/Common/helper.hpp"
#include "./Library/TurboCode/RecursiveConvolutionalCode.hpp"


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

    bvec bits;

    RecursiveConvolutionalCode convolutionalCode;
    convolutionalCode.set_generator_polynomials(ivec("15 11"), 4);
    convolutionalCode.encode(randb(16));

    return 0;
}
