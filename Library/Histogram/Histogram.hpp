#ifndef Histogram_hpp
#define Histogram_hpp

#include <fstream>
#include <itpp/itbase.h>

class Histogram
{
    private:
        const size_t DATA_SIZE;                                     // Data size of axis vectors.
        double min;                                                 // Minimum value of axis.
        double max;                                                 // Maximum value of axis.
        double stride;                                              // Stride of axis.
        bool isValid = false;                                       // Whether range is set or not.
        int underCount = 0;                                         // Count that are under of range.
        int totalCount = 0;                                         // Total Count of histogram.
        itpp::ivec counts;                                          // Counts of existance of data.

        int get_index(double);                                      // Get index corresponded to the given value.
        void increment(int);                                        // Increment count of index corresponded.
        itpp::vec normalize_in_discrete();                          // Normalize to be that integration is 1.
        itpp::vec normalize_in_continous();                         // Normalize to be that integration is 1.

    public:
        Histogram(int size = 1001);                                 // Constructor.
        Histogram(const itpp::vec&);                                // Constructor.
        Histogram(double, double, int size = 1001);                 // Constructor.
        ~Histogram();                                               // Deconstructor.
        void clear();                                               // Clear counts.
        void set_range(double, double);                             // Set min and max if still not.
        void set_range(const itpp::vec&);                           // Set min and max if still not.
        void read_amplitude(const itpp::cvec&);                     // Make instantaneous amplitude spectrum.
        void read_phase(const itpp::cvec&);                         // Make instantaneous phase spectrum.
        void read_power(const itpp::cvec&);                         // Make instantaneous power spectrum.
        void read_normalized_power(const itpp::cvec&);              // Make normalized instantaneous power spectrum.
        void print_discrete_histogram(std::ofstream&);              // Print histogram on a file.
        void print_continuous_histogram(std::ofstream&);            // Print histogram on a file.
        void print_CCDF(std::ofstream&);                            // Print CCDF on a file.
};

#endif
