#ifndef TrellisShaping_hpp
#define TrellisShaping_hpp

#include <itpp/itcomm.h>
#include "../../Library/TrellisShaping/TrellisCode.hpp"
#include "../../Library/TrellisShaping/QAM.hpp"

enum class TSMetricType {Autocorrelation, Clipping};
enum class TSShapingMethod {Sign, Tail};
enum class TSTerminationType {Unterminated, Tailbite};

class TrellisShaping
{
    private:
        TSMetricType metric;
        TSMappingType mappingType;
        TSShapingMethod shapingMethod;
        TSTerminationType terminationType;
        int numberOfShapingBits;
        double rate;
        bool isValid = false;
        bool isSetUp = false;
        TrellisCode generator;
        TrellisCode inverseSymdrome;
        TrellisCode parityChecker;
        TrellisCode compoundDecoder;
        QAM QAMModulator;
        itpp::cvec referenceSymbols;

        double get_optimum_clipping_ratio(double);
        void set_constellation(const int, const TSMappingType mappingType);
        void set_shaping_method(const TSMetricType);
        void set_termination_type(const TSTerminationType);
        void set_convolutional_code();
        itpp::bvec add_shaping_bits(const itpp::bvec&, const itpp::bvec&);
        itpp::cvec encode_unterminated(const itpp::bvec&);
        itpp::cvec encode_tailbite(const itpp::bvec&);
        double calculate_metric(const itpp::cvec&, const itpp::cvec&, int);
        double calculate_transit_probability(const int, const int, const std::complex<double>&, const double);
        itpp::vec decode_by_SOVA(const itpp::cvec&, const double);
        itpp::vec decode_by_BCJR(const itpp::cvec&, const double);

    public:
        TrellisShaping(const int m, const TSMetricType, const TSMappingType, const TSTerminationType = TSTerminationType::Unterminated);
        ~TrellisShaping();
        double get_rate();
        void make_reference_symbols(const itpp::bvec&, double, double = -1);
        itpp::bvec form_syndrome(const itpp::bvec&);
        itpp::cvec trellis_coded_modulate(const itpp::bvec&);
        itpp::bvec demodulate_bits(const itpp::cvec&, const double);
        itpp::vec demodulate_soft_bits(const itpp::cvec&, const double);
};

#endif
