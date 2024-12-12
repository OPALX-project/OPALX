#ifndef BINHISTO_HPP
#define BINHISTO_HPP

namespace ParticleBinning {


    template <typename size_type, typename bin_index_type, typename value_type, bool UseDualView, class... Properties>
    void Histogram<size_type, bin_index_type, value_type, UseDualView, Properties...>::copyFields(const Histogram& other) {
        debug_name_m    = other.debug_name_m;
        numBins_m       = other.numBins_m;
        totalBinWidth_m = other.totalBinWidth_m;

        histogram_m = other.histogram_m;
        binWidths_m = other.binWidths_m;
        postSum_m   = other.postSum_m;
    }

    template <typename size_type, typename bin_index_type, typename value_type, bool UseDualView, class... Properties>
    void Histogram<size_type, bin_index_type, value_type, UseDualView, Properties...>::mergeBins(value_type maxBinRatio) {
        // TODO 
        // Should merge neighbouring bins such that the width/N_part ratio is roughly maxBinRatio.
        // TODO: Find algorithm for that
        std::cout << "Warning: mergeBins not implemented yet!" << std::endl;
    }


} // namespace ParticleBinning

#endif // BINHISTO_HPP