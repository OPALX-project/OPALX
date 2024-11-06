#ifndef PARALLEL_REDUCE_TOOLS_H
#define PARALLEL_REDUCE_TOOLS_H

#include <array>
#include <variant>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace ParticleBinning {  
    template<typename SizeType, typename IndexType, IndexType N>
    struct ArrayReduction {
        SizeType the_array[N];

        KOKKOS_INLINE_FUNCTION 
        ArrayReduction() { 
            for (IndexType i = 0; i < N; i++ ) { the_array[i] = 0; }
        }
        KOKKOS_INLINE_FUNCTION  
        ArrayReduction(const ArrayReduction& rhs) { 
            for (IndexType i = 0; i < N; i++ ){ the_array[i] = rhs.the_array[i]; }
        }
        KOKKOS_INLINE_FUNCTION
        ArrayReduction& operator=(const ArrayReduction& rhs) {
            if (this != &rhs) {
                for (IndexType i = 0; i < N; ++i) { the_array[i] = rhs.the_array[i]; }
            }
            return *this;
        }
        KOKKOS_INLINE_FUNCTION
        ArrayReduction& operator+=(const ArrayReduction& src) {
            for (IndexType i = 0; i < N; i++ ) { the_array[i] += src.the_array[i]; }
            return *this;
        }
    };

    /*
    Define logic for maxArrSize different reducer array types where N \in [1, ..., maxArrSize] 
    */
    // Define maxArrSize as a constexpr value based on IndexType
    template<typename IndexType>
    constexpr IndexType maxArrSize = 128;

    // Generate an array of values from 1 to maxArrSize for the given IndexType
    template<typename IndexType, IndexType... I>
    constexpr std::array<IndexType, sizeof...(I)> generate_sequence(std::integer_sequence<IndexType, I...>) {
        return { static_cast<IndexType>(I + 1)... };  // Generates [1, 2, ..., maxArrSize]
    }

    // Generate the bin_sizes array using maxArrSize and IndexType
    template<typename IndexType>
    constexpr auto bin_sizes = generate_sequence<IndexType>(std::make_integer_sequence<IndexType, maxArrSize<IndexType>>{});



    // Helper to create a variant of ArrayReduction types based on bin_sizes
    template<typename SizeType, typename IndexType, IndexType... I>
    using ReductionVariantHelper = std::variant<ArrayReduction<SizeType, IndexType, bin_sizes<IndexType>[I]>...>;

    // Define ReductionVariant by directly expanding `std::make_integer_sequence`
    template<typename SizeType, typename IndexType>
    using ReductionVariant = ReductionVariantHelper<SizeType, IndexType, std::make_integer_sequence<IndexType, maxArrSize<IndexType>>{}>;

    // Function to create the correct ArrayReduction type based on binCount
    template<typename SizeType, typename IndexType>
    ReductionVariant<SizeType, IndexType> selectReductionType(IndexType binCount) {
        for (IndexType i = 0; i < maxArrSize<IndexType>; ++i) {
            if (bin_sizes<IndexType>[i] == binCount) {
                // Return a properly constructed ArrayReduction type in a variant
                return ReductionVariant<SizeType, IndexType>(ArrayReduction<SizeType, IndexType, bin_sizes<IndexType>[i]>());
            }
        }
        throw std::runtime_error("Unsupported bin count");
    }




}

namespace Kokkos {  
    template<typename SizeType, typename IndexType, IndexType N>
    struct reduction_identity<ParticleBinning::ArrayReduction<SizeType, IndexType, N>> {
        KOKKOS_FORCEINLINE_FUNCTION static ParticleBinning::ArrayReduction<SizeType, IndexType, N> sum() {
            return ParticleBinning::ArrayReduction<SizeType, IndexType, N>();
        }
    };
}

#endif // PARALLEL_REDUCE_TOOLS_H