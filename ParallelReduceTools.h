#ifndef PARALLEL_REDUCE_TOOLS_H
#define PARALLEL_REDUCE_TOOLS_H

#include <variant> // for std::variant
#include <memory>
#include <utility> // for std::index_sequence

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
    Test 09.11.2024 --> only a few possible hardcoded ArrayReduction types in the factory function (for testing)
    

    // Define a variant type that can hold any of the Reduction<N> types
    template<typename SizeType, typename IndexType>
    using ReductionVariant = std::variant<
        ArrayReduction<SizeType, IndexType, 9>, ArrayReduction<SizeType, IndexType, 10>, ArrayReduction<SizeType, IndexType, 11>, ArrayReduction<SizeType, IndexType, 12>
    >;

    // Factory function to create the appropriate reduction object
    template<typename SizeType, typename IndexType>
    ReductionVariant<SizeType, IndexType> createReductionObject(int binCount) {
        switch (binCount) {
            case 9: return ArrayReduction<SizeType, IndexType, 9>();
            case 10: return ArrayReduction<SizeType, IndexType, 10>();
            case 11: return ArrayReduction<SizeType, IndexType, 11>();
            case 12: return ArrayReduction<SizeType, IndexType, 12>();
            default: throw std::out_of_range("binCount is out of the allowed range");
        }
    }*/

    /*
    Define logic for maxArrSize different reducer array types where N \in [1, ..., maxArrSize] 
    */

    // Set max array size as a constexpr
    //template<typename IndexType>
    //constexpr IndexType maxArrSize = 50;
    template<typename IndexType>
    constexpr IndexType maxArrSize = 15; // 128 needs a few minutes to compile. Good in between is 30. Fast compilation with 15

    // Primary template for ReductionVariantHelper (not used directly)
    template<typename SizeType, typename IndexType, typename Sequence>
    struct ReductionVariantHelper;

    // Specialization of ReductionVariantHelper that accepts std::integer_sequence and expands it
    template<typename SizeType, typename IndexType, IndexType... Sizes>
    struct ReductionVariantHelper<SizeType, IndexType, std::integer_sequence<IndexType, Sizes...>> {
        using type = std::variant<ArrayReduction<SizeType, IndexType, Sizes + 1>...>;
    };

    // Define the ReductionVariant type alias using the helper with std::make_integer_sequence
    template<typename SizeType, typename IndexType>
    using ReductionVariant = typename ReductionVariantHelper<SizeType, IndexType, std::make_integer_sequence<IndexType, maxArrSize<IndexType>>>::type;

    template<typename SizeType, typename IndexType, IndexType N>
    ReductionVariant<SizeType, IndexType> createReductionObjectHelper(IndexType binCount) {
        if constexpr (N > maxArrSize<IndexType>) {
            throw std::out_of_range("binCount is out of the allowed range");
        } else if (binCount == N) {
            return ArrayReduction<SizeType, IndexType, N>();
        } else {
            return createReductionObjectHelper<SizeType, IndexType, N + 1>(binCount);
        }
    }

    template<typename SizeType, typename IndexType>
    ReductionVariant<SizeType, IndexType> createReductionObject(IndexType binCount) {
        return createReductionObjectHelper<SizeType, IndexType, 1>(binCount);
    }



    // Factory function to create the appropriate reduction object
    /*template<typename SizeType, typename IndexType>
    ReductionVariant<SizeType, IndexType> createReductionObject(IndexType binCount) {
        switch (binCount) {
            case 9: return ArrayReduction<SizeType, IndexType, 9>();
            case 10: return ArrayReduction<SizeType, IndexType, 10>();
            case 11: return ArrayReduction<SizeType, IndexType, 11>();
            case 12: return ArrayReduction<SizeType, IndexType, 12>();
            default: throw std::out_of_range("binCount is out of the allowed range");
        }
    }*/



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