#ifndef PARALLEL_REDUCE_TOOLS_H
#define PARALLEL_REDUCE_TOOLS_H

#include <variant> // for std::variant
// #include <memory>
#include <utility> // for std::index_sequence

namespace ParticleBinning {

    // Possibility to manually choose between different reduction methods! 
    enum class HistoReductionMode {
        Standard,          // Default auto selection
        ParallelReduce,    // Force usage of parallel_reduce if binCount <= maxArrSize
        TeamBased,         // Force team-based/atomic reduction if gpu enabled
        HostOnly
    };

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

    // Set max array size as a constexpr
    template<typename IndexType>
    constexpr IndexType maxArrSize = 5; // 128 needs a few minutes to compile. Good in between magic number is 30. Fast compilation with 15

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
    

    /*
    The following struct is only used for the host version of the code. It is not used in the CUDA version.
    ...since dynamically allocated arrays work only on host (but still faster than team based!)
     */
    template<typename SizeType, typename IndexType>
    struct HostArrayReduction {
        SizeType* the_array;

        // Static variable to define array size
        static IndexType binCountStatic;

        // Only compile without CUDA, since it is not needed otherwise and would not compile with the reduction identity below!
        #ifndef KOKKOS_ENABLE_CUDA
        HostArrayReduction() { 
            the_array = new SizeType[binCountStatic];
            for (IndexType i = 0; i < binCountStatic; i++ ) { the_array[i] = 0; }
        }

        HostArrayReduction(const HostArrayReduction& rhs) { 
            the_array = new SizeType[binCountStatic];
            for (IndexType i = 0; i < binCountStatic; i++ ){ the_array[i] = rhs.the_array[i]; }
        }
        
        ~HostArrayReduction() { delete[] the_array; }
        
        HostArrayReduction& operator=(const HostArrayReduction& rhs) {
            the_array = new SizeType[binCountStatic];
            if (this != &rhs) {
                for (IndexType i = 0; i < binCountStatic; ++i) { the_array[i] = rhs.the_array[i]; }
            }
            return *this;
        }
        
        HostArrayReduction& operator+=(const HostArrayReduction& src) {
            for (IndexType i = 0; i < binCountStatic; i++ ) { the_array[i] += src.the_array[i]; }
            return *this;
        }
        #else
        KOKKOS_INLINE_FUNCTION
        HostArrayReduction& operator+=(const HostArrayReduction& src) {
            // throw an error if this function is called on CUDA (shouldn't happen)
            Kokkos::abort("Error: HostArrayReduction is not supported on CUDA!\n       Note: It exists only for compilation compatibility.");
            return *this;
        }
        #endif
    };

    // Initialize the static variable outside the struct
    template<typename SizeType, typename IndexType>
    IndexType HostArrayReduction<SizeType, IndexType>::binCountStatic = 10;


    /**
     * The following are some helper functions for debugging
     */

    /**
    * @brief Computes the \(p\)-norm of a vector field (for debugging purpose).
    * 
    * @param field The input vector field.
    * @param p The order of the norm (default is 2 for Euclidean norm).
    * 
    * @return The computed \(p\)-norm of the vector field.
    */
    template<typename T, unsigned Dim>
    T vnorm(const VField_t<T, Dim>& field, int p = 2) {
        T sum = 0;
        ippl::parallel_reduce("VectorFieldNormReduce", field.getFieldRangePolicy(),
            KOKKOS_LAMBDA(const ippl::RangePolicy<Dim>::index_array_type& idx, T& loc_sum) {
                ippl::Vector<T, Dim> e = apply(field, idx);
                loc_sum += std::pow(e.dot(e), p/2.0);
            }, Kokkos::Sum<T>(sum));
        return std::pow(sum, 1.0/p);
    }

}

namespace Kokkos {  
    // This one is for usage on GPU
    template<typename SizeType, typename IndexType, IndexType N>
    struct reduction_identity<ParticleBinning::ArrayReduction<SizeType, IndexType, N>> {
        KOKKOS_FORCEINLINE_FUNCTION static ParticleBinning::ArrayReduction<SizeType, IndexType, N> sum() {
            return ParticleBinning::ArrayReduction<SizeType, IndexType, N>();
        }
    };
    
    // This one is for usage on host
    template<typename SizeType, typename IndexType>
    struct reduction_identity<ParticleBinning::HostArrayReduction<SizeType, IndexType>> {
        KOKKOS_FORCEINLINE_FUNCTION static ParticleBinning::HostArrayReduction<SizeType, IndexType> sum() {
            return ParticleBinning::HostArrayReduction<SizeType, IndexType>();
        }
    };
}

#endif // PARALLEL_REDUCE_TOOLS_H