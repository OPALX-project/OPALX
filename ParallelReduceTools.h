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


    /**
     * @brief A templated structure for performing array-based reductions in parallel computations.
     * 
     * @tparam SizeType The type used for the elements of the array.
     * @tparam IndexType The type used for indexing the array.
     * @tparam N The fixed size of the array.
     * 
     * This structure provides functionality to initialize, copy, assign, and perform element-wise
     * addition on arrays of fixed size. It is designed to be used in parallel reduction operations
     * with Kokkos.
     * 
     * Member Functions:
     * - **ArrayReduction()**: Default constructor that initializes all elements of the array to zero.
     * - **ArrayReduction(const ArrayReduction& rhs)**: Copy constructor that performs a deep copy of the array.
     * - **ArrayReduction& operator=(const ArrayReduction& rhs)**: Assignment operator that performs a deep copy of the array.
     * - **ArrayReduction& operator+=(const ArrayReduction& src)**: Element-wise addition operator that adds the elements
     *   of another `ArrayReduction` instance to the current instance.
     */
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

    /**
     * @brief Maximum allowed array size for compile-time array reduction types.
     * 
     * This constexpr variable defines the upper limit for statically-sized array reductions
     * that are used in parallel histogram reductions.
     * For small histograms, this is prefereable to pure atomics. However, for large values,
     * compile time is increased significantly. 
     * 
     * @tparam IndexType The type used for array indexing and size specification.
     * @note 5 is empirically chosen, since in real world applications the team-based reduction
     *       is usually enough. 
     */
    template<typename IndexType>
    constexpr IndexType maxArrSize = 5;

    /**
     * @brief Primary template for generating std::variant types containing ArrayReduction specializations.
     * 
     * This template is not used directly but serves as the base for the specialized version.
     * It works in conjunction with `std::integer_sequence` to generate a variant containing
     * ArrayReduction types with sizes from `1` to `maxArrSize`.
     * 
     * This is necessary to pre-compile types for histogram reduction on device. 
     * 
     * @tparam SizeType The type used for array elements in the reduction operations.
     * @tparam IndexType The type used for array indexing.
     * @tparam Sequence A std::integer_sequence type that will be expanded in the specialization.
     */
    template<typename SizeType, typename IndexType, typename Sequence>
    struct ReductionVariantHelper;

    /**
     * @brief Specialized template that expands `std::integer_sequence` to create a variant of `ArrayReduction` types.
     * 
     * This specialization takes a `std::integer_sequence` and expands it into a `std::variant`
     * containing `ArrayReduction` types with consecutive sizes. Each size in the sequence
     * is incremented by 1 to create `ArrayReduction<SizeType, IndexType, Sizes + 1>`.
     * 
     * @tparam SizeType The type used for array elements in the reduction operations.
     * @tparam IndexType The type used for array indexing.
     * @tparam Sizes... Parameter pack of integer values from the `std::integer_sequence`.
     */
    template<typename SizeType, typename IndexType, IndexType... Sizes>
    struct ReductionVariantHelper<SizeType, IndexType, std::integer_sequence<IndexType, Sizes...>> {
        using type = std::variant<ArrayReduction<SizeType, IndexType, Sizes + 1>...>;
    };

    /**
     * @brief Type alias for a `std::variant` containing all possible `ArrayReduction` types up to `maxArrSize`.
     * 
     * This alias creates a variant that can hold `ArrayReduction` objects with sizes from `1` to `maxArrSize`.
     * It uses `ReductionVariantHelper` with `std::make_integer_sequence` to generate all variants
     * type at compile time.
     * 
     * @tparam SizeType The type used for array elements in the reduction operations.
     * @tparam IndexType The type used for array indexing and size specification.
     */
    template<typename SizeType, typename IndexType>
    using ReductionVariant = typename ReductionVariantHelper<SizeType, IndexType, std::make_integer_sequence<IndexType, maxArrSize<IndexType>>>::type;

    /**
     * @brief Recursive helper function to create `ArrayReduction` objects with compile-time size matching.
     * 
     * This function uses compile-time recursion to find and instantiate the correct `ArrayReduction`
     * type based on the requested `binCount`. It recursively increments the template parameter `N`
     * until it matches `binCount`, then returns the appropriate `ArrayReduction` instance.
     * 
     * @tparam SizeType The type used for array elements in the reduction operations.
     * @tparam IndexType The type used for array indexing.
     * @tparam N Current array size being tested in the recursive search.
     * @param binCount The desired number of bins for the reduction array.
     * 
     * @return ReductionVariant containing the `ArrayReduction` object with size matching `binCount`.
     * 
     * @throws `std::out_of_range` If `binCount` exceeds `maxArrSize`.
     * 
     * @note This function should not be called directly. Use `createReductionObject()` instead and use 
     *       `std::visit` to access the reducer.
     */
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

    /**
     * @brief Factory function to create `ArrayReduction` objects with runtime-specified size.
     * 
     * This function provides an interface for creating `ArrayReduction` objects where the 
     * array size is determined at runtime. It internally uses compile-time template recursion 
     * to match the runtime `binCount` with the template specialization. 
     * 
     * This allows using "runtime sized arrays" as reducer objects inside a Kokkos kernel
     * on device.
     * 
     * @tparam SizeType The type used for array elements in the reduction operations.
     * @tparam IndexType The type used for array indexing.
     * @param binCount The desired number of bins for the reduction array (must be $\leq$ `maxArrSize`).
     * 
     * @return ReductionVariant containing the `ArrayReduction` object with the specified size.
     * 
     * @throws `std::out_of_range` If `binCount` exceeds `maxArrSize` or is less than 1.
     * 
     * Example usage (also see implementation of `AdaptBins::initLocalHisto()`):
     * ```
     * std::visit(
     *     [&](auto& reducer_arr) {
     *         executeInitLocalHistoReduction(reducer_arr);
     *     }, createReductionObject<size_type, bin_index_type>(binCount)
     * );
     * ```
     */
    template<typename SizeType, typename IndexType>
    ReductionVariant<SizeType, IndexType> createReductionObject(IndexType binCount) {
        return createReductionObjectHelper<SizeType, IndexType, 1>(binCount);
    }
    

    /**
     * @brief Host-only array reduction structure for dynamic array sizes in `Kokkos::parallel_reduce`.
     * 
     * This structure provides array-based reduction functionality specifically only for host execution.
     * Unlike `ArrayReduction` which uses compile-time fixed sizes, HostArrayReduction uses dynamically
     * allocated arrays with runtime-determined sizes. This approach offers better performance than
     * team-based reductions on host systems (low concurrency) while maintaining flexibility for 
     * arbitrary bin counts (way more efficient than pure atomics).
     * 
     * @tparam SizeType The type used for the elements of the array.
     * @tparam IndexType The type used for indexing and array size specification.
     * 
     * @note This structure is conditionally compiled and only available when KOKKOS_ENABLE_CUDA is not defined.
     * @note All instances of this class share the same array size through the static binCountStatic member.
     * 
     * Usage Example (`ReducerType` is a specialization of `HostArrayReduction`):
     * ```
     * Kokkos::parallel_reduce("initLocalHist", bunch_m->getLocalNum(), 
     *     KOKKOS_LAMBDA(const size_type& i, ReducerType& update) {
     *           bin_index_type ndx = binIndex(i);  // Determine the bin index for this particle
     *           update.the_array[ndx]++;           // Increment the corresponding bin count in the reduction array
     *     }, Kokkos::Sum<ReducerType>(to_reduce)
     * );
     * ```
     */
    template<typename SizeType, typename IndexType>
    struct HostArrayReduction {

        /**
         * @brief Pointer to the dynamically allocated array for reduction operations.
         */
        SizeType* the_array;

        /**
         * @brief Static variable defining the array size for all instances of this template specialization.
         * 
         * @note This value is set outside when the classed is used, typically before any instances are created.
         * @warning Changing this value after object creation may lead to inconsistent behavior.
         */
        static IndexType binCountStatic;

#ifndef KOKKOS_ENABLE_CUDA
        /**
         * @brief Default constructor that allocates and zero-initializes the array.
         */
        HostArrayReduction() { 
            the_array = new SizeType[binCountStatic];
            for (IndexType i = 0; i < binCountStatic; i++ ) { the_array[i] = 0; }
        }

        /**
         * @brief Copy constructor that performs a deep copy of the array from another instance.
         * 
         * @param rhs The instance to copy from.
         */
        HostArrayReduction(const HostArrayReduction& rhs) { 
            the_array = new SizeType[binCountStatic];
            for (IndexType i = 0; i < binCountStatic; i++ ){ the_array[i] = rhs.the_array[i]; }
        }
        
        /**
         * @brief Destructor that deallocates the dynamically allocated array.
         */
        ~HostArrayReduction() { delete[] the_array; }
        
        /**
         * @brief Assignment operator that performs a deep copy of the array from another instance.
         * 
         * @param rhs The instance to copy from.
         * @return Reference to this instance after assignment.
         */
        HostArrayReduction& operator=(const HostArrayReduction& rhs) {
            the_array = new SizeType[binCountStatic];
            if (this != &rhs) {
                for (IndexType i = 0; i < binCountStatic; ++i) { the_array[i] = rhs.the_array[i]; }
            }
            return *this;
        }
        
        /**
         * @brief Element-wise addition operator that adds the elements of another `HostArrayReduction` instance.
         * 
         * @param src The source instance to add to this instance.
         * @return Reference to this instance after addition.
         */
        HostArrayReduction& operator+=(const HostArrayReduction& src) {
            for (IndexType i = 0; i < binCountStatic; i++ ) { the_array[i] += src.the_array[i]; }
            return *this;
        }
#else
        /**
         * @brief CUDA-disabled version of the addition operator that throws an error.
         * 
         * @note This constructor is not intended to be used on CUDA, as `HostArrayReduction` is not supported there.
         *       Instead, it throws an error to indicate that this functionality is not available.
         */
        KOKKOS_INLINE_FUNCTION
        HostArrayReduction& operator+=(const HostArrayReduction& src) {
            Kokkos::abort("Error: HostArrayReduction is not supported on CUDA!\n       Note: It exists only for compilation compatibility.");
            return *this;
        }
#endif
    };

    /**
     * @brief Static member initialization for the array size variable.
     * 
     * This line provides the definition and default initialization for the static member
     * `binCountStatic`. The default value of 10 serves as a place holder and should
     * be overridden based on the required histogram size.
     * 
     * @tparam SizeType The type used for array elements.
     * @tparam IndexType The type used for indexing and size specification.
     * 
     * @note This must be defined outside the class template for proper static member initialization.
     * @note The value should be set appropriately before creating any `HostArrayReduction` instances.
     * @note This is only used in the host version of the code. 
     */
    template<typename SizeType, typename IndexType>
    IndexType HostArrayReduction<SizeType, IndexType>::binCountStatic = 10;


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
    /**
     * @brief Kokkos reduction identity specialization for custom reduction types.
     * 
     * These specializations of `Kokkos::reduction_identity` are required to enable the use of
     * custom `ArrayReduction` types with Kokkos built-in reducers. The reduction identity 
     * function provides the neutral element for the reduction operation, which serves as the 
     * initial value for thread-private reduction variables.
     * 
     * @note The GPU version (`ArrayReduction`) uses compile-time fixed arrays, while the host 
     *       version (`HostArrayReduction`) uses dynamically allocated arrays for flexibility 
     *       with runtime-determined sizes.
     */
    template<typename SizeType, typename IndexType, IndexType N>
    struct reduction_identity<ParticleBinning::ArrayReduction<SizeType, IndexType, N>> {
        KOKKOS_FORCEINLINE_FUNCTION static ParticleBinning::ArrayReduction<SizeType, IndexType, N> sum() {
            return ParticleBinning::ArrayReduction<SizeType, IndexType, N>();
        }
    };
    
    /**
     * @brief Kokkos reduction identity specialization for host-only `HostArrayReduction` types.
     * @see reduction_identity<ParticleBinning::ArrayReduction> for detailed explanation of reduction identities.
     */
    template<typename SizeType, typename IndexType>
    struct reduction_identity<ParticleBinning::HostArrayReduction<SizeType, IndexType>> {
        KOKKOS_FORCEINLINE_FUNCTION static ParticleBinning::HostArrayReduction<SizeType, IndexType> sum() {
            return ParticleBinning::HostArrayReduction<SizeType, IndexType>();
        }
    };
}

#endif // PARALLEL_REDUCE_TOOLS_H