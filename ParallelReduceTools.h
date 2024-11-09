#ifndef PARALLEL_REDUCE_TOOLS_H
#define PARALLEL_REDUCE_TOOLS_H

#include <iostream>
#include <variant>
#include <array>
#include <memory>

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
    Test 09.11.2024 --> only a few possible hardcoded ArrayReduction types in the factory function
    */

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
    }

    /*template<typename SizeType, typename IndexType>
    void performReduction(Kokkos::View<IndexType*> bins, Kokkos::View<SizeType*> localBinHisto,
                          ReductionVariant<SizeType, IndexType>& to_reduce, SizeType localNum, IndexType binCount) {   

        std::visit([&](auto& arg) {
            using T = std::decay_t<decltype(arg)>;
            std::cout << "Performing reduction with type: " << typeid(T).name() << std::endl;
            std::cout << "Size of the_array: " << arg.n_elements << std::endl;
            executeReduction<SizeType, IndexType>(bins, localBinHisto, arg, localNum, binCount);
        }, to_reduce);
    }*/





    /*
    Define logic for maxArrSize different reducer array types where N \in [1, ..., maxArrSize] 
    */

    // Set max array size as a constexpr
    constexpr int maxArrSize = 128;

    /*// Function to execute the Kokkos reduction, separate from recursive template
    template<typename BunchType, typename SizeType, typename IndexType, IndexType N>
    void executeKokkosReduction(BunchType* bunch_m, ArrayReduction<SizeType, IndexType, N>& to_reduce, Kokkos::View<SizeType*> localBinHisto) {
        Kokkos::View<IndexType*> binIndex = bunch_m->bin.getView();  
        // Perform the Kokkos reduction
        Kokkos::parallel_reduce("initLocalHist", bunch_m->getLocalNum(), 
            KOKKOS_LAMBDA(const SizeType& i, ArrayReduction<SizeType, IndexType, N>& update) {
                IndexType ndx = binIndex(i);  // Determine the bin index for this particle
                update.the_array[ndx]++;      // Increment the corresponding bin count in the reduction array
            }, Kokkos::Sum<ArrayReduction<SizeType, IndexType, N>>(to_reduce)
        );

        // Copy the reduced results to the final histogram
        Kokkos::parallel_for("finalize_histogram", N, 
            KOKKOS_LAMBDA(const SizeType& i) {
                localBinHisto(i) = to_reduce.the_array[i];
            }
        );
    }

    // Recursive template function to select the appropriate ArrayReduction<N> type
    template <typename BunchType, typename SizeType, typename IndexType, IndexType N = 1>
    void selectReductionType(std::shared_ptr<BunchType> bunch_m, IndexType binCount, Kokkos::View<SizeType*> localBinHisto) {
        if constexpr (N > maxArrSize) {
            throw std::out_of_range("binCount is out of range");
        } else if (binCount == N) {
            ArrayReduction<SizeType, IndexType, N> to_reduce;
            executeKokkosReduction(bunch_m, to_reduce, localBinHisto);
        } else {
            selectReductionType<BunchType, SizeType, IndexType, N + 1>(bunch_m, binCount, localBinHisto);
        }
    }*/

    /*
    // Helper to create an integer sequence from 1 to N
    template<typename IndexType, IndexType... I>
    constexpr auto make_integer_sequence_from_1(std::integer_sequence<IndexType, I...>) {
        return std::integer_sequence<IndexType, (I + 1)...>{};
    }

    // Main template that generates the sequence from 1 to N
    template<typename IndexType, IndexType N>
    using make_integer_sequence_1_to_N = decltype(make_integer_sequence_from_1(std::make_integer_sequence<IndexType, N>{}));

    // Helper to convert a tuple of types to a variant
    template<typename... Ts>
    struct tuple_to_variant {
        using type = std::variant<Ts...>;
    };

    template<typename Tuple>
    using tuple_to_variant_t = typename tuple_to_variant<Tuple>::type;

    // Helper alias to create a tuple of ArrayReduction types for each element in the sequence
    template<typename SizeType, typename IndexType, IndexType... BinSizes>
    auto make_array_reduction_tuple(std::integer_sequence<IndexType, BinSizes...>) {
        return std::tuple<ArrayReduction<SizeType, IndexType, BinSizes>...>{};
    }

    // Generate tuple of ArrayReduction types from 1 to maxArrSize
    template<typename SizeType, typename IndexType, IndexType maxArrSize>
    using bin_sizes_tuple = decltype(make_array_reduction_tuple<SizeType>(
        make_integer_sequence_1_to_N<IndexType, maxArrSize>{}));

    // Use helper to generate variant from tuple
    template<typename SizeType, typename IndexType, IndexType maxArrSize>
    using ArrayReductionVariant = tuple_to_variant_t<bin_sizes_tuple<SizeType, IndexType, maxArrSize>>;

    // Function to access the correct ArrayReduction type in variant by binCount
    template<typename SizeType, typename IndexType, IndexType maxArrSize>
    auto get_array_reduction_variant(IndexType binCount) {
        using ReductionTuple = bin_sizes_tuple<SizeType, IndexType, maxArrSize>;
        using VariantType = ArrayReductionVariant<SizeType, IndexType, maxArrSize>;

        if (binCount < 1 || binCount > maxArrSize) {
            throw std::out_of_range("binCount is out of range");
        }

        VariantType result;
        std::apply([&](const auto&... reductions) {
            bool matched = ((binCount == decltype(reductions)::size ? (result = reductions, true) : false) || ...);
            if (!matched) throw std::out_of_range("binCount index was invalid");
        }, ReductionTuple{});

        return result;
    }*/


    /*// Shifted integer sequence generator
    template<typename IndexType, IndexType... I>
    constexpr std::integer_sequence<IndexType, (I + 1)...> make_integer_sequence_from_1(std::integer_sequence<IndexType, I...>) {
        return {};
    }

    // Main template that generates the sequence from 1 to N
    template<typename IndexType, IndexType N>
    using make_integer_sequence_1_to_N = decltype(make_integer_sequence_from_1(std::make_integer_sequence<IndexType, N>{})); 

    // Helper alias to create a tuple of ArrayReduction types for each element in the sequence
    template<typename SizeType, typename IndexType, IndexType... BinSizes>
    auto make_array_reduction_tuple(std::integer_sequence<IndexType, BinSizes...>) {
        return std::tuple<ArrayReduction<SizeType, IndexType, BinSizes>...>{};
    }

    // Main template that uses maxArrSize to generate the tuple of types
    template<typename SizeType, typename IndexType>
    using bin_sizes_tuple = decltype(make_array_reduction_tuple<SizeType>(
        make_integer_sequence_1_to_N<IndexType, maxArrSize<IndexType>>{}
    ));

    // Convert tuple of ArrayReduction types to a variant
    template<typename Tuple, std::size_t... Is>
    auto tuple_to_variant_impl(Tuple, std::index_sequence<Is...>) {
        return std::variant<std::tuple_element_t<Is, Tuple>...>{};
    }

    template<typename Tuple>
    using tuple_to_variant_t = decltype(tuple_to_variant_impl(std::declval<Tuple>(), std::make_index_sequence<std::tuple_size_v<Tuple>>{}));

    // Helper function to access tuple element by runtime index
    template <typename Tuple, typename Func, std::size_t... Is>
    auto tuple_at_runtime_index(Tuple&& tuple, std::size_t index, Func&& func, std::index_sequence<Is...>) {
        bool matched = false;
        (..., (index == Is && (func(std::get<Is>(std::forward<Tuple>(tuple))), matched = true)));
        if (!matched) throw std::out_of_range("binCount is out of the allowed range");
    }

    // Function to get ArrayReduction as variant based on binCount (1-based index)
    template<typename SizeType, typename IndexType>
    auto get_array_reduction_variant(IndexType binCount) {
        constexpr auto reductions = bin_sizes_tuple<SizeType, IndexType>{};
        using VariantType = tuple_to_variant_t<decltype(reductions)>;

        if (binCount < 1 || binCount > maxArrSize<IndexType>) {
            throw std::out_of_range("binCount is out of range");
        }

        // Select the correct type by using the tuple index
        VariantType result;
        tuple_at_runtime_index(reductions, binCount - 1, [&](auto&& element) {
            result = element;
        }, std::make_index_sequence<std::tuple_size_v<decltype(reductions)>>{});

        return result;
    }*/




    /*
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
    }*/

    /*
    template<int... Is>
    struct TypeList {
        using Type = std::variant<MyClass<Is>...>;
    };

    using MyVariant = typename TypeList<std::make_integer_sequence<int, 101>{}>::Type;

    MyVariant createInstance(int binCount) {
        if (binCount < 1 || binCount > 100) {
            throw std::out_of_range("binCount must be between 1 and 100");
        }
        
        return std::visit([](auto& instance) -> MyVariant {
            using T = std::decay_t<decltype(instance)>;
            return T();
        }, MyVariant{ std::in_place_index<binCount> });
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