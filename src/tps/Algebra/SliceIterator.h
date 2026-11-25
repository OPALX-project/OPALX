#ifndef TPS_SliceIterator_HH
#define TPS_SliceIterator_HH
#include "Utilities/CLRangeError.h"
// Copyright: see NOTICE and LICENSE

// Template class: SliceIterator

#include <compare>
#include <iterator>

namespace tps {

    // Class SliceIterator.
    /// Iterator for array slice.
    ///  An iterator permitting to iterate on a non-constant array
    ///  with a stride different from 1.

    template <class T>
    class SliceIterator {
    public:
        /// The iterator tag, taken from the standard template library.
        typedef std::random_access_iterator_tag iterator_category;

        /// The value type.
        typedef T value_type;

        /// The pointer difference type.
        typedef std::ptrdiff_t difference_type;

        /// The pointer type.
        typedef T* pointer;

        /// The reference type.
        typedef T& reference;

        /// Constructor.
        ///  Construct from the slice origin, nonzero stride and logical position.
        ///  Begin and end of a range share origin and stride; end has position=size.
        ///  Advancing/comparing positions does not form a pointer outside storage.
        ///  Dereference/subscript require a valid element, as for standard iterators.
        SliceIterator() = default;
        SliceIterator(T* array, std::ptrdiff_t stride, std::ptrdiff_t position = 0);

        SliceIterator(const SliceIterator<T>&);
        SliceIterator<T>& operator=(const SliceIterator<T>&);
        bool operator==(const SliceIterator<T>& rhs) const;
        bool operator!=(const SliceIterator<T>& rhs) const;
        std::strong_ordering operator<=>(const SliceIterator& rhs) const {
            if (cursor != rhs.cursor || stride != rhs.stride)
                throw ::CLRangeError(
                        "SliceIterator::compare", "Iterators belong to different slices.");
            return position <=> rhs.position;
        }
        friend SliceIterator operator+(difference_type n, const SliceIterator& iterator) {
            return iterator + n;
        }

        /// Increment (iterate forward).
        SliceIterator<T>& operator++();

        /// Increment (iterate forward).
        SliceIterator<T> operator++(int);

        /// Decrement (iterate backward).
        SliceIterator<T>& operator--();

        /// Decrement (iterate backward).
        SliceIterator<T> operator--(int);

        /// Increment by multiple stride.
        SliceIterator<T>& operator+=(std::ptrdiff_t);

        /// Decrement by multiple stride.
        SliceIterator<T>& operator-=(std::ptrdiff_t);

        /// Add multiple stride.
        SliceIterator<T> operator+(std::ptrdiff_t) const;

        /// Subtract multiple stride.
        SliceIterator<T> operator-(std::ptrdiff_t) const;

        /// Difference.
        difference_type operator-(const SliceIterator<T>&) const;

        /// Dereference.
        T& operator*() const;

        /// Delegate.
        T& operator[](std::ptrdiff_t) const;

    private:
        // Origin stays within storage; end positions never form a pointer.
        T* cursor               = nullptr;
        std::ptrdiff_t position = 0;
        std::ptrdiff_t stride   = 1;
    };

    // Class ConstSliceIterator.
    /// Iterator for array slice.
    ///  An iterator permitting to iterate on a constant array
    ///  with a stride different from 1.

    template <class T>
    class ConstSliceIterator {
    public:
        /// The iterator tag, taken from the standard template library.
        typedef std::random_access_iterator_tag iterator_category;

        /// The value type.
        typedef T value_type;

        /// The pointer difference type.
        typedef std::ptrdiff_t difference_type;

        /// The pointer type.
        typedef const T* pointer;

        /// The reference type.
        typedef const T& reference;

        /// Constructor.
        ///  Construct from the slice origin, nonzero stride and logical position.
        ///  Begin and end of a range share origin and stride; end has position=size.
        ///  Advancing/comparing positions does not form a pointer outside storage.
        ///  Dereference/subscript require a valid element, as for standard iterators.
        ConstSliceIterator() = default;
        ConstSliceIterator(const T* array, std::ptrdiff_t stride, std::ptrdiff_t position = 0);

        ConstSliceIterator(const ConstSliceIterator<T>& rhs);
        ConstSliceIterator<T>& operator=(const ConstSliceIterator<T>&);
        bool operator==(const ConstSliceIterator<T>& rhs) const;
        bool operator!=(const ConstSliceIterator<T>& rhs) const;
        std::strong_ordering operator<=>(const ConstSliceIterator& rhs) const {
            if (cursor != rhs.cursor || stride != rhs.stride)
                throw ::CLRangeError(
                        "ConstSliceIterator::compare", "Iterators belong to different slices.");
            return position <=> rhs.position;
        }
        friend ConstSliceIterator operator+(difference_type n, const ConstSliceIterator& iterator) {
            return iterator + n;
        }

        /// Increment (iterate forward).
        ConstSliceIterator<T>& operator++();

        /// Increment (iterate forward).
        ConstSliceIterator<T> operator++(int);

        /// Decrement (iterate backward).
        ConstSliceIterator<T>& operator--();

        /// Decrement (iterate backward).
        ConstSliceIterator<T> operator--(int);

        /// Add multiple stride and assign.
        ConstSliceIterator<T>& operator+=(std::ptrdiff_t);

        /// Subtract multiple stride and assign.
        ConstSliceIterator<T>& operator-=(std::ptrdiff_t);

        /// Add multiple stride.
        ConstSliceIterator<T> operator+(std::ptrdiff_t) const;

        /// Subtract multiple stride.
        ConstSliceIterator<T> operator-(std::ptrdiff_t) const;

        /// Get pointer difference.
        difference_type operator-(const ConstSliceIterator<T>&) const;

        /// Dereference.
        const T& operator*() const;

        /// Delegate.
        const T& operator[](std::ptrdiff_t) const;

    private:
        const T* cursor         = nullptr;
        std::ptrdiff_t position = 0;
        std::ptrdiff_t stride   = 1;
    };

    // Template function implementation for SliceIterator.

    template <class T>
    inline SliceIterator<T>::SliceIterator(T* array, std::ptrdiff_t n, std::ptrdiff_t p)
        : cursor(array), position(p), stride(n) {
        if (n == 0) throw ::CLRangeError("SliceIterator", "Zero stride.");
    }

    template <class T>
    inline SliceIterator<T>::SliceIterator(const SliceIterator<T>& rhs)
        : cursor(rhs.cursor), position(rhs.position), stride(rhs.stride) {}

    template <class T>
    inline SliceIterator<T>& SliceIterator<T>::operator=(const SliceIterator<T>& rhs) {
        cursor   = rhs.cursor;
        position = rhs.position;
        stride   = rhs.stride;
        return *this;
    }

    template <class T>
    inline bool SliceIterator<T>::operator==(const SliceIterator<T>& rhs) const {
        return cursor == rhs.cursor && stride == rhs.stride && position == rhs.position;
    }

    template <class T>
    inline bool SliceIterator<T>::operator!=(const SliceIterator<T>& rhs) const {
        return !(*this == rhs);
    }

    template <class T>
    inline SliceIterator<T>& SliceIterator<T>::operator++() {
        ++position;
        return *this;
    }

    template <class T>
    inline SliceIterator<T> SliceIterator<T>::operator++(int) {
        SliceIterator<T> tmp(*this);
        ++position;
        return tmp;
    }

    template <class T>
    inline SliceIterator<T>& SliceIterator<T>::operator--() {
        --position;
        return *this;
    }

    template <class T>
    inline SliceIterator<T> SliceIterator<T>::operator--(int) {
        SliceIterator<T> tmp(*this);
        --position;
        return tmp;
    }

    template <class T>
    inline SliceIterator<T>& SliceIterator<T>::operator+=(std::ptrdiff_t n) {
        position += n;
        return *this;
    }

    template <class T>
    inline SliceIterator<T>& SliceIterator<T>::operator-=(std::ptrdiff_t n) {
        position -= n;
        return *this;
    }

    template <class T>
    inline SliceIterator<T> SliceIterator<T>::operator+(std::ptrdiff_t n) const {
        SliceIterator<T> tmp(*this);
        return tmp += n;
    }

    template <class T>
    inline SliceIterator<T> SliceIterator<T>::operator-(std::ptrdiff_t n) const {
        SliceIterator<T> tmp(*this);
        return tmp -= n;
    }

    template <class T>
    inline typename SliceIterator<T>::difference_type SliceIterator<T>::operator-(
            const SliceIterator<T>& rhs) const {
        if (cursor != rhs.cursor || stride != rhs.stride)
            throw ::CLRangeError(
                    "SliceIterator::difference", "Iterators have different origins or strides.");
        return position - rhs.position;
    }

    template <class T>
    inline T& SliceIterator<T>::operator*() const {
        return cursor[position * stride];
    }

    template <class T>
    inline T& SliceIterator<T>::operator[](std::ptrdiff_t n) const {
        return cursor[(position + n) * stride];
    }

    // Template function implementation for ConstSliceIterator.

    template <class T>
    inline ConstSliceIterator<T>::ConstSliceIterator(
            const T* array, std::ptrdiff_t n, std::ptrdiff_t p)
        : cursor(array), position(p), stride(n) {
        if (n == 0) throw ::CLRangeError("ConstSliceIterator", "Zero stride.");
    }

    template <class T>
    inline ConstSliceIterator<T>::ConstSliceIterator(const ConstSliceIterator<T>& rhs)
        : cursor(rhs.cursor), position(rhs.position), stride(rhs.stride) {}

    template <class T>
    inline ConstSliceIterator<T>& ConstSliceIterator<T>::operator=(
            const ConstSliceIterator<T>& rhs) {
        cursor   = rhs.cursor;
        position = rhs.position;
        stride   = rhs.stride;
        return *this;
    }

    template <class T>
    inline bool ConstSliceIterator<T>::operator==(const ConstSliceIterator<T>& rhs) const {
        return cursor == rhs.cursor && stride == rhs.stride && position == rhs.position;
    }

    template <class T>
    inline bool ConstSliceIterator<T>::operator!=(const ConstSliceIterator<T>& rhs) const {
        return !(*this == rhs);
    }

    template <class T>
    inline ConstSliceIterator<T>& ConstSliceIterator<T>::operator++() {
        ++position;
        return *this;
    }

    template <class T>
    inline ConstSliceIterator<T> ConstSliceIterator<T>::operator++(int) {
        ConstSliceIterator<T> tmp(*this);
        ++position;
        return tmp;
    }

    template <class T>
    inline ConstSliceIterator<T>& ConstSliceIterator<T>::operator--() {
        --position;
        return *this;
    }

    template <class T>
    inline ConstSliceIterator<T> ConstSliceIterator<T>::operator--(int) {
        ConstSliceIterator<T> tmp(*this);
        --position;
        return tmp;
    }

    template <class T>
    inline ConstSliceIterator<T>& ConstSliceIterator<T>::operator+=(std::ptrdiff_t n) {
        position += n;
        return *this;
    }

    template <class T>
    inline ConstSliceIterator<T>& ConstSliceIterator<T>::operator-=(std::ptrdiff_t n) {
        position -= n;
        return *this;
    }

    template <class T>
    inline ConstSliceIterator<T> ConstSliceIterator<T>::operator+(std::ptrdiff_t n) const {
        ConstSliceIterator<T> tmp(*this);
        return tmp += n;
    }

    template <class T>
    inline ConstSliceIterator<T> ConstSliceIterator<T>::operator-(std::ptrdiff_t n) const {
        ConstSliceIterator<T> tmp(*this);
        return tmp -= n;
    }

    template <class T>
    inline typename ConstSliceIterator<T>::difference_type ConstSliceIterator<T>::operator-(
            const ConstSliceIterator<T>& rhs) const {
        if (cursor != rhs.cursor || stride != rhs.stride)
            throw ::CLRangeError(
                    "SliceIterator::difference", "Iterators have different origins or strides.");
        return position - rhs.position;
    }

    template <class T>
    inline const T& ConstSliceIterator<T>::operator*() const {
        return cursor[position * stride];
    }

    template <class T>
    inline const T& ConstSliceIterator<T>::operator[](std::ptrdiff_t n) const {
        return cursor[(position + n) * stride];
    }

}  // namespace tps

#endif  // TPS_SliceIterator_HH
