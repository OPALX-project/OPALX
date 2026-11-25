#ifndef TPS_Array1D_HH
#define TPS_Array1D_HH
#include "Utilities/CLRangeError.h"
// Copyright: see NOTICE and LICENSE

// Template class: Array1D

#include <algorithm>
#include <ostream>
#include <utility>
#include "tps/HostStorage.hpp"

namespace tps {

    // Template class Array1D<T>
    /// One-dimensional array.
    ///  A templated representation for one-dimensional arrays.
    ///  This class implements storage management and component access,
    ///  but contains no arithmetic operations.

    template <class T>
    class Array1D {
    public:
        /// The value type of this array.
        typedef T value_type;

        /// The iterator type for the array.
        typedef T* iterator;

        /// The iterator type for constant array.
        typedef const T* const_iterator;

        /// Default constructor.
        ///  Constructs an array of zero length.
        Array1D();

        /// Constructor.
        ///  Reserves space for `n` elements and value-initializes them.
        explicit Array1D(int n);

        /// Constructor.
        ///  Reserves space for `n` elements and store `n` copies of
        ///  `t`.
        Array1D(int n, const T& t);

        Array1D(const Array1D<T>&);
        ~Array1D() = default;
        Array1D<T>& operator=(const Array1D<T>&);

        /// Get reference to element.
        ///  Return a reference to the element in position `n`.
        ///  Throw ::CLRangeError, if `n` is out of range.
        T& operator()(int n);

        /// Get value of element.
        ///  Return the value of the element in position `n`.
        ///  Throw ::CLRangeError, if `n` is out of range.
        const T& operator()(int n) const;

        /// Get reference to element.
        ///  Return a reference to the element in position `n`.
        ///  Throws ::CLRangeError, if `n` is out of range.
        T& operator[](int);

        /// Get value of element.
        ///  Return the value of the element in position `n`.
        ///  Throws ::CLRangeError, if `n` is out of range.
        const T& operator[](int) const;

        /// Get beginning of data.
        ///  Return pointer to beginning of array.
        ///  Version for non-constant objects.
        iterator begin();

        /// Get end of data.
        ///  Return pointer past end of array.
        ///  Version for non-constant objects.
        iterator end();

        /// Get beginning of data.
        ///  Return pointer to beginning of array.
        ///  Version for constant objects.
        const_iterator begin() const;

        /// Get end of data.
        ///  Return pointer past end of array.
        ///  Version for constant objects.
        const_iterator end() const;

        /// Get array size.
        int size() const;

        /// Change array size.
        ///  Elements newly allocated are value-initialized.
        void resize(int size);

        /** @brief Borrow mutable HostSpace storage without transferring ownership.
         * @details The View is unmanaged; finish all uses before resize, assignment
         * or destruction of the owning array. Copies of the array own separate storage.
         */
        auto view() {
            return Kokkos::View<T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                    data, len);
        }
        /// Borrow read-only storage, with the same lifetime constraints as view().
        auto view() const {
            return Kokkos::View<
                    const T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                    data, len);
        }

    protected:
        static int checkedSize(int n) {
            if (n < 0) throw ::CLRangeError("Array1D::size", "Negative array size.");
            return n;
        }

        // The array size.
        int len;

        // The array data.
        using Storage = detail::HostStorage<T>;
        Storage storage;
        T* data;  // non-owning view, synchronized with storage
    };

    // Template function implementation.

    template <class T>
    Array1D<T>::Array1D() : len(0), data(0) {}

    template <class T>
    Array1D<T>::Array1D(const Array1D<T>& array)
        : len(array.len),
          storage(len ? Storage("tps_Array1D", len) : Storage()),
          data(storage.data()) {
        if (len) std::copy(array.begin(), array.end(), begin());
    }

    template <class T>
    Array1D<T>::Array1D(int n)
        : len(checkedSize(n)),
          storage(len ? Storage("tps_Array1D", len) : Storage()),
          data(storage.data()) {}

    template <class T>
    Array1D<T>::Array1D(int n, const T& val)
        : len(checkedSize(n)),
          storage(len ? Storage("tps_Array1D", len) : Storage()),
          data(storage.data()) {
        std::fill(begin(), end(), val);
    }

    template <class T>
    Array1D<T>& Array1D<T>::operator=(const Array1D<T>& rhs) {
        if (this != &rhs) {
            Array1D copy(rhs);
            std::swap(len, copy.len);
            std::swap(storage, copy.storage);
            std::swap(data, copy.data);
        }
        return *this;
    }

    template <class T>
    inline T& Array1D<T>::operator[](int i) {
        return (*this)(i);
    }

    template <class T>
    inline const T& Array1D<T>::operator[](int i) const {
        return (*this)(i);
    }

    template <class T>
    T& Array1D<T>::operator()(int i) {
        if (i < 0 || i >= size()) {
            throw ::CLRangeError("Array1D::operator()", "Index out of range.");
        }
        return data[i];
    }

    template <class T>
    const T& Array1D<T>::operator()(int i) const {
        if (i < 0 || i >= size()) {
            throw ::CLRangeError("Array1D::operator()", "Index out of range.");
        }
        return data[i];
    }

    template <class T>
    typename Array1D<T>::iterator Array1D<T>::begin() {
        return data;
    }

    template <class T>
    typename Array1D<T>::iterator Array1D<T>::end() {
        return len ? data + len : data;
    }

    template <class T>
    typename Array1D<T>::const_iterator Array1D<T>::begin() const {
        return data;
    }

    template <class T>
    typename Array1D<T>::const_iterator Array1D<T>::end() const {
        return len ? data + len : data;
    }

    template <class T>
    int Array1D<T>::size() const {
        return len;
    }

    template <class T>
    void Array1D<T>::resize(int n) {
        checkedSize(n);
        if (len < n) {
            auto replacement = Storage("tps_Array1D_resize", n);
            if (len) std::copy(data, data + len, replacement.data());
            storage = std::move(replacement);
            data    = storage.data();
        }
        len = n;
    }

    template <class T>
    std::ostream& operator<<(std::ostream& os, const Array1D<T>& v) {
        for (int i = 0; i < v.size(); ++i) {
            os << v[i] << " ";
        }

        os << std::endl;
        return os;
    }

}  // namespace tps

#endif  // TPS_Array1D_HH
