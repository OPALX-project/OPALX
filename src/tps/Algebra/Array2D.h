#ifndef TPS_Array2D_HH
#define TPS_Array2D_HH
#include "Utilities/CLRangeError.h"
#include "Utilities/LogicalError.h"
// Copyright: see NOTICE and LICENSE

// Template class: Array2D

#include "tps/Algebra/Array1D.h"
#include "tps/Algebra/SliceIterator.h"
#include <algorithm>
#include "tps/HostStorage.hpp"
#include <utility>
#include <climits>
#include <iosfwd>

namespace tps {

// Template class Array2D<T>
/// Two-dimensional array.
///  A templated representation for 2-dimensional arrays.
///  This class implements storage management and component access,
///  as well as access to rows and columns, but contains no arithmetic
///  operations.

template <class T> class Array2D {

  public:
    /// The value type of the array.
    typedef T value_type;

    /// The iterator type for sequential access to all elements.
    typedef T* iterator;

    /// The iterator type for sequential access to all elements
    ///  of a constant array.
    typedef const T* const_iterator;

    /// The iterator type for access by rows.
    typedef T* row_iterator;

    /// The iterator type for access by rows for a constant array.
    typedef const T* const_row_iterator;

    /// The iterator type for access by columns.
    typedef SliceIterator<T> col_iterator;

    /// The iterator type for access by columns for a constant array.
    typedef ConstSliceIterator<T> const_col_iterator;

    // Checked row views keep both subscripts of a[r][c] bounds checked.
    class row_view {
      public:
        row_view(Array2D<T>& owner, int row) : owner(owner), row(row) {}
        T& operator[](int column) const {
            return owner(row, column);
        }

      private:
        Array2D<T>& owner;
        int row;
    };
    class const_row_view {
      public:
        const_row_view(const Array2D<T>& owner, int row) : owner(owner), row(row) {}
        const T& operator[](int column) const {
            return owner(row, column);
        }

      private:
        const Array2D<T>& owner;
        int row;
    };

    /// Default constructor.
    ///  Constructs array of zero rows and zero columns.
    Array2D();

    /// Constructor.
    ///  Reserves `rows` times `cols` elements and leaves
    ///  them value-initialized.
    Array2D(int rows, int cols);

    /// Constructor.
    ///  Reserves `rows` x `cols` elements and sets them to `t`.
    Array2D(int rows, int cols, const T& t);

    Array2D(const Array2D<T>&);
    ~Array2D() = default;
    Array2D<T>& operator=(const Array2D<T>&);

    /// Get element reference.
    ///  Return a reference to element in row `r` and column `c`.
    T& operator()(int r, int c);

    /// Get element value.
    ///  Return the value of element in row `r` and column `c`.
    const T& operator()(int r, int c) const;

    /// Get number of rows.
    int nrows() const;

    /// Get number of columns.
    int ncols() const;

    /// Get total size (rows times columns).
    int size() const;

    /// Get pointer to beginning of data.
    ///  Treat the array as one-dimensional.
    ///  Version for non-constant objects.
    iterator begin();

    /// Get pointer past end of data.
    ///  Treat the array as one-dimensional.
    ///  Version for non-constant objects.
    iterator end();

    /// Get pointer to beginning of data.
    ///  Treat the array as one-dimensional.
    ///  Version for constant objects.
    const_iterator begin() const;

    /// Get pointer past end of data.
    ///  Treat the array as one-dimensional.
    ///  Version for constant objects.
    const_iterator end() const;

    /// Get row iterator.
    ///  Return pointer to beginning of row r.
    ///  Throw ::CLRangeError, if row is out of range.
    ///  Version for non-constant objects.
    row_iterator row_begin(int r);

    /// Get row iterator.
    ///  Return pointer past end of row r.
    ///  Throw ::CLRangeError, if row is out of range.
    ///  Version for non-constant objects.
    row_iterator row_end(int r);

    /// Get row iterator.
    ///  Return pointer to beginning of constant row r.
    ///  Throws ::CLRangeError, if row is out of range.
    ///  Version for constant objects.
    const_row_iterator row_begin(int r) const;

    /// Get row iterator.
    ///  Return pointer past end of constant row r.
    ///  Throws ::CLRangeError, if row is out of range.
    ///  Version for constant objects.
    const_row_iterator row_end(int r) const;

    /// Get checked row view.
    ///  Return a view supporting checked column subscripting.
    ///  Throws ::CLRangeError, if either subscript is out of range.
    ///  Version for non-constant objects.
    row_view operator[](int r);

    /// Get checked row view.
    ///  Return a view supporting checked column subscripting.
    ///  Throws ::CLRangeError, if either subscript is out of range.
    ///  Version for constant objects.
    const_row_view operator[](int r) const;

    /// Get column iterator.
    ///  Return pointer to beginning of column c.
    ///  Throw ::CLRangeError if c is out of range.
    ///  Version for non-constant objects.
    col_iterator col_begin(int c);

    /// Get column iterator.
    ///  Return pointer past end of column c.
    ///  Throw ::CLRangeError if c is out of range.
    ///  Version for non-constant objects.
    col_iterator col_end(int c);

    /// Get column iterator.
    ///  Return pointer to beginning of column c.
    ///  Throw ::CLRangeError if c is out of range.
    ///  Version for constant objects.
    const_col_iterator col_begin(int c) const;

    /// Get column iterator.
    ///  Return pointer past end of column c.
    ///  Throw ::CLRangeError if c is out of range.
    ///  Version for constant objects.
    const_col_iterator col_end(int c) const;

    /// Fetch column.
    ///  Store the column `c` into the array `toArray`.
    ///  Throw ::CLRangeError if c is out of range.
    void getColumn(Array1D<T>& toArray, int c) const;

    /// Fetch row.
    ///  Store the row `r` into the array `toArray`.
    ///  Throw ::CLRangeError if r is out of range.
    void getRow(Array1D<T>& toArray, int r) const;

    /// Store column.
    ///  Fill the column `c` from the array `fromArray`.
    ///  Throw ::CLRangeError if c is out of range.
    void putColumn(const Array1D<T>& fromArray, int c);

    /// Store row.
    ///  Fill the row `r` from the array `fromArray`.
    ///  Throw ::CLRangeError if r is out of range.
    void putRow(const Array1D<T>& fromArray, int r);

    /// Exchange columns.
    ///  Exchange the columns `c1` and `c2`.

    //  Throw ::CLRangeError, if either index is out of range.
    void swapColumns(int c1, int c2);

    /// Exchange rows.
    ///  Exchange the rows `r1` and `r2`.
    ///  Throw ::CLRangeError, if either index is out of range.
    void swapRows(int row1, int row2);

    /** @brief Borrow mutable HostSpace storage without transferring ownership.
     * @details The View is unmanaged; finish all uses before resize, assignment
     * or destruction of the owning array. Copies of the array own separate storage.
     */
    auto view() {
        return Kokkos::View<T**, Kokkos::LayoutRight, Kokkos::HostSpace,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>(data, rows, cols);
    }
    /// Borrow read-only storage, with the same lifetime constraints as view().
    auto view() const {
        return Kokkos::View<const T**, Kokkos::LayoutRight, Kokkos::HostSpace,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>(data, rows, cols);
    }

    /// Borrow all entries in row-major order; the owning array must remain unchanged/live.
    auto flatView() {
        return Kokkos::View<T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(data,
                                                                                            len);
    }
    /// Borrow read-only entries in row-major order; no ownership is transferred.
    auto flatView() const {
        return Kokkos::View<const T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            data, len);
    }

  protected:
    static int checkedSize(int r, int c) {
        if (r < 0 || c < 0 || (c && r > INT_MAX / c))
            throw ::CLRangeError("Array2D::size", "Negative or overflowing dimensions.");
        return r * c;
    }

    // Array dimensions.
    int rows;
    int cols;
    int len; // rows * cols

    // Array data.
    using Storage = detail::HostStorage<T>;
    Storage storage;
    T* data; // non-owning view, synchronized with storage
};

// Template implementation.

template <class T> inline Array2D<T>::Array2D() : rows(0), cols(0), len(0), data(0) {}

template <class T>
inline Array2D<T>::Array2D(const Array2D<T>& array)
    : rows(array.rows), cols(array.cols), len(array.len),
      storage(len ? Storage("tps_Array2D", len) : Storage()), data(storage.data()) {
    if (len)
        std::copy(array.begin(), array.end(), begin());
}

template <class T>
inline Array2D<T>::Array2D(int r, int c)
    : rows(r), cols(c), len(checkedSize(r, c)),
      storage(len ? Storage("tps_Array2D", len) : Storage()), data(storage.data()) {}

template <class T>
inline Array2D<T>::Array2D(int r, int c, const T& val)
    : rows(r), cols(c), len(checkedSize(r, c)),
      storage(len ? Storage("tps_Array2D", len) : Storage()), data(storage.data()) {
    std::fill(begin(), end(), val);
}

template <class T> inline Array2D<T>& Array2D<T>::operator=(const Array2D<T>& rhs) {
    if (this != &rhs) {
        Array2D copy(rhs);
        std::swap(rows, copy.rows);
        std::swap(cols, copy.cols);
        std::swap(len, copy.len);
        std::swap(storage, copy.storage);
        std::swap(data, copy.data);
    }
    return *this;
}

template <class T> inline T& Array2D<T>::operator()(int r, int c) {
    if (r < 0 || r >= rows || c < 0 || c >= ncols()) {
        throw ::CLRangeError("Array2d::operator()", "Index out of range.");
    }
    return data[cols * r + c];
}

template <class T> inline const T& Array2D<T>::operator()(int r, int c) const {
    if (r < 0 || r >= rows || c < 0 || c >= ncols()) {
        throw ::CLRangeError("Array2d::operator()", "Index out of range.");
    }
    return data[cols * r + c];
}

template <class T> inline int Array2D<T>::nrows() const {
    return rows;
}

template <class T> inline int Array2D<T>::ncols() const {
    return cols;
}

template <class T> inline int Array2D<T>::size() const {
    return len;
}

template <class T> inline typename Array2D<T>::iterator Array2D<T>::begin() {
    return data;
}

template <class T> inline typename Array2D<T>::iterator Array2D<T>::end() {
    return len ? data + len : data;
}

template <class T> inline typename Array2D<T>::const_iterator Array2D<T>::begin() const {
    return data;
}

template <class T> inline typename Array2D<T>::const_iterator Array2D<T>::end() const {
    return len ? data + len : data;
}

template <class T> typename Array2D<T>::row_iterator Array2D<T>::row_begin(int r) {
    if (r < 0 || r >= rows) {
        throw ::CLRangeError("Array2D::row_begin()", "Row index out of range.");
    }
    return cols ? data + cols * r : data;
}

template <class T> typename Array2D<T>::row_iterator Array2D<T>::row_end(int r) {
    if (r < 0 || r >= rows) {
        throw ::CLRangeError("Array2D::row_end()", "Row index out of range.");
    }
    return cols ? data + cols * (r + 1) : data;
}

template <class T> typename Array2D<T>::const_row_iterator Array2D<T>::row_begin(int r) const {
    if (r < 0 || r >= rows) {
        throw ::CLRangeError("Array2D::row_begin()", "Row index out of range.");
    }
    return cols ? data + cols * r : data;
}

template <class T> typename Array2D<T>::const_row_iterator Array2D<T>::row_end(int r) const {
    if (r < 0 || r >= rows) {
        throw ::CLRangeError("Array2D::row_end()", "Row index out of range.");
    }
    return cols ? data + cols * (r + 1) : data;
}

template <class T> typename Array2D<T>::row_view Array2D<T>::operator[](int r) {
    if (r < 0 || r >= rows)
        throw ::CLRangeError("Array2D::operator[]", "Row index out of range.");
    return row_view(*this, r);
}

template <class T> typename Array2D<T>::const_row_view Array2D<T>::operator[](int r) const {
    if (r < 0 || r >= rows)
        throw ::CLRangeError("Array2D::operator[]", "Row index out of range.");
    return const_row_view(*this, r);
}

template <class T> typename Array2D<T>::col_iterator Array2D<T>::col_begin(int c) {
    if (c < 0 || c >= cols) {
        throw ::CLRangeError("Array2D::col_begin()", "Column index out of range.");
    }
    return col_iterator(rows ? data + c : data, cols);
}

template <class T> typename Array2D<T>::col_iterator Array2D<T>::col_end(int c) {
    if (c < 0 || c >= cols) {
        throw ::CLRangeError("Array2D::col_end()", "Column index out of range.");
    }
    return col_iterator(rows ? data + c : data, cols, rows);
}

template <class T> typename Array2D<T>::const_col_iterator Array2D<T>::col_begin(int c) const {
    if (c < 0 || c >= cols) {
        throw ::CLRangeError("Array2D::col_begin()", "Column index out of range.");
    }
    return const_col_iterator(rows ? data + c : data, cols);
}

template <class T> typename Array2D<T>::const_col_iterator Array2D<T>::col_end(int c) const {
    if (c < 0 || c >= cols) {
        throw ::CLRangeError("Array2D::col_end()", "Column index out of range.");
    }
    return const_col_iterator(rows ? data + c : data, cols, rows);
}

template <class T> void Array2D<T>::getColumn(Array1D<T>& toArray, int c) const {
    if (c < 0 || c >= cols)
        throw ::CLRangeError("Array2D::getColumn()", "Column index out of range.");
    if (toArray.size() != rows)
        toArray = Array1D<T>(rows);
    std::copy(col_begin(c), col_end(c), toArray.begin());
}

template <class T> void Array2D<T>::getRow(Array1D<T>& toArray, int r) const {
    if (r < 0 || r >= rows)
        throw ::CLRangeError("Array2D::getRow()", "Row index out of range.");
    if (toArray.size() != cols)
        toArray = Array1D<T>(cols);
    std::copy(row_begin(r), row_end(r), toArray.begin());
}

template <class T> void Array2D<T>::putColumn(const Array1D<T>& fromArray, int c) {
    if (fromArray.size() != rows) {
        throw ::LogicalError("Array2D::putColumn()", "Dimensions inconsistent.");
    }
    std::copy(fromArray.begin(), fromArray.end(), col_begin(c));
}

template <class T> void Array2D<T>::putRow(const Array1D<T>& fromArray, int r) {
    if (fromArray.size() != cols) {
        throw ::LogicalError("Array2D::putRow()", "Dimensions inconsistent.");
    }
    std::copy(fromArray.begin(), fromArray.end(), row_begin(r));
}

template <class T> void Array2D<T>::swapColumns(int col1, int col2) {
    std::swap_ranges(col_begin(col1), col_end(col1), col_begin(col2));
}

template <class T> void Array2D<T>::swapRows(int row1, int row2) {
    std::swap_ranges(row_begin(row1), row_end(row1), row_begin(row2));
}

template <class T> std::ostream& operator<<(std::ostream& os, const Array2D<T>& v) {
    for (int i = 0; i < v.nrows(); ++i) {
        for (int j = 0; j < v.ncols(); ++j) {
            os << v(i, j) << " ";
        }
        os << std::endl;
    }
    return os;
}

} // namespace tps

#endif // TPS_Array2D_HH
