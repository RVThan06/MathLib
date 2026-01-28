#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <algorithm>
#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <initializer_list>
#include <iostream>
#include <mathlib\core\Vector.hpp>
#include <stdexcept>

namespace mathlib::core {

/**
 * @brief A 2D matrix class (R x C)
 * This class template implements 2D matrix where the matrix elements
 * are restricted to floating point numbers.
 * @tparam T is the flaoting point type for matrix elements
 * @tparam Row is the number of rows in matrix
 * @tparam Col is the number of columns in matrix
 */

template <std::floating_point T, std::size_t Row, std::size_t Col>
struct Matrix {
    static_assert(Row > 0 && Col > 0, "Matrix dimension cannot be 0.");

    /// @brief data member to store matrix as 1D array (Row-Major) - zero initialised by default
    std::array<T, Row * Col> m_array{};

    /// @brief data members to store matrix dimension
    static constexpr std::size_t m_rows = Row;
    static constexpr std::size_t m_cols = Col;
    static constexpr std::size_t m_size = Row * Col;

    // --- Accessor for static data member ---

    static constexpr std::size_t get_rows() { return m_rows; }
    static constexpr std::size_t get_cols() { return m_cols; }
    static constexpr std::size_t get_size() { return m_size; }

    // --- Constructors ---

    /// @brief default constructor --> zero intialised matrix
    constexpr Matrix() {};

    /**
     * @brief constructor to intiliased all elements with a given value
     * @param value the value to be used to initialise entire matrix.
     */
    constexpr explicit Matrix(T value) { std::fill(m_array.begin(), m_array.end(), value); }

    /**
     * @brief List constructor to initialise 2D matrix
     * @param list a 2D list with matrix elements
     * @note if 2D list provided is jagged, vector construction will fail, list should match
     * matrix dimension.
     */
    constexpr Matrix(std::initializer_list<std::initializer_list<T>> list) {
        // debug build check
        assert(list.size() == m_rows && "Intilializer list size mismatch with matrix rows.");

        // release build
        if (list.size() != m_rows) {
            throw std::length_error("Intilializer list size mismatch with matrix rows.");
        }

        // COpy values from intializer list after safety check is done
        for (std::size_t i = 0; i < m_rows; i++) {
            // debug check
            assert((list.begin() + i)->size() == m_cols && "Intilializer list size mismatch with matrix columns.");

            // release check
            if ((list.begin() + i)->size() != m_cols) {
                throw std::length_error("Intilializer list size mismatch with matrix columns.");
            }

            for (std::size_t j = 0; j < m_cols; j++) {
                (*this)(i, j) = *((list.begin() + i).begin() + j);
            }
        }
    }

    // --- Paranthesis operator overload ---

    /**
     * @brief Access elements at (row, col)
     * @warning there is no bound checking done
     * @return reference to underlying m_array element of type T
     */
    constexpr T& operator()(std::size_t row, std::size_t col) { return m_array[(row * col) + col]; }
    constexpr const T& operator()(std::size_t row, std::size_t col) const { return m_array[(row * col) + col]; }

    // --- Matrix elements mutator functions ---

    /**
     * @brief allows to set new value at (row, col).
     * @param row: row index
     * @param col: column index
     * @param value: value to be stored at given matrix location
     * @throw std::out_of_range if row/col index is beyond matrix bound
     */
    constexpr void set_value(std::size_t row, std::size_t col, T value) {
        // debug build check
        assert((row < Row && col < Col) && "Row or column is out of matrix index bound.");

        // release build check
        if (row >= Row && col >= Col) {
            throw std::out_of_range("Row or column is out of matrix index bound.");
        }

        // safely set element value
        m_array[(row * col) + col] = value;
    }

    /**
     * @brief Allows to set value for an entire row
     * @param row: the row index
     * @param row_vector: a vector same size as the num of cols
     */
    constexpr void set_row(std::size_t row, Vector<T, m_cols> row_vector) {
        // debug build check
        assert((row < Row) && "Matrix::set_row row exceeds matrix index.");

        // release build check
        if (row >= Row) {
            throw std::out_of_range("Matrix::set_row row exceeds matrix index.");
        }
        for (std::size_t i = 0; i < m_cols; i++) {
            m_array[(row * i) + i] = row_vector[i];
        }
    }

    /**
     * @brief Allows to set value for an entire column
     * @param col: column index
     * @param column_vector: a vector same size as the num of rows
     */
    constexpr void set_column(std::size_t col, Vector<T, m_rows> column_vector) {
        // debug build check
        assert((col < m_cols) && "Matrix::set_column column exceeds matrix index.");

        // release build check
        if (col >= m_cols) {
            throw std::out_of_range("Matrix::set_column column exceeds matrix index.");
        }
        for (std::size_t i = 0; i < m_rows; i++) {
            m_array[(i * col) + col] = column_vector[i];
        }
    }

    /**
     * @brief To extract an entire row of elements for any given row.
     * @param row: the index of row to be extracted
     * @return all elements in a row returned as Mathlib::Vector object
     */
    constexpr Vector<T, Col> get_row_vector(std::size_t row) {
        // debug build check
        assert((row < Row) && "Matrix::get_row_vector row exceeds matrix index.");

        // release build check
        if (row >= Row) {
            throw std::out_of_range("Matrix::get_row_vector row exceeds matrix index.");
        }
        Vector<T, Col> row_vector;
        for (std::size_t i = 0; i < m_cols; i++) {
            row_vector[i] = (*this)(row, i);
        }
        return row_vector;
    }

    /**
     * @brief To extract an entire column of elements for any given column.
     * @param column: the index of column to be extracted
     * @return all elements in a column returned as Mathlib::Vector object
     */
    constexpr Vector<T, Row> get_column_vector(std::size_t column) {
        // debug build check
        assert((column < m_cols) && "Matrix::get_column_vector column exceeds matrix index.");

        // release build check
        if (column >= m_cols) {
            throw std::out_of_range("Matrix::get_column_vector column exceeds matrix index.");
        }
        Vector<T, Row> column_vector;
        for (std::size_t i = 0; i < m_rows; i++) {
            column_vector[i] = (*this)(i, column);
        }
        return column_vector;
    }

    /**
     * @brief TO return a identity matrix for square matrix.
     * @note If called by non-square matrix, compilation will fail due to static_assert
     */
    constexpr Matrix get_identity() const {
        static_assert(m_rows == m_cols, "Matrix is not a square matrix.");
        Matrix identity_matrix;
        for (std::size_t i = 0; i < m_rows; i++) {
            identity_matrix(i, i) = static_cast<T>(1);
        }
        return identity_matrix;
    }

    // --- Compound assignments ---

    /**
     * @brief Compound addition of two matrix, modified LHS matrix.
     * @return reference to LHS matrix.
     */
    constexpr Matrix& operator+=(const Matrix& rhs) {
        for (std::size_t i = 0; i < m_size; i++) {
            m_array[i] += rhs.m_array[i];
        }
        return *this;
    }

    /**
     * @brief Compound subtraction of two matrix, modified LHS matrix.
     * @return reference to LHS matrix.
     */
    constexpr Matrix& operator-=(const Matrix& rhs) {
        for (std::size_t i = 0; i < m_size; i++) {
            m_array[i] -= rhs.m_array[i];
        }
        return *this;
    }

    /**
     * @brief Compound multiplication of matrix with scalar, modifies LHS matrix.
     * @return reference to LHS matrix.
     */
    constexpr Matrix& operator*=(const T& scalar) {
        for (std::size_t i = 0; i < m_size; i++) {
            m_array[i] *= scalar;
        }
        return *this;
    }

    /**
     * @brief Compound division of matrix with scalar, modifies LHS matrix.
     * @throw std::overflow_error if scalar value is smaller than epsilon.
     * @return reference to LHS matrix.
     */
    constexpr Matrix& operator/=(const T& scalar) {
        if (std::abs(scalar) <= EPSILON<T>) {
            throw std::overflow_error("Division by zero.");
        }
        for (std::size_t i = 0; i < m_size; i++) {
            m_array[i] /= scalar;
        }
        return *this;
    }

    // --- Friend arithmetic operators (Canonical implementation) ---

    /**
     * @brief Canonical implementation of arithmetic operators by calling compound assignment.
     * @return a copy of matrix object with the result of the operation.
     */
    friend constexpr Matrix operator+(Matrix lhs, const Matrix& rhs) { return (lhs += rhs); }
    friend constexpr Matrix operator-(Matrix lhs, const Matrix& rhs) { return (lhs -= rhs); }
    friend constexpr Matrix operator*(Matrix lhs, const T& scalar) { return (lhs += scalar); }
    friend constexpr Matrix operator*(const T& scalar, Matrix rhs) { return (rhs += scalar); }
    friend constexpr Matrix operator/(Matrix lhs, const T& scalar) { return (lhs /= scalar); }

    /**
     * @brief Matrix multiplication (R x C) x (C x H) = (R x H)
     * Uses member templates - fried function that is also a template
     * allows for another template arguement besides class template arguements.
     * @return A copy of matrix object of dimension (R x H)
     */
    template <std::size_t SecondCols>
    friend constexpr Matrix<T, Row, SecondCols> operator*(const Matrix& lhs, const Matrix<T, Col, SecondCols>& rhs) {
        Matrix<T, Row, SecondCols> result;

        for (std::size_t x = 0; x < m_rows; x++) {
            for (std::size_t y = 0; y < SecondCols; y++) {
                T value{};
                for (std::size_t z = 0; z < m_cols; z++) {
                    value = value + (lhs(x, z) * rhs(z, y));
                }
                result(x, y) = value;
            }
        }
        return result;
    }

    /**
     * @brief Matrix multiplication with column vector (R x C) x (C x 1) = (R x 1)
     * @return A copy of Vector object of size R.
     */
    constexpr Vector<T, Row> multiply_column(Vector<T, Col> column_vector) const {
        Vector<T, Row> result;

        for (std::size_t x = 0; x < m_rows; x++) {
            T value{};
            for (std::size_t z = 0; z < m_cols; z++) {
                value = value + ((*this)(x, z) * column_vector[z]);
            }
            result(x, y) = value;
        }
        return result;
    }

    // --- Ostream overload ---

    friend std::ostream& operator<<(std::ostream& out, const Matrix& matrix) {
        for (std::size_t i = 0; i < Row; i++) {
            out << "[ ";
            for (std::size_t j = 0; j < Col; j++) {
                out << matrix(i, j) << (j < Col - 1 ? ", " : " ");
            }
            out << "]\n";
        }

        return out;
    }

    // --- Matrix operation ---

    /**
     * @brief Hadamart product which performs element wise multiplication.
     * @note Matrix dimension should match.
     * @return A copy of matrix object with result of multiplication.
     */
    constexpr Matrix Hadamard_product(const Matrix& rhs) const {
        Matrix result;
        for (std::size_t i = 0; i < Row; i++) {
            for (std::size_t j = 0; j < Col; j++) {
                result.m_array[i] = m_array[i] * rhs.m_array[i];
            }
        }
        return result;
    }

    /**
     * @brief transpose to swap row and column of current matrix
     * @return a copy of matrix object with rows and columns swapped.
     */
    constexpr Matrix<T, Col, Row> transpose() const {
        Matrix<T, Col, Row> result;

        for (std::size_t i = 0; i < Row; i++) {
            for (std::size_t j = 0; j < Col; j++) {
                result(i, j) = (*this)(j, i);
            }
        }
        return result;
    }

    /**
     * @brief Reshapes matrix into a new valid shape.
     * @note reshape dimension has to preserve total matrix elements count.
     * @return a reshaped matrix with new dimension but same element count.
     */
    template <std::size_t New_Row, std::size_t New_Col>
    constexpr Matrix<T, New_Row, New_Col> reshape() const {
        // compile time check
        static_assert(New_Row * New_Col == m_size, "Matrix::reshape must preserve total elements count.");
        Matrix<T, New_Row, New_Col> result;
        result.m_array = m_array;
        return result;
    }

    /**
     * @brief Extracts submatrix from current matrix.
     * @return a submatrix of the same dimension or smaller than current matrix.
     * @throw std::out_of_range
     */
    template <std::size_t sub_row, std::size_t sub_column>
    constexpr Matrix<T, sub_row, sub_column> sub_matrix(std::size_t row_start, std::size_t column_start) const {
        // compile time check of submatrix size
        static_assert(sub_row <= Row && sub_column <= Col,
                      "Submatrix dimension must be within original matrix dimension.");
        // debug build check
        assert(sub_row + row_start <= Row && sub_column + column_start <= Col);

        // release build check
        if (sub_row + row_start > Row || sub_column + column_start > Col) {
            throw std::out_of_range("Matrix::sub_matrix index exceeds matrix dimension.");
        }

        Matrix<T, sub_row, sub_column> sub_matrix;
        for (std::size_t i = 0; i < sub_row; i++) {
            for (std::size_t j = 0; j < sub_column; j++) {
                sub_matrix(i, j) = (*this)(i + row_start, j + column_start);
            }
        }
        return sub_matrix;
    }
};
} // namespace mathlib::core

#endif
