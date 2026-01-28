#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <algorithm>
#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <initializer_list>
#include <iostream>
#include <mathlib\Utilities.hpp>
#include <stdexcept>

namespace mathlib::core {

/**
 * @brief Multidimensional Vector class.
 * Rhis class template implements multidimensional vector where the vector elements
 * are restricted to floating point numbers.
 * @tparam T is the type for vector elements, restricted to floating point.
 * @tparam N the dimension of the vector, at least 1
 */

template <std::floating_point T, std::size_t N>
struct Vector {
    static_assert(N > 0, "Vector dimension cannot be zero.");

    /// @brief data members to store vector and its dimension
    std::array<T, N> m_array{};
    static constexpr std::size_t m_dimensions = N;

    /// @brief access for static data memeber
    static constexpr std::size_t get_dimension() { return m_dimensions; }

    // --- Constructors ---

    /// @brief Default constructor - all magnitude intialised to zero
    constexpr Vector() {};

    /**
     * @brief Constructor with arguement
     * This constructor will initialise all vector dimension with the same value.
     * @param magnitude scalar value for vector dimensions.
     */
    constexpr explicit Vector(T magnitude) { std::fill(m_array.begin(), m_array.end(), magnitude); }

    /**
     * @brief List constructor
     * @note If a list shorter than vector dimension is provided, the higher dimensions will
     * have zero magnitude since the underlyig array is zero intialised before copying.
     * @param list the intialiser list with vector magnitude for each dimension.
     */
    constexpr Vector(std::initializer_list<T> list) {
        // debug build check
        assert(list.size() <= m_dimensions && "Initializer list is larger than vector dimension.");

        // release build check
        if (list.size() <= m_dimensions) {
            throw std::length_error("Initializer list is larger than vector dimension.");
        }
        std::copy(list.begin(), list.end(), m_array.begin());
    }

    // --- Subscript overloads ---

    /**
     * @brief subscript overload for const and non-const vector objects
     * @warning there is no bound checking done
     * @return reference to underlying m_array element
     */
    constexpr T& operator[](std::size_t index) { return m_array[index]; }
    constexpr const T& operator[](std::size_t index) const { return m_array[index]; }

    // --- Unary operators ---

    /**
     * @brief Negation operator to invert sign of all vector dimensions.
     * @return a new vector object with sign of each dimension inverted.
     */
    constexpr Vector operator-() const {
        Vector<T, N> temp;
        for (std::size_t i = 0; i < N; i++) {
            temp[i] = -m_array[i];
        }
        return temp;
    }

    constexpr Vector operator+() const { return *this; }

    // --- Compound Assignments ---

    /**
     * @brief compound addition of two vector objects.
     * Element wise addition of each dimension to modify LHS operand.
     * @return a reference to the LHS operand.
     */
    constexpr Vector& operator+=(const Vector& rhs) {
        for (std::size_t i = 0; i < N; i++) {
            m_array[i] += rhs[i];
        }
        return *this;
    }

    /**
     * @brief compound subtraction of two vector objects.
     * Element wise subtraction of each dimension to modify LHS operand.
     * @return a reference to the LHS operand.
     */
    constexpr Vector& operator-=(const Vector& rhs) {
        for (std::size_t i = 0; i < N; i++) {
            m_array[i] -= rhs[i];
        }
        return *this;
    }

    /**
     * @brief compound multiplication of two vector objects.
     * Element wise multiplication of each dimension to modify LHS operand.
     * @return a reference to the LHS operand.
     */
    constexpr Vector& operator*=(const Vector& rhs) {
        for (std::size_t i = 0; i < N; i++) {
            m_array[i] *= rhs[i];
        }
        return *this;
    }

    /**
     * @brief compound multiplication of vector and scalar
     * Element wise multiplication of each dimension with scalar to modify LHS operand.
     * @return a reference to the LHS operand.
     */
    constexpr Vector& operator*=(const T& scalar) {
        for (std::size_t i = 0; i < N; i++) {
            m_array[i] *= scalar;
        }
        return *this;
    }

    /**
     * @brief compound division of two vectors
     * Element wise division of each dimension to modify LHS operand.
     * @return a reference to the LHS operand.
     */
    constexpr Vector& operator/=(const Vector& rhs) {
        for (std::size_t i = 0; i < N; i++) {
            if (std::abs(rhs[i]) <= EPSILON<T>) {
                throw std::overflow_error("Division by zero.");
            }
            m_array[i] /= rhs[i];
        }
        return *this;
    }

    /**
     * @brief compound division of vector with scalar value.
     * Element wise division of each dimension to modify LHS operand.
     * @return a reference to the LHS operand.
     */
    constexpr Vector& operator/=(const T& scalar) {
        if (std::abs(scalar) <= EPSILON<T>) {
            throw std::overflow_error("Division by zero.");
        }
        for (std::size_t i = 0; i < N; i++) {
            m_array[i] /= scalar;
        }
        return *this;
    }

    // --- Arithmetic operators (canonical implementation) ---

    /**
     * @brief Canonical implementation of arithmetic operators by calling compound assignment.
     * @return a copy of vector object with the result of the operation.
     */
    friend constexpr Vector operator+(Vector lhs, const Vector& rhs) { return (lhs += rhs); }
    friend constexpr Vector operator-(Vector lhs, const Vector& rhs) { return (lhs -= rhs); }
    friend constexpr Vector operator*(Vector lhs, const Vector& rhs) { return (lhs *= rhs); }
    friend constexpr Vector operator/(Vector lhs, const Vector& rhs) { return (lhs /= rhs); }
    friend constexpr Vector operator*(Vector lhs, const T& scalar) { return (lhs *= scalar); }
    friend constexpr Vector operator*(const T& scalar, Vector rhs) { return (rhs *= scalar); }
    friend constexpr Vector operator/(Vector lhs, const T& scalar) { return (lhs /= scalar); }

    // --- Equality check ---

    /**
     * @brief equality check of two vector objects.
     * Checks if each dimension of vector is the same.
     * @return true - vector values match/ false = otherwise
     */
    constexpr bool operator==(const Vector& rhs) const {
        // debug build check
        assert(m_array.size() == rhs.m_array.size() && "vector dimensions should match for comparison.");

        // release build check
        if (m_array.size() != rhs.m_array.size()) {
            throw std::out_of_range("vector dimensions should match for comparison.");
        }

        for (std::size_t i = 0; i < N; i++) {
            if (!mathlib::equals(m_array[i], rhs.m_array[i])) {
                return false;
            }
        }
        return true;
    }

    constexpr bool operator!=(const Vector& rhs) const { return !(*this == rhs); }

    // --- Ostream overload for printing ---

    friend std::ostream& operator<<(std::ostream& out, const Vector& my_vector) {
        out << "[ ";
        for (std::size_t i = 0; i < N; i++) {
            out << my_vector[i] << (i == N - 1 ? "" : ", ");
        }
        out << " ]";
        return out;
    }

    // --- Vector mathematical operations ---

    /**
     * @brief Computes dot product of two vectors
     * @return a floating point result of dot product.
     */
    constexpr T dot(const Vector& rhs) const {
        if (m_array.size() != rhs.m_array.size()) {
            throw std::out_of_range("vector dimensions should match for dot product.");
        }

        T total = 0;
        for (std::size_t i = 0; i < N; i++) {
            total += (m_array[i] * rhs[i]);
        }
        return total;
    }

    /**
     * @brief computes the unit vector
     * @returns the unit vector as a Vector object by copy.
     */
    Vector get_unit_vector() const {
        if (magnitude() <= EPSILON<T>) {
            return 0;
        }
        return (*this / magnitude());
    }

    /**
     * @brief Computes the angle between two vectors.
     * @return the angle computed in radians.
     */
    T angle_between(const Vector& rhs) const {
        T denominator = magnitude() * rhs.magnitude();

        if (denominator <= EPSILON<T>) {
            return T{0};
        }

        T cos_theta = dot(rhs) / denominator;
        // Clamp value to [-1, 1] to avoid NaN from acos due to precision errors
        // e.g. dot product might result in 1.0000001 which breaks acos
        cos_theta = std::clamp(cos_theta, static_cast<T>(-1.0), static_cast<T>(1.0));
        return std::acos(cos_theta);
    }

    /**
     * @brief Computes the projection of LHS vector onto RHS vector
     * @returns a projection vector by copy.
     */
    Vector projection(const Vector& rhs) const {
        if (rhs.magnitude() <= EPSILON<T>) {
            return Vector(0);
        }
        T mag = magnitude();
        return ((dot(rhs) / (mag * mag)) * rhs);
    }

    /**
     * @brief Computes the distance between two vectors.
     * @return the magitude of the distance as floating point value.
     */
    T distance_between(const Vector& rhs) const { return (*this - rhs).magnitude(); }

    /// @brief Computes the magnitude of a vector and returns it.
    T magnitude() const { return std::sqrt((*this)->dot(*this)); }

    /// @brief To perform inplace normalisation of vector
    void normalise() { *this = get_unit_vector(); }

    /// @brief Computes the rejection of a vector --> opposite of projection
    Vector rejection(const Vector& rhs) const { return *this - projection(rhs); }
};

} // namespace mathlib::core

#endif
