#ifndef COMPLEX_HPP
#define COMPLEX_HPP

#include <cassert>
#include <cmath>
#include <concepts>
#include <iostream>
#include <mathlib\Common.hpp>
#include <stdexcept>

namespace mathlib::core {

/**
 * @brief A templated complex number class (a + bi)
 * * This class implements a complex number class where
 * it restricts the implementation to floating-point (float, double, long double)
 * @tparam T The underlying floating-point type. Must satisfy `std::floating_point`.
 */

template <std::floating_point T>
struct Complex {
    // --- Constructors ---

    /** @brief Default constructor. Initialises 0 + 0i*/
    constexpr Complex() : m_real{0}, m_imag{0} {}

    /**
     * @brief Constructs a complex number with only real component
     * @param real is the real component and imaginary is set to 0.
     */
    constexpr Complex(T real) : m_real{real}, m_imag{0} {}

    /**
     * @brief Constructs both real and imaginary using user provided value
     * @param real is the real component.
     * @param imag is the imaginary component.
     */
    constexpr Complex(T real, T imag) : m_real{real}, m_imag{imag} {}

    // --- Equality check ---

    /**
     * @brief Equality operators to compare two complex number.
     * @note uses fuzzy comparison to ensure safe comparison of two complex number
     * within an acceptable margin of difference
     * @return true/false
     */
    bool operator==(const Complex& rhs) const {
        if (mathlib::equals(m_real, rhs.m_real) && mathlib::equals(m_imag, rhs.m_imag)) {
            return true;
        }
        return false;
    }

    bool operator!=(const Complex& rhs) const { return !(*this == rhs); }

    // --- Unary Operators ---

    /** @brief Returns the negation of a complex number (-a, -bI) */
    constexpr Complex operator-() const { return Complex(-m_real, -m_imag); }
    constexpr Complex operator+() const { return Complex(m_real, m_imag); }

    // --- Compound Assignment Operators ---

    /**
     * @brief Adds a complex number to this number.
     * This causes the modification of the left hand side (LHS) operand.
     * @return a reference to the LHS object
     */
    constexpr Complex& operator+=(const Complex rhs_complx) {
        m_real += rhs_complx.m_real;
        m_imag += rhs_complx.m_imag;
        return *this;
    }

    /**
     *  @brief Subtracts a complex number from this. Similar implementation to addition.
     */
    constexpr Complex& operator-=(const Complex rhs_complx) {
        m_real -= rhs_complx.m_real;
        m_imag -= rhs_complx.m_imag;
        return *this;
    }

    /**
     * @brief Multiplies this complex number by another.
     * @note Uses the mathematical formula: : (a+bi)(c+di) = (ac-bd) + (ad+bc)i.
     * Modifies the feft hand side (LHS) object to store the result value.
     * @return Aa reference to the LHS object.
     */
    constexpr Complex& operator*=(const Complex rhs_complx) {
        T temp_real = (m_real * rhs_complx.m_real) - (m_imag * rhs_complx.m_imag);
        T temp_img = (m_real * rhs_complx.m_imag) + (m_imag * rhs_complx.m_real);
        m_real = temp_real;
        m_imag = temp_img;
        return *this;
    }

    /**
     * @brief Multiples this complex number with a scalar value.
     */
    constexpr Complex& operator*=(const T& scalar) {
        m_real *= scalar;
        m_imag *= scalar;
        return *this;
    }

    /**
     * @brief Divides this complex number by another.
     * @note Performs division by multiplying by the conjugate denominator. Also
     * checks for zero division error
     * @throws std::runtime_error if tried to divide by 0
     */
    constexpr Complex& operator/=(const Complex& rhs_complex) {
        // debug build check
        assert(rhs_complex.m_real > EPSILON<T> && rhs_complex.m_imag > EPSILON<T>);
        // release build check
        if (rhs_complex.m_real < EPSILON<T> || rhs_complex.m_imag < EPSILON<T>) {
            throw std::overflow_error("Division by zero");
        }

        Complex conj(rhs_complex.m_real, -(rhs_complex.m_imag));
        Complex numerator = Complex(m_real, m_imag) * conj;
        Complex denominator = rhs_complex * conj;

        m_real = numerator.m_real / denominator.m_real;
        m_imag = numerator.m_imag / denominator.m_real;
        return *this;
    }

    /**
     * @brief Divides this complex number with a scalar value.
     * @note performs division by simply diving real and imaginary part by scalar.
     * And performs zero division check.
     * @throws std::runtime_error if tried to divide by 0
     */
    constexpr Complex& operator/=(T scalar) {
        // debug build check
        assert(scalar > EPSILON<T>);
        // release build check
        if (scalar < EPSILON<T>) {
            throw std::overflow_error("Division by zero");
        }

        m_real /= scalar;
        m_imag /= scalar;
        return *this;
    }

    // --- Binary Operators (Hidden Friends) ---

    friend constexpr Complex operator+(Complex lhs, const Complex& rhs) { return (lhs += rhs); }
    friend constexpr Complex operator-(Complex lhs, const Complex& rhs) { return (lhs -= rhs); }
    friend constexpr Complex operator*(Complex lhs, const Complex& rhs) { return (lhs *= rhs); }
    friend constexpr Complex operator*(Complex lhs, T scalar) { return (lhs *= scalar); }
    friend constexpr Complex operator*(T scalar, Complex rhs) { return (rhs *= scalar); }
    friend constexpr Complex operator/(Complex lhs, const Complex& rhs) { return (lhs /= rhs); }
    friend constexpr Complex operator/(Complex lhs, T scalar) { return (lhs /= scalar); }

    // --- Ostream output operator overload ---

    friend std::ostream& operator<<(std::ostream& out, const Complex& number) {
        out << number.m_real << (number.m_imag >= 0 ? " + " : " - ") << std::abs(number.m_imag) << " i";
        return out;
    }

    // --- Mathematical Utilities ---

    /**
     * @brief Returns the complex conjugate.
     * @return A new complex number (real, -imag).
     */
    constexpr Complex conjugate() const { return Complex(m_real, -m_imag); }

    /**
     * @brief Calculates the magnitude (or modulus) of the complex number.
     * @return sqrt{real^2 + imag^2}
     * Use cmath function std::sqrt so it is not constexpr function.
     */
    T magnitude() const { return std::sqrt((m_real * m_real) + (m_imag * m_imag)); }

    /**
     * @brief Calculates the phase angle (argument).
     * @return The angle in **radians** in the range (-\pi, \pi].
     */
    T angle() const { return std::atan2(m_imag, m_real); }

    // --- Data memebers ---

    /// @brief The real part of the complex number.
    T m_real;

    /// @brief The real part of the complex number.
    T m_imag;
};

// --- Type Aliases ---

/// @brief Alias for single precision complex numbers.
using Complexf = Complex<float>;

/// @brief Alias for double precision complex numbers.
using Complexd = Complex<double>;

} // namespace mathlib::core

#endif
