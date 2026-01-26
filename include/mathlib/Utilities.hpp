#ifndef COMMON_HPP
#define COMMON_HPP

#include <cmath>
#include <concepts>
#include <limits>
#include <numbers>

namespace mathlib {

/**
 * @brief A math constant and utilities function module to assist
 * with mathematical calculation.
 */

/**
 * @brief Concept to restrict parameters to only be of arithmetic type.
 * Only integral and floating point parameters can be used.
 */
template <typename T>
concept arithmetic = std::is_arithmetic_v<T>;

// --- Constants ---

/**
 * @brief Mathematical constant Pi.
 * It is a variable template based on precision of the parameter type.
 * @tparam Floating point type, default = double
 */
template <std::floating_point T = double>
constexpr T PI = std::numbers::pi_v<T>;

/**
 * @brief Mathematical constant Two Pi.
 * Used in diameter calculation where 2 * pi is needed.
 * @tparam Floating point type, default = double
 */
template <std::floating_point T = double>
constexpr T PI_TWO = static_cast<T>(2.0) * PI<T>;

/**
 * @brief Mathematical constant Half of Pi.
 * Used in diameter calculation where pi / 2 is needed.
 * @tparam Floating point type, default = double
 */
template <std::floating_point T = double>
constexpr T PI_HALF = static_cast<T>(0.5) * PI<T>;

/**
 * @brief Mathematical constant Two Pi.
 * Used in diameter calculation where 2 * pi is needed.
 * @tparam Floating point type, default = double
 */
template <std::floating_point T = double>
constexpr T EPSILON = std::numeric_limits<T>::epsilon();

// --- Utilities ---

/**
 * @brief Convertor of angle from degrees to radians.
 * @tparam Floating point type, default = double
 */
template <std::floating_point T = double>
constexpr T to_radians(T degrees) {
    return (degrees * (PI<T> / static_cast<T>(180.0)));
}

/**
 * @brief Convertor of angle from radians to degree.
 * @tparam Floating point type, default = double
 */
template <std::floating_point T = double>
constexpr T to_degree(T radians) {
    return (radians * (static_cast<T>(180.0) / PI<T>));
}

/**
 * @brief Flaoting point numbers comparison.
 * Compares two floating point number within the error margin (epsilon).
 * @tparam Floating point type, default = double
 * @param obj1 first number/obj
 * @param obj2 second number/obj
 * @return True if (obj1 - obj2) less than margin of error (epsilon).
 */
template <std::floating_point T = double>
constexpr bool equals(T obj1, T obj2, T epsilon = EPSILON<T>) {
    return std::abs(obj1 - obj2) <= (epsilon * static_cast<T>(100.0));
}
} // namespace mathlib

#endif
