#ifndef DIFFERENTIATION_H
#define DIFFERENTIATION_H

#include <concepts>
#include <iostream>
#include <mathlib/Utilities.hpp>

namespace mathlib::numerical {
/**
 * @brief Differentiation of function using central difference method.
 * Formula: f(x + h) - f(x - h) / 2h.
 * @tparam T -> the floating point type for function parameter.
 * @tparam func -> the function to be differentiated.
 * @param function -> the function f(x).
 * @param x -> the point on the function where the derivative is computed.
 * @param h -> step size of the function.
 * @return returns then calculated derivative at x.
 */
template <std::floating_point T, mathlib::ActualFucntion<T> func>
T central_difference(func function, T x, T h = static_cast<T>(1e-6)) {
    return (function(x + h) - function(x - h)) / (h * static_cast<T>(2.0));
}

/**
 * @brief Richardson extrapolation for more accurate differentiation.
 * The step size h is much larger since the truncation error from approximating the gradient
 * of tangent is cancelled out, so larger step size is to prevent floating point round off error.
 * Formula: 4D(h/2) - D(h)/3.
 * @tparam T -> the floating point type for function parameter.
 * @tparam func -> the function to be differentiated.
 * @param function -> the function f(x).
 * @param x -> the point on the function where the derivative is computed.
 * @param h -> step size of the function.
 * @return returns then calculated derivative at x
 */
template <std::floating_point T, mathlib::ActualFucntion<T> func>
T richardson_extrapolation(func function, T x, T h = static_cast<T>(0.1)) {
    // approximation with intiial step h
    T deriv = (function(x + h) - function(x - h)) / (h * static_cast<T>(2.0));

    // aprroximation with h/2
    T h_half = h * static_cast<T>(0.5);
    T deriv_2 = (function(x + h_half) - function(x - h_half)) / (h_half * static_cast<T>(2.0));

    return ((static_cast<T>(4.0) * deriv_2) - deriv) / static_cast<T>(3.0);
}

/**
 * @brief Calculates the second derivative of the function.
 * Formula: (f(x + h) - 2f(x) + f(x - h))/ h^2.
 * @return the second derivative at point x
 */
template <std::floating_point T, mathlib::ActualFucntion<T> func>
T second_derivative(func function, T x, T h = static_cast<T>(1e-6)) {
    return (function(x + h) - (static_cast<T>(2.0) * function(x)) + function(x - h)) / (h * h);
}
} // namespace mathlib::numerical

#endif
