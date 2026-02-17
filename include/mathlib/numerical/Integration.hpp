#ifndef INTEGRATION_HPP
#define INTEGRATION_HPP

#include <cmath>
#include <concepts>
#include <iostream>
#include <mathlib/Utilities.hpp>

namespace mathlib::numerical {

/**
 * @brief Trapezoidal integration method.
 * @tparam  T -> the floating point type for function parameter.
 * @tparam func -> the function to be differentiated.
 * @param function -> the function f(x).
 * @param start -> the lower limit of integration
 * @param end -> the upper limit of integration
 * @param intervals -> number of sub-intervals for integration, higher value for better resolution.
 * @return -> the evaluated value of integration of the function at point x.
 */
template <std::floating_point T, mathlib::ActualFucntion<T> func>
T trapezoidal_integral(func function, T start, T end, int intervals) {
    T height = (end - start) / static_cast<T>(intervals);

    T sum = (function(start) + function(end)) * static_cast<T>(0.5);
    for (int i = 1; i < intervals; ++i) {
        T x = start + (i * height);
        sum += function(x);
    }

    return sum * height;
}

/**
 * @brief Simpson 1/3 integration method.
 * @tparam  T -> the floating point type for function parameter.
 * @tparam func -> the function to be differentiated.
 * @param function -> the function f(x).
 * @param start -> the lower limit of integration
 * @param end -> the upper limit of integration
 * @param intervals -> number of sub-intervals for integration, higher value for better resolution.
 * @note The number of intervals must be even or else will be incremented by one to make it even.
 * @return -> the evaluated value of integration of the function at point x.
 */
template <std::floating_point T, mathlib::ActualFucntion<T> func>
T simpson_one_three(func function, T start, T end, int intervals) {
    // check if interval is even
    if (intervals % 2 != 0) {
        intervals += 1;
    }

    T height = (end - start) / static_cast<T>(intervals);
    T sum = (function(start) + function(end));

    for (int i = 1; i < intervals; ++i) {
        T x = start + (i * height);

        if (i % 2 == 0) {
            sum += (static_cast<T>(4) * function(x));
            continue;
        }
        sum += (static_cast<T>(2) * function(x));
    }

    return (sum * (height / static_cast<T>(3)));
}

/**
 * @brief Simpson 3/8 integration method.
 * @tparam  T -> the floating point type for function parameter.
 * @tparam func -> the function to be differentiated.
 * @param function -> the function f(x).
 * @param start -> the lower limit of integration
 * @param end -> the upper limit of integration
 * @param intervals -> number of sub-intervals for integration, higher value for better resolution.
 * @note The number of intervals must be multiple of 3or else will be incremented to make it a multiple of 3.
 * @return -> the evaluated value of integration of the function at point x.
 */
template <std::floating_point T, mathlib::ActualFucntion<T> func>
T simpson_three_eight(func function, T start, T end, int intervals) {
    // check if interval is multiple of 3
    if (intervals % 3 != 0) {
        intervals += (3 - (intervals % 3));
    }

    // get the height
    T height = (end - start) / static_cast<T>(intervals);

    // sum the endpoints
    T sum = (function(start) + function(end));

    for (int i = 1; i < intervals; ++i) {
        T x = start + (i * height);

        if (i % 3 == 0) {
            sum += (static_cast<T>(2) * function(x));
            continue;
        }
        sum += (static_cast<T>(3) * function(x));
    }
    return (sum * (static_cast<T>(3) * height / static_cast<T>(8)));
}

} // namespace mathlib::numerical

#endif
