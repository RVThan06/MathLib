#ifndef ROOTS_HPP
#define ROOTS_HPP

#include <cmath>
#include <concepts>
#include <iostream>
#include <mathlib/Utilities.hpp>
#include <optional>

namespace mathlib::numerical {

/**
 * @brief Newton Rhapson algorithm that checks for roots
 * requires the fucntion and its first derivative.
 * @tparam T is the floating point type for the function value
 * @tparam func a concept restricted type for function to call
 * @tparam derv a concept restricted type for the derivate of the function
 * @param first_guess is the initial guess x coordinate
 * @param function is the function to solve
 * @param derivative is the derivative of the function
 * @param iterations the max iterations to perform, default value = 100
 * @param tolerance the smallest margin needed to accept a root value
 * @return std::optional<T> with the root value if found, else std::nullopt
 */
template <std::floating_point T, mathlib::ActualFucntion<T> func, mathlib::ActualFucntion<T> derv>
std::optional<T> Newton_Rhapson(func function, derv derivative, T first_guess, int iterations = 100,
                                T tolerance = static_cast<T>(1e-6)) {
    T current_guess = first_guess;

    for (int i = 0; i < iterations; ++i) {
        T solution_y = function(current_guess);

        // check if guess is right
        if (std::abs(solution_y) < tolerance) {
            return current_guess;
        }

        T slope_dy = derivative(current_guess);

        // check if the gradient is flat to prevent zero division
        if (std::abs(slope_dy) < EPSILON<T>) {
            return std::nullopt;
        }

        // Newton Rhapson to update current guess
        T step_size = (solution_y / slope_dy);
        current_guess -= step_size;

        // check if convergence is happening
        if (std::abs(step_size) < tolerance) {
            return current_guess;
        }
    }

    // max iterations and failed to converge
    return std::nullopt;
}

/**
 * @brief Bisection method that uses bracketing.
 * Slower than Newton Rhapson but guarantees root if one exist.
 * @tparam T is the floating point type for the function value
 * @tparam func a concept restricted type for function to call
 * @param guess_A is the first x coordinate guess, smaller than guess B
 * @param guess_B is the second x coordinate guess
 * @param function the function to be solved
 * @param iterations the max iterations to perform, default value = 100
 * @param tolerance the smallest margin needed to accept a root value
 * @return std::optional<T> with the root value if found, else std::nullopt
 */

template <std::floating_point T, mathlib::ActualFucntion<T> func>
std::optional<T> Bisection(T guess_A, T guess_B, func function, int iterations = 100,
                           T tolerance = static_cast<T>(1e-6)) {
    T fa = function(guess_A);
    T fb = function(guess_B);

    // both guesses are positive / negative
    if (fa * fb > 0) {
        return std::nullopt;
    }

    for (int i = 0; i < iterations; ++i) {
        T mid_point = guess_A + (guess_B - guess_A) * 0.5; // safer prevents overflow by not adding two numbers
        T fmid = function(mid_point);

        // if we find the root or the interval is very narrow
        if (std::abs(fmid) <= tolerance || (guess_B - guess_A) * 0.5 < tolerance) {
            return mid_point;
        }

        // decide which half to keep
        if (fmid * fa > 0) {
            guess_A = mid_point;
            fa = fmid;
        } else {
            guess_B = mid_point;
        }
    }

    // sufficient iterations, return best guess
    return guess_A + (guess_B - guess_A) * 0.5;
}

/**
 * @brief Secant method, does not require first derivative.
 * @tparam T is the floating point type for the function value
 * @tparam func a concept restricted type for function to call
 * @param x0 is the first x for secant line
 * @param x1 is the second x for secant line
 * @param function the function to be solved
 * @param iterations the max iterations to perform, default value = 100
 * @param tolerance the smallest margin needed to accept a root value
 * @return std::optional<T> with the root value if found, else std::nullopt
 */

template <std::floating_point T, mathlib::ActualFucntion<T> Func>
std::optional<T> Secant(Func function, T x0, T x1, int iterations = 100, T tolerance = static_cast<T>(1e-6)) {
    for (int i = 0; i < iterations; ++i) {
        T f_x0 = function(x0);
        T f_x1 = function(x1);

        // prevent division by zero
        if (std::abs(f_x1 - f_x0) < EPSILON<T>) {
            return std::nullopt;
        }

        // Secant step: x = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        T step = f_x1 * ((x1 - x0) / (f_x1 - f_x0));

        // check for step size convergence
        if (std::abs(step) < tolerance) {
            return x1 - step;
        }

        // update next point
        x0 = x1;
        x1 = x1 - step;
    }

    // no convergence
    return std::nullopt;
}

/**
 * @brief False position algorithm for root finding.
 * @tparam T is the floating point type for the function value
 * @tparam func a concept restricted type for function to call
 * @param guess_A is the lower bound x coordinate (min)
 * @param guess_B is the upper bound x coordinate (max)
 * @param function the function to be solved
 * @param iterations the max iterations to perform, default value = 100
 * @param tolerance the smallest margin needed to accept a root value
 * @return std::optional<T> with the root value if found, else std::nullopt
 */
template <std::floating_point T, mathlib::ActualFucntion<T> func>
std::optional<T> False_Position(T guess_A, T guess_B, func function, int iterations = 100,
                                T tolerance = static_cast<T>(1e-6)) {
    T fa = function(guess_A);
    T fb = function(guess_B);

    // both guesses are positive / negative
    if (fa * fb > 0) {
        return std::nullopt;
    }

    for (int i = 0; i < iterations; ++i) {
        // prevent division by zero
        if (std::abs(fb - fa) < EPSILON<T>) {
            return std::nullopt;
        }

        T mid_point = guess_B - fb * ((guess_B - guess_A) / (fb - fa));
        T fmid = function(mid_point);

        // check if we found the root
        if (std::abs(fmid) < tolerance) {
            return mid_point;
        }

        // decide which half to keep
        if (fmid * fa > 0) {
            guess_A = mid_point;
            fa = fmid;
        } else {
            guess_B = mid_point;
            fb = fmid;
        }
    }

    // sufficient iterations, return best guess
    return (guess_B - fb * ((guess_B - guess_A) / (fb - fa)));
}

} // namespace mathlib::numerical

#endif
