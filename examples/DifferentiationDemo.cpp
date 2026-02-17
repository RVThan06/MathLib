#include <cstddef>
#include <iomanip>
#include <iostream>
#include <mathlib/Utilities.hpp>
#include <mathlib/numerical/Differentiation.hpp>

int main() {
    /*******************************************************************************************
     * DIFFERENTIATION
     * *****************************************************************************************
     */
    auto poly_func = [](double x) {
        return std::pow(x, 4) + (static_cast<double>(3.0) * std::pow(x, 3)) + static_cast<double>(5.0) * x;
    };

    std::cout << "Regular polynomial differentiation\n\n";
    double solution = mathlib::numerical::central_difference(poly_func, 1.5);
    std::cout << "Central difference solution: " << std::setprecision(10) << solution << "\n";

    solution = mathlib::numerical::richardson_extrapolation(poly_func, 1.5);
    std::cout << "Richardson extrapolation solution: " << std::setprecision(10) << solution << "\n";

    solution = mathlib::numerical::second_derivative(poly_func, 1.5);
    std::cout << "Second derivative solution: " << std::setprecision(10) << solution << "\n\n";

    auto trigo_func = [](double x) { return std::sin(x * x) + std::cos(std::pow(x, 3)); };

    std::cout << "Trogonometric function differentiation\n\n";
    solution = mathlib::numerical::central_difference(trigo_func, 1.5);
    std::cout << "Central difference solution: " << std::setprecision(10) << solution << "\n";

    solution = mathlib::numerical::richardson_extrapolation(trigo_func, 1.5);
    std::cout << "Richardson extrapolation solution: " << std::setprecision(10) << solution << "\n";

    solution = mathlib::numerical::second_derivative(trigo_func, 1.5);
    std::cout << "Second derivative solution: " << std::setprecision(10) << solution << "\n";

    return 0;
}
