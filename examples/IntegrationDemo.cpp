#include <cstddef>
#include <iomanip>
#include <iostream>
#include <mathlib/Utilities.hpp>
#include <mathlib/numerical/Integration.hpp>

int main() {
    /*******************************************************************************************
     * INTEGRATION
     * *****************************************************************************************
     */

    auto poly_func = [](double x) {
        return std::pow(x, 4) + (static_cast<double>(3.0) * std::pow(x, 3)) + static_cast<double>(5.0) * x;
    };

    std::cout << "Regular polynomial integration\n\n";
    double solution = mathlib::numerical::trapezoidal_integral(poly_func, 2.0, 5.0, 1000);
    std::cout << "Trapezoidal integration solution: " << std::setprecision(10) << solution << "\n";

    solution = mathlib::numerical::simpson_one_three(poly_func, 2.0, 5.0, 1000);
    std::cout << "Simpson's 1/3 rule solution: " << std::setprecision(10) << solution << "\n";

    solution = mathlib::numerical::simpson_three_eight(poly_func, 2.0, 5.0, 1000);
    std::cout << "Simpson's 3/8 rule solution: " << std::setprecision(10) << solution << "\n\n";

    auto trigo_func = [](double x) { return std::sin(x * x) + std::cos(std::pow(x, 3)); };

    std::cout << "Trogonometric function Integration\n\n";
    solution = mathlib::numerical::trapezoidal_integral(trigo_func, 0.5, 1.2, 1000);
    std::cout << "Trapezoidal integration solution: " << std::setprecision(10) << solution << "\n";

    solution = mathlib::numerical::simpson_one_three(trigo_func, 0.5, 1.2, 1000);
    std::cout << "Simpson's 1/3 rule solution: " << std::setprecision(10) << solution << "\n";

    solution = mathlib::numerical::simpson_three_eight(trigo_func, 0.5, 1.2, 1000);
    std::cout << "Simpson's 3/8 rule solution: " << std::setprecision(10) << solution << "\n\n";
    return 0;
}
