#include <cstddef>
#include <iomanip>
#include <iostream>
#include <mathlib\Utilities.hpp>
#include <mathlib\numerical\Roots.hpp>
#include <optional>

int main() {
    // /*******************************************************************************************
    //  * ROOT FINDING
    //  * *****************************************************************************************
    //  */

    // auto poly_func = [](double x) {
    //     return std::pow(x, 4) + (static_cast<double>(3.0) * std::pow(x, 3)) + static_cast<double>(5.0) * x;
    // };

    // auto first_deriv = [](double x) {
    //     return static_cast<double>(4.0) * std::pow(x, 3) + (static_cast<double>(9.0) * std::pow(x, 2)) +
    //            static_cast<double>(5.0);
    // };

    auto poly_func_2 = [](double x) { return x * x + (static_cast<double>(2.0) * x) - static_cast<double>(15.0); };

    auto first_deriv_2 = [](double x) { return static_cast<double>(2.0) * x + 2; };

    std::cout << "Regular polynomial integration\n\n";
    std::optional<double> solution = mathlib::numerical::Newton_Rhapson(poly_func_2, first_deriv_2, 6.0, 200);
    if (solution) {
        std::cout << "Newton_Rhapson solution: " << std::setprecision(10) << *solution << "\n";
    } else {
        std::cout << "Newton_Rhapson solution:: No roots found\n";
    }

    solution = mathlib::numerical::Bisection(-1.0, 6.0, poly_func_2);
    if (solution) {
        std::cout << "Bisection solution: " << std::setprecision(10) << *solution << "\n";
    } else {
        std::cout << "Bisection solution:: No roots found\n";
    }

    solution = mathlib::numerical::Secant(poly_func_2, 6.0, -1.0, 200);
    if (solution) {
        std::cout << "Secant solution: " << std::setprecision(10) << *solution << "\n";
    } else {
        std::cout << "Secant solution:: No roots found\n";
    }

    solution = mathlib::numerical::False_Position(-1.0, 6.0, poly_func_2);
    if (solution) {
        std::cout << "Secant solution: " << std::setprecision(10) << *solution << "\n";
    } else {
        std::cout << "Secant solution:: No roots found\n";
    }
    return 0;
}
