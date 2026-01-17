#include <iomanip>
#include <iostream>

#include "mathlib\Common.hpp"
#include "mathlib\core\Complex.hpp"

int main() {
    // Create complex
    constexpr mathlib::core::Complexd cplx_num1(1.5, 2.5);
    constexpr mathlib::core::Complexd cplx_num2(2.0, -0.8);

    std::cout << "Number 1: " << cplx_num1 << "\n";
    std::cout << "Number 2: " << cplx_num2 << "\n";

    // Test out mathematical operators
    std::cout << "Number 1 + Number 2: " << cplx_num1 + cplx_num2 << "\n";
    std::cout << "Number 1 - Number 2: " << cplx_num1 - cplx_num2 << "\n";
    std::cout << "Number 1 * Number 2: " << cplx_num1 * cplx_num2 << "\n";
    std::cout << "Number 1 / Number 2: " << cplx_num1 / cplx_num2 << "\n";

    std::cout << "Number 1 * 6: " << cplx_num1 * 6 << "\n";
    std::cout << "Number 1 / 6: " << cplx_num1 / 6 << "\n";

    // mathematical operations
    std::cout << "Conjugate of Number 1: " << cplx_num1.conjugate() << "\n";
    std::cout << "Magnitude of Number 1: " << cplx_num1.magnitude() << "\n";
    std::cout << "Angle(degree) of Number 1: " << cplx_num1.angle() << "\n"; // in radians
    std::cout << "Angle(radian) of Number 1: " << mathlib::to_degree<double>(cplx_num1.angle()) << "\n";

    // division by zero attempt
    try {
        mathlib::core::Complexd empty_num;
        std::cout << cplx_num1 / empty_num;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
    }

    try {
        std::cout << cplx_num1 / 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n\n";
    }

    std::cout << "PI (float): " << std::setprecision(20) << mathlib::PI<float> << "\n";
    std::cout << "PI (double): " << std::setprecision(20) << mathlib::PI<double> << "\n";

    return 0;
}
