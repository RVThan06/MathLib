#include <iostream>

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
    std::cout << "Conjugate of Number 1: " << cplx_num1.magnitude() << "\n";
    std::cout << "Conjugate of Number 1: " << cplx_num1.angle() << "\n"; // in radians

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
        std::cerr << "Error: " << e.what() << '\n';
    }

    return 0;
}
