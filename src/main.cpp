#include <cstddef>
#include <iomanip>
#include <iostream>
#include <mathlib\Utilities.hpp>
#include <mathlib\core\Complex.hpp>
#include <mathlib\core\Matrix.hpp>
#include <mathlib\core\Vector.hpp>
#include <mathlib\numerical\Differentiation.hpp>
#include <mathlib\numerical\Integration.hpp>
#include <mathlib\numerical\Roots.hpp>
#include <optional>

int main() {
    /*******************************************************************************************
     * COMPLEX NUMBER
     * *****************************************************************************************
     */
    constexpr mathlib::core::Complexd cplx_num1(1.5, 2.5);
    constexpr mathlib::core::Complexd cplx_num2(2.0, -0.8);

    std::cout << "Number 1: " << cplx_num1 << "\n";
    std::cout << "Number 2: " << cplx_num2 << "\n";

    /*******************************************************************************************
     * VECTOR
     * *****************************************************************************************
     */

    // default initialisation
    mathlib::core::Vector<float, 4> vec_1{};
    std::cout << "Vector 1 dimension: " << vec_1.get_dimension() << "\n";
    std::cout << "vector 1 --> " << vec_1 << "\n\n";

    // direct initialisation
    mathlib::core::Vector<float, 4> vec_2(3.0f);
    std::cout << "Vector 2 dimension: " << vec_2.get_dimension() << "\n";
    std::cout << "vector 2 --> " << vec_2 << "\n\n";

    // list initialisation
    mathlib::core::Vector<float, 4> vec_3{2, 3, 4, 5};
    std::cout << "Vector 3 dimension: " << vec_3.get_dimension() << "\n";
    std::cout << "vector 3 --> " << vec_3 << "\n\n";

    // constant vector initialisation
    mathlib::core::Vector<float, 4> vec_4_const{1, 2, 3, 4};
    std::cout << "Vector 4 dimension: " << vec_4_const.get_dimension() << "\n";
    std::cout << "vector 4 --> " << vec_4_const << "\n";
    std::cout << "Vector 4 element 2: " << vec_4_const[1] << "\n\n";

    // /*******************************************************************************************
    //  * MATRIX
    //  * *****************************************************************************************
    //  */

    // create zero matrix
    mathlib::core::Matrix<double, 3, 3> mat_1{};
    std::cout << "--- Matrix 1 ---\n";
    std::cout << mat_1 << "\n";

    std::cout << "Mat 1 rows: " << mat_1.get_rows() << "\n";
    std::cout << "Mat 1 cols: " << mat_1.get_cols() << "\n";
    std::cout << "Mat 1 size: " << mat_1.get_size() << "\n\n";

    // create matrix with value
    mathlib::core::Matrix<double, 3, 3> mat_2(6);
    std::cout << "--- Matrix 2 ---\n";
    std::cout << mat_2 << "\n";

    // create matrix with initializer list
    mathlib::core::Matrix<double, 3, 3> mat_3{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    std::cout << "--- Matrix 3 ---\n";
    std::cout << mat_3 << "\n";

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

    /*******************************************************************************************
     * INTEGRATION
     * *****************************************************************************************
     */

    std::cout << "Regular polynomial integration\n\n";
    solution = mathlib::numerical::trapezoidal_integral(poly_func, 2.0, 5.0, 1000);
    std::cout << "Trapezoidal integration solution: " << std::setprecision(10) << solution << "\n";

    solution = mathlib::numerical::simpson_one_three(poly_func, 2.0, 5.0, 1000);
    std::cout << "Simpson's 1/3 rule solution: " << std::setprecision(10) << solution << "\n";

    solution = mathlib::numerical::simpson_three_eight(poly_func, 2.0, 5.0, 1000);
    std::cout << "Simpson's 3/8 rule solution: " << std::setprecision(10) << solution << "\n\n";

    // /*******************************************************************************************
    //  * ROOT FINDING
    //  * *****************************************************************************************
    //  */

    auto poly_func_2 = [](double x) { return x * x + (static_cast<double>(2.0) * x) - static_cast<double>(15.0); };

    auto first_deriv_2 = [](double x) { return static_cast<double>(2.0) * x + 2; };

    std::cout << "Regular polynomial integration\n\n";
    std::optional<double> solu = mathlib::numerical::Newton_Rhapson(poly_func_2, first_deriv_2, 6.0, 200);
    if (solution) {
        std::cout << "Newton_Rhapson solution: " << std::setprecision(10) << *solu << "\n";
    } else {
        std::cout << "Newton_Rhapson solution:: No roots found\n";
    }

    return 0;
}
