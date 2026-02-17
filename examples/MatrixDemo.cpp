#include <cstddef>
#include <iomanip>
#include <iostream>
#include <mathlib/Utilities.hpp>
#include <mathlib/core/Matrix.hpp>
#include <mathlib/core/Vector.hpp>

int main() {
    /*******************************************************************************************
     * MATRIX
     * *****************************************************************************************
     */

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

    // accessing using ()
    for (std::size_t i = 0; i < mat_2.m_rows; ++i) {
        mat_2(i, i) = 1;
    }
    std::cout << "--- Matrix 2 ---\n";
    std::cout << mat_2 << "\n";

    // set entire row
    mathlib::core::Vector<double, 3> new_row{4, 4, 4};
    mat_2.set_row(0, new_row);
    std::cout << "--- Matrix 2 ---\n";
    std::cout << mat_2 << "\n";

    // set entire column
    mathlib::core::Vector<double, 3> new_col{7, 7, 7};
    mat_2.set_column(1, new_col);
    std::cout << "--- Matrix 2 ---\n";
    std::cout << mat_2 << "\n";

    // retrieve column and verify it
    auto column_1 = mat_2.get_column_vector(1);
    if (column_1 == new_col) {
        std::cout << column_1 << "\n";
    }

    // retrieve row
    auto row_1 = mat_2.get_row_vector(1);
    std::cout << row_1 << "\n\n";

    // get identity
    std::cout << "Matrix 2 identity \n";
    std::cout << mat_2.get_identity() << "\n";

    // mathematical operations
    std::cout << "Matrix operations \n\n";
    std::cout << "--- Matrix 2 ---\n";
    std::cout << mat_2 << "\n";
    std::cout << "--- Matrix 3 ---\n";
    std::cout << mat_3 << "\n";

    std::cout << "Mat 2 + Mat 3 \n";
    std::cout << mat_2 + mat_3 << "\n";

    std::cout << "Mat 2 - Mat 3 \n";
    std::cout << mat_2 - mat_3 << "\n";

    std::cout << "Mat 2 * 2 \n";
    std::cout << mat_2 * 2 << "\n";

    std::cout << "Mat 2 / 2 \n";
    std::cout << mat_2 / 2 << "\n";

    std::cout << "Mat 2 * Mat 3 \n";
    std::cout << mat_2 * mat_3 << "\n";

    // special matrix operations
    std::cout << "Mat 2 * Mat 3 (hadamard Product)\n";
    std::cout << mat_2.Hadamard_product(mat_3) << "\n";

    std::cout << "Mat 2 traspose\n";
    std::cout << mat_2.transpose() << "\n";

    std::cout << "Mat 2 reshape\n";
    std::cout << mat_2.reshape<1, 9>() << "\n";

    std::cout << "Mat 2 submatrix\n";
    std::cout << mat_2.sub_matrix<2, 2>(0, 0) << "\n";
    return 0;
}
