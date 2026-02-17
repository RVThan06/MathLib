#include <gtest/gtest.h>

#include <mathlib\core\Matrix.hpp>

TEST(MatrixTests, Construction) {
    // 1. Default intialisation
    mathlib::core::Matrix<double, 2, 2> mat_1{};
    EXPECT_DOUBLE_EQ(mat_1(0, 0), 0);
    EXPECT_DOUBLE_EQ(mat_1(0, 1), 0);
    EXPECT_DOUBLE_EQ(mat_1(1, 0), 0);
    EXPECT_DOUBLE_EQ(mat_1(1, 1), 0);

    // 1. Direct intialisation
    mathlib::core::Matrix<double, 2, 2> mat_2(1);
    EXPECT_DOUBLE_EQ(mat_2(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(mat_2(0, 1), 1.0);
    EXPECT_DOUBLE_EQ(mat_2(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(mat_2(1, 1), 1.0);

    // 1. List intialisation
    mathlib::core::Matrix<double, 2, 2> mat_3{{1, 2}, {3, 4}};
    EXPECT_DOUBLE_EQ(mat_3(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(mat_3(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(mat_3(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(mat_3(1, 1), 4.0);
}

TEST(MatrixTests, Arithmetic) {
    mathlib::core::Matrix<double, 2, 2> mat_1(1);
    mathlib::core::Matrix<double, 2, 2> mat_2(2);

    // matrix addition
    mathlib::core::Matrix<double, 2, 2> result = mat_1 + mat_2;
    EXPECT_DOUBLE_EQ(result(0, 0), 3.0);
    EXPECT_DOUBLE_EQ(result(0, 1), 3.0);
    EXPECT_DOUBLE_EQ(result(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(result(1, 1), 3.0);

    // matrix subtraction
    result = mat_1 - mat_2;
    EXPECT_DOUBLE_EQ(result(0, 0), -1.0);
    EXPECT_DOUBLE_EQ(result(0, 1), -1.0);
    EXPECT_DOUBLE_EQ(result(1, 0), -1.0);
    EXPECT_DOUBLE_EQ(result(1, 1), -1.0);

    // matrix multiplication
    result = mat_1 * mat_2;
    EXPECT_DOUBLE_EQ(result(0, 0), 4.0);
    EXPECT_DOUBLE_EQ(result(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(result(1, 0), 4.0);
    EXPECT_DOUBLE_EQ(result(1, 1), 4.0);

    // matrix multiplication (scalar)
    result = mat_1 * 4;
    EXPECT_DOUBLE_EQ(result(0, 0), 4.0);
    EXPECT_DOUBLE_EQ(result(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(result(1, 0), 4.0);
    EXPECT_DOUBLE_EQ(result(1, 1), 4.0);

    // matrix division scalar
    result = mat_1 / 2.0;
    EXPECT_DOUBLE_EQ(result(0, 0), 0.5);
    EXPECT_DOUBLE_EQ(result(0, 1), 0.5);
    EXPECT_DOUBLE_EQ(result(1, 0), 0.5);
    EXPECT_DOUBLE_EQ(result(1, 1), 0.5);
}

TEST(MatrixTests, Operations) {
    mathlib::core::Matrix<double, 2, 2> mat_1(1);
    mathlib::core::Matrix<double, 2, 2> mat_2(2);

    // Hadamard Product
    mathlib::core::Matrix<double, 2, 2> result = mat_1.Hadamard_product(mat_2);
    EXPECT_DOUBLE_EQ(result(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(result(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(result(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(result(1, 1), 2.0);

    // transpose
    mathlib::core::Matrix<double, 2, 2> mat_3{{1, 2}, {3, 4}};
    result = mat_3.transpose();
    EXPECT_DOUBLE_EQ(result(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(result(0, 1), 3.0);
    EXPECT_DOUBLE_EQ(result(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(result(1, 1), 4.0);

    // reshape
    auto temp = mat_3.reshape<1, 4>();
    EXPECT_EQ(temp.m_rows, 1);
    EXPECT_EQ(temp.m_cols, 4);

    // submatrix
    EXPECT_EQ((mat_3.sub_matrix<1, 2>(0, 0)).m_rows, 1);
    EXPECT_EQ((mat_3.sub_matrix<1, 2>(0, 0)).m_cols, 2);
}

TEST(MatrixTests, SafetyCheck) {
    mathlib::core::Matrix<double, 2, 2> mat_1(1);

    // zero division check
    EXPECT_THROW(mat_1 / 0.0, std::overflow_error);

    // jagged intialiser list check
    EXPECT_THROW((mathlib::core::Matrix<double, 2, 2>{{1, 2}, {2}}), std::length_error);
}
