#include <gtest/gtest.h>

#include <mathlib\core\Complex.hpp>

TEST(ComplexTests, Construction) {
    // 1. Default initialisation
    mathlib::core::Complex<double> cplx_1;
    EXPECT_DOUBLE_EQ(cplx_1.m_real, 0.0);
    EXPECT_DOUBLE_EQ(cplx_1.m_imag, 0.0);

    // 2. Direct initialisation
    mathlib::core::Complex<double> cplx_2(2.0);
    EXPECT_DOUBLE_EQ(cplx_2.m_real, 2.0);
    EXPECT_DOUBLE_EQ(cplx_2.m_imag, 0.0);

    // 3. Direct initialisation (2 arguements)
    mathlib::core::Complex<double> cplx_3(2.0, -1.0);
    EXPECT_DOUBLE_EQ(cplx_3.m_real, 2.0);
    EXPECT_DOUBLE_EQ(cplx_3.m_imag, -1.0);
}

TEST(ComplexTests, Arithmetic) {
    mathlib::core::Complexd cplx_num1(1.5, 2.5);
    mathlib::core::Complexd cplx_num2(2.0, -0.8);

    // Addition
    mathlib::core::Complexd result = cplx_num1 + cplx_num2;
    EXPECT_DOUBLE_EQ(result.m_real, 3.5);
    EXPECT_DOUBLE_EQ(result.m_imag, 1.7);

    // Subtraction
    result = cplx_num1 - cplx_num2;
    EXPECT_DOUBLE_EQ(result.m_real, -0.5);
    EXPECT_DOUBLE_EQ(result.m_imag, 3.3);

    // multiplication
    result = cplx_num1 * cplx_num2;
    EXPECT_DOUBLE_EQ(result.m_real, 5);
    EXPECT_DOUBLE_EQ(result.m_imag, 3.8);

    // division
    result = cplx_num1 / cplx_num2;
    EXPECT_NEAR(result.m_real, 0.215517, 1e-4);
    EXPECT_NEAR(result.m_imag, 1.33621, 1e-4);

    // multiplication (scalar)
    result = cplx_num1 * 6.0;
    EXPECT_DOUBLE_EQ(result.m_real, 9.0);
    EXPECT_DOUBLE_EQ(result.m_imag, 15.0);

    // division (scalar)
    result = cplx_num1 / 6.0;
    EXPECT_NEAR(result.m_real, 0.25, 1e-4);
    EXPECT_NEAR(result.m_imag, 0.416667, 1e-4);
}

TEST(ComplexTests, Utilities) {
    mathlib::core::Complexd cplx_num1(1.5, 2.5);

    // conjugate check
    mathlib::core::Complexd result = cplx_num1.conjugate();
    EXPECT_DOUBLE_EQ(result.m_real, 1.5);
    EXPECT_DOUBLE_EQ(result.m_imag, -2.5);

    // magnitude
    EXPECT_NEAR(cplx_num1.magnitude(), 2.91548, 1e-4);

    // angle
    EXPECT_NEAR(cplx_num1.angle(), 1.03038, 1e-4);
}

TEST(ComplexTests, SafetyChecks) {
    mathlib::core::Complexd cplx_num1(1.5, 2.5);
    mathlib::core::Complexd cplx_num2(2.0, -0.8);
    mathlib::core::Complexd cplx_num3{};

    // check division by zero
    EXPECT_THROW(cplx_num1 / 0.0, std::overflow_error);
    EXPECT_THROW(cplx_num1 / cplx_num3, std::overflow_error);

    // no throws to check for safe divison
    EXPECT_NO_THROW(cplx_num1 / 2.0);
    EXPECT_NO_THROW(cplx_num1 / cplx_num2);
}
