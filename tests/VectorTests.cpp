#include <gtest/gtest.h>

#include <mathlib/core/Vector.hpp>

// --- Construction Tests ---
TEST(VectorTests, Construction) {
    // 1. Default initialisation (zero initialised)
    mathlib::core::Vector<float, 3> v1;
    EXPECT_FLOAT_EQ(v1[0], 0.0);
    EXPECT_FLOAT_EQ(v1[1], 0.0);
    EXPECT_FLOAT_EQ(v1[2], 0.0);

    // 2. Scalar initialisation
    mathlib::core::Vector<float, 3> v2(5.0);
    EXPECT_FLOAT_EQ(v2[0], 5.0);
    EXPECT_FLOAT_EQ(v2[1], 5.0);
    EXPECT_FLOAT_EQ(v2[2], 5.0);

    // 3. List initialisation
    mathlib::core::Vector<float, 3> v3 = {7.0, 8.0, 9.0};
    EXPECT_FLOAT_EQ(v3[0], 7.0);
    EXPECT_FLOAT_EQ(v3[1], 8.0);
    EXPECT_FLOAT_EQ(v3[2], 9.0);
}

// --- Arithmetic Tests ---
TEST(VectorTests, Arithmetic) {
    mathlib::core::Vector<float, 3> v1 = {1.0, 7.0, 3.0};
    mathlib::core::Vector<float, 3> v2 = {4.0, 12.0, -2.0};

    // Addition
    mathlib::core::Vector<float, 3> result = v1 + v2;
    EXPECT_FLOAT_EQ(result[0], 5.0);
    EXPECT_FLOAT_EQ(result[1], 19.0);
    EXPECT_FLOAT_EQ(result[2], 1.0);

    // subtraction
    result = v1 - v2;
    EXPECT_FLOAT_EQ(result[0], -3.0);
    EXPECT_FLOAT_EQ(result[1], -5.0);
    EXPECT_FLOAT_EQ(result[2], 5.0);

    // Scalar Multiplication
    result = v1 * 2.0;
    EXPECT_FLOAT_EQ(result[0], 2.0);
    EXPECT_FLOAT_EQ(result[1], 14.0);
    EXPECT_FLOAT_EQ(result[2], 6.0);

    // Scalar Division
    result = v1 / 2.0;
    EXPECT_FLOAT_EQ(result[0], 0.5);
    EXPECT_FLOAT_EQ(result[1], 3.5);
    EXPECT_FLOAT_EQ(result[2], 1.5);

    // element wise multiplication (hadamard)
    result = v1 * v2;
    EXPECT_FLOAT_EQ(result[0], 4.0);
    EXPECT_FLOAT_EQ(result[1], 84.0);
    EXPECT_FLOAT_EQ(result[2], -6.0);

    // element wise division (hadamard)
    result = v1 / v2;
    EXPECT_FLOAT_EQ(result[0], 0.25);
    EXPECT_FLOAT_EQ(result[1], 0.5833333f);
    EXPECT_FLOAT_EQ(result[2], -1.5);
}

// --- Advanced Vector operation Tests ---
TEST(VectorTests, VectorAlgebra) {
    mathlib::core::Vector<float, 4> v1(3.0f);
    mathlib::core::Vector<float, 4> v2{2, 3, 4, 5};

    // utility function test
    EXPECT_FLOAT_EQ(v1.magnitude(), 6.0);
    EXPECT_FLOAT_EQ(v1.dot(v2), 42);
    EXPECT_NEAR(v1.angle_between(v2), 0.309193f, 1e-4f);
    EXPECT_FLOAT_EQ(v1.distance_between(v2), 2.44949f);
}

// --- Vector Utilities Tests ---
TEST(VectorTests, Utilities) {
    mathlib::core::Vector<float, 4> v1{2, 3, 4, 5};

    // check dimension
    EXPECT_EQ(v1.get_dimension(), 4);

    // check unit vector - sqrt() so use expect near
    mathlib::core::Vector<float, 4> v_unit = v1.get_unit_vector();
    EXPECT_NEAR(v_unit[0], 0.272166f, 1e-5f);
    EXPECT_NEAR(v_unit[1], 0.408248f, 1e-5f);
    EXPECT_NEAR(v_unit[2], 0.544331f, 1e-5f);
    EXPECT_NEAR(v_unit[3], 0.680414f, 1e-5f);

    // vector projection
    mathlib::core::Vector<float, 4> v2(3.0f);
    mathlib::core::Vector<float, 4> v_proj = v2.projection(v1);
    EXPECT_NEAR(v_proj[0], 1.55556f, 1e-5f);
    EXPECT_NEAR(v_proj[1], 2.33333f, 1e-5f);
    EXPECT_NEAR(v_proj[2], 3.11111f, 1e-5f);
    EXPECT_NEAR(v_proj[3], 3.88889f, 1e-5f);
}

// --- Exception Tests ---
TEST(VectorTests, SafetyChecks) {
    mathlib::core::Vector<float, 4> v1{2, 3, 4, 5};

    // scalar divison by zero
    EXPECT_THROW(v1 / 0.0, std::overflow_error);

    // vector division by zero
    mathlib::core::Vector<float, 4> vec_zero{};
    EXPECT_THROW(v1 / vec_zero, std::overflow_error);
}
