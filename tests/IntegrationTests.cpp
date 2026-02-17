#include <gtest/gtest.h>

#include <mathlib/numerical/Integration.hpp>

TEST(IntegrationTests, Trapezoidal) {
    auto poly_func = [](double x) {
        return std::pow(x, 4) + (static_cast<double>(3.0) * std::pow(x, 3)) + static_cast<double>(5.0) * x;
    };

    double solution = mathlib::numerical::trapezoidal_integral(poly_func, 2.0, 5.0, 1000);
    EXPECT_NEAR(solution, 1127.850493, 1e-5);
}

TEST(IntegrationTests, Simpson_1_3) {
    auto poly_func = [](double x) {
        return std::pow(x, 4) + (static_cast<double>(3.0) * std::pow(x, 3)) + static_cast<double>(5.0) * x;
    };

    double solution = mathlib::numerical::simpson_one_three(poly_func, 2.0, 5.0, 1000);
    EXPECT_NEAR(solution, 1126.775985, 1e-5);
}

TEST(IntegrationTests, Simpson_3_8) {
    auto poly_func = [](double x) {
        return std::pow(x, 4) + (static_cast<double>(3.0) * std::pow(x, 3)) + static_cast<double>(5.0) * x;
    };

    double solution = mathlib::numerical::simpson_three_eight(poly_func, 2.0, 5.0, 1000);
    EXPECT_NEAR(solution, 1127.85, 1e-5);
}
