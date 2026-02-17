#include <gtest/gtest.h>

#include <mathlib/numerical/Differentiation.hpp>

TEST(DifferentiationTests, CentralDifference) {
    auto poly_func = [](double x) {
        return std::pow(x, 4) + (static_cast<double>(3.0) * std::pow(x, 3)) + static_cast<double>(5.0) * x;
    };

    double solution = mathlib::numerical::central_difference(poly_func, 1.5);
    EXPECT_NEAR(solution, 38.7499999, 1e-5);
}

TEST(DifferentiationTests, RichardsonExtra) {
    auto poly_func = [](double x) {
        return std::pow(x, 4) + (static_cast<double>(3.0) * std::pow(x, 3)) + static_cast<double>(5.0) * x;
    };

    double solution = mathlib::numerical::richardson_extrapolation(poly_func, 1.5);
    EXPECT_NEAR(solution, 38.75, 1e-5);
}

TEST(DifferentiationTests, SecondDerivative) {
    auto poly_func = [](double x) {
        return std::pow(x, 4) + (static_cast<double>(3.0) * std::pow(x, 3)) + static_cast<double>(5.0) * x;
    };

    double solution = mathlib::numerical::second_derivative(poly_func, 1.5);
    EXPECT_NEAR(solution, 54.00124792, 1e-5);
}
