#include <gtest/gtest.h>

#include <mathlib\numerical\Roots.hpp>

TEST(RootTests, NewtonRhapson) {
    auto poly_func_2 = [](double x) { return x * x + (static_cast<double>(2.0) * x) - static_cast<double>(15.0); };
    auto first_deriv_2 = [](double x) { return static_cast<double>(2.0) * x + 2; };

    // Newton Rhapson method
    std::optional<double> solution = mathlib::numerical::Newton_Rhapson(poly_func_2, first_deriv_2, 6.0, 200);
    EXPECT_NEAR(*solution, 3.0, 1e-3);
}

TEST(RootTests, Bisection) {
    auto poly_func_2 = [](double x) { return x * x + (static_cast<double>(2.0) * x) - static_cast<double>(15.0); };

    // Bisection method
    std::optional<double> solution = mathlib::numerical::Bisection(-1.0, 6.0, poly_func_2);
    EXPECT_NEAR(*solution, 3.0, 1e-3);
}

TEST(RootTests, Secant) {
    auto poly_func_2 = [](double x) { return x * x + (static_cast<double>(2.0) * x) - static_cast<double>(15.0); };

    // Secant method
    std::optional<double> solution = mathlib::numerical::Secant(poly_func_2, 6.0, -1.0, 200);
    EXPECT_NEAR(*solution, 3.0, 1e-3);
}

TEST(RootTests, FalsePosition) {
    auto poly_func_2 = [](double x) { return x * x + (static_cast<double>(2.0) * x) - static_cast<double>(15.0); };

    // False position method
    std::optional<double> solution = mathlib::numerical::False_Position(-1.0, 6.0, poly_func_2);
    EXPECT_NEAR(*solution, 3.0, 1e-3);
}
