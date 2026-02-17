## MathLib

A modern, high-performance, header-only C++20 Math Library designed for numerical analysis and scientific computation.

### 🚀 Features

Module: Core

Vectors: N-dimensional templated vectors.

Matrices: Templated $R \times C$ matrices with support for Row-Major order storage.

Complex Numbers: Full support for complex arithmetic ($a + bi$).

Module: Numerical

Root Finding: Newton-Raphson, Secant, Bisection, and False Position methods.

Integration: Numerical approximation via Trapezoidal and Simpson's rules.

Differentiation: Central difference and Richardson Extrapolation for first and second derivatives.

### 🛠️ Requirements

Compiler: C++20 compliant compiler (GCC 10+, Clang 10+, or MSVC 19.29+).

Build System: CMake 3.20 or higher.

### 📦 Installation

MathLib is header-only. You can simply copy the include/ directory to your project, or integrate it via CMake.

Integration via CMake

Add this to your CMakeLists.txt:

```
add_subdirectory(path/to/MathLib)
target_link_libraries(YourProject PRIVATE MathLib)
```

Or using FetchContent:

```
include(FetchContent)
FetchContent_Declare(
  mathlib
  GIT_REPOSITORY [https://github.com/RVThan06/MathLib](https://github.com/RVThan06/MathLib)
  GIT_TAG main
)
FetchContent_MakeAvailable(mathlib)
```

### 💻 Quick Start

Vector operations

```
#include <mathlib/core/Vector.hpp>

void example() {
    mathlib::core::Vector<float, 3> vec_1{2, 3, 4};
    mathlib::core::Vector<float, 3> vec_2{1, 1, 1};

    // Dot product
    double dot_value = vec_1.dot(vec_2);

    // Angle between two vector
    double angle =  vec_1.angle_between(vec_2)

    // Addition of vector
    mathlib::core::Vector<float, 3> temp =  vec_1 + vec_2;
}
```

Numerical (Root Finding)

```
#include <mathlib/numerical/Roots.hpp>

auto func = [](double x) { return x * x - 2.0; };
auto deriv = [](double x) { return 2.0 * x; };

auto root = mathlib::numerical::Newton_Raphson(func, deriv, 1.0);
if (root) {
    std::cout << "Square root of 2 is approx: " << *root << std::endl;
}
```

### 🧪 Running Tests

This project uses GoogleTest for unit testing. To build and run the suite:

```
cmake -B build -DMATHLIB_BUILD_TESTS=ON
cmake --build build
cd build && ctest --output-on-failure
```
