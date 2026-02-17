include(FetchContent)

message(STATUS "Fetching GoogleTest...")

# 1. Download GoogleTest from GitHub
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
  # [FIX] Force extracted files to have the current timestamp.
  # This suppresses the CMP0135 warning and ensures correct rebuilding.
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

# 2. Prevent GoogleTest from overwriting compiler settings (Critical for Windows)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# 3. Make it available
FetchContent_MakeAvailable(googletest)

# --- Define the Test Executable ---

# Create one big executable that runs all tests
add_executable(MathLibTests
    VectorTests.cpp
    ComplexTests.cpp
    MatrixTests.cpp
    DifferentiationTests.cpp
    IntegrationTests.cpp
)

# Link MathLib (code to test) and GTest (testing framework)
target_link_libraries(MathLibTests
    PRIVATE
    MathLib
    GTest::gtest_main
)

# Register with CTest (CMake's test runner)
include(GoogleTest)
gtest_discover_tests(MathLibTests)
