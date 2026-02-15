#include <cstddef>
#include <iomanip>
#include <iostream>
#include <mathlib\Utilities.hpp>
#include <mathlib\core\Vector.hpp>

int main() {
    /*******************************************************************************************
     * VECTOR DEMO
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

    // negation of vector
    std::cout << "Negation of vector 2 --> " << -vec_2 << "\n\n";

    // mathematical operations
    std::cout << "--- Mathematical operations of vectors ---\n";
    std::cout << "vector 2 --> " << vec_2 << " ||| " << "vector 3 --> " << vec_3 << "\n\n";

    std::cout << "Vector 2 + Vector 3 --> " << vec_2 + vec_3 << "\n";
    std::cout << "Vector 2 - Vector 3 --> " << vec_2 - vec_3 << "\n";
    std::cout << "Vector 2 * Vector 3 --> " << vec_2 * vec_3 << "\n";
    std::cout << "Vector 2 / Vector 3 --> " << vec_2 / vec_3 << "\n";
    std::cout << "Vector 2 * -6 --> " << vec_2 * -6 << "\n";
    std::cout << "Vector 2 / -6 --> " << vec_2 / -6 << "\n\n";

    // vector operations
    std::cout << "--- Vector operations ---\n";
    std::cout << "Vector 2 magnitude --> " << vec_2.magnitude() << "\n";
    std::cout << "Vector 2 unit vector --> " << vec_2.get_unit_vector() << "\n";
    std::cout << "Vector 3 magnitude --> " << vec_3.magnitude() << "\n";
    std::cout << "Vector 3 unit vector --> " << vec_3.get_unit_vector() << "\n";
    std::cout << "Vector 2 . Vector 3 --> " << vec_2.dot(vec_3) << "\n";
    std::cout << "Angle between Vector 2 & Vector 3 --> " << vec_2.angle_between(vec_3) << " radians \n";
    std::cout << "Distance between Vector 2 & Vector 3 --> " << vec_2.distance_between(vec_3) << "\n";
    std::cout << "Projection of vector 2 onto vector 3 --> " << vec_2.projection(vec_3) << " \n";

    return 0;
}
