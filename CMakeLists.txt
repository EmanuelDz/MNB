cmake_minimum_required(VERSION 3.10)
project(numlib VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Flagi kompilacji
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -O2)
endif()

# Biblioteka statyczna numlib
add_library(numlib STATIC
    src/linear_systems.cpp
    src/interpolation.cpp
    src/approximation.cpp
    src/integration.cpp
    src/ode.cpp
    src/nonlinear.cpp
)
target_include_directories(numlib PUBLIC include)
if(NOT MSVC)
    target_link_libraries(numlib m)
endif()

# Testy
enable_testing()

foreach(test_name
        test_linear
        test_interpolation
        test_integration
        test_ode
        test_nonlinear
        test_approximation)
    add_executable(${test_name} tests/${test_name}.cpp)
    target_link_libraries(${test_name} numlib)
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()

# Przyklady
add_executable(example1 examples/example_integration_nonlinear.cpp)
target_link_libraries(example1 numlib)

add_executable(example2 examples/example_ode_interpolation.cpp)
target_link_libraries(example2 numlib)
