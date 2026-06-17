/**
 * @file test_nonlinear.cpp
 * @brief Testy jednostkowe dla metod rozwiazywania rownan nieliniowych.
 */
#include "../include/numlib.h"
#include <iostream>
#include <cmath>
#include <cassert>

static bool close(long double a, long double b, long double tol = 1e-10L) {
    return fabsl(a - b) < tol;
}

// f(x) = x^2 - 2,  pierwiastek: sqrt(2) ~ 1.41421356...
static long double f_sqrt2(long double x)  { return x * x - 2.0L; }
static long double df_sqrt2(long double x) { return 2.0L * x; }

// f(x) = cos(x) - x,  pierwiastek: ~ 0.7390851332151607
static long double f_cosx(long double x)  { return cosl(x) - x; }
static long double df_cosx(long double x) { return -sinl(x) - 1.0L; }

// f(x) = x^3 - x - 2,  pierwiastek: ~ 1.5213797...
static long double f_cubic(long double x)  { return x * x * x - x - 2.0L; }
static long double df_cubic(long double x) { return 3.0L * x * x - 1.0L; }

// Test 1: Newton analityczny - sqrt(2)
void test_newton_analytic_sqrt2() {
    long double r = numlib::newton_analytic(f_sqrt2, df_sqrt2, 1.5L);
    assert(close(r, sqrtl(2.0L)));
    std::cout << "[PASS] test_newton_analytic_sqrt2 (r=" << r << ")\n";
}

// Test 2: Newton numeryczny - sqrt(2)
void test_newton_numeric_sqrt2() {
    long double r = numlib::newton_numeric(f_sqrt2, 1.5L);
    assert(close(r, sqrtl(2.0L)));
    std::cout << "[PASS] test_newton_numeric_sqrt2 (r=" << r << ")\n";
}

// Test 3: Sieczne - sqrt(2), dwa przyblizenia poczatkowe
void test_secant_sqrt2() {
    long double r = numlib::secant(f_sqrt2, 1.0L, 2.0L);
    assert(close(r, sqrtl(2.0L)));
    std::cout << "[PASS] test_secant_sqrt2 (r=" << r << ")\n";
}

// Test 4: Newton analityczny - cos(x) = x
void test_newton_cosx() {
    long double r = numlib::newton_analytic(f_cosx, df_cosx, 0.5L);
    assert(close(r, 0.7390851332151607L));
    std::cout << "[PASS] test_newton_cosx (r=" << r << ")\n";
}

// Test 5: Sieczne - wielomian trzeciego stopnia
void test_secant_cubic() {
    long double r = numlib::secant(f_cubic, 1.0L, 2.0L);
    assert(fabsl(f_cubic(r)) < 1e-10L);
    std::cout << "[PASS] test_secant_cubic (r=" << r << ", f(r)=" << f_cubic(r) << ")\n";
}

// Test 6: Newton analityczny i numeryczny daja zblizone wyniki
void test_newton_analytic_vs_numeric() {
    long double rA = numlib::newton_analytic(f_cubic, df_cubic, 1.5L);
    long double rN = numlib::newton_numeric(f_cubic, 1.5L);
    assert(close(rA, rN, 1e-10L));
    std::cout << "[PASS] test_newton_analytic_vs_numeric\n";
}

// Test 7: Pochodna numeryczna - sprawdzenie dokladnosci
void test_numerical_derivative() {
    // d/dx (x^3) w x=2 = 12
    long double d = numlib::numerical_derivative(
        [](long double x){ return x * x * x; }, 2.0L);
    assert(close(d, 12.0L, 1e-7L));
    std::cout << "[PASS] test_numerical_derivative (d=" << d << ")\n";
}

int main() {
    std::cout << "Testy: rownania nieliniowe\n";
    test_newton_analytic_sqrt2();
    test_newton_numeric_sqrt2();
    test_secant_sqrt2();
    test_newton_cosx();
    test_secant_cubic();
    test_newton_analytic_vs_numeric();
    test_numerical_derivative();
    std::cout << "Wszystkie testy rownan nieliniowych: OK\n\n";
    return 0;
}
