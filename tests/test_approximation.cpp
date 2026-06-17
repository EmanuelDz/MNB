/**
 * @file test_approximation.cpp
 * @brief Testy jednostkowe dla aproksymacji sredniokwadratowej.
 */
#include "../include/numlib.h"
#include <iostream>
#include <cmath>
#include <cassert>

static bool close(long double a, long double b, long double tol) {
    return fabsl(a - b) < tol;
}

static long double f_linear(long double x) { return 2.0L * x + 1.0L; }
static long double f_quadratic(long double x) { return x * x - x + 3.0L; }
static long double f_sin(long double x) { return sinl(x); }
static long double f_const(long double x) { (void)x; return 5.0L; }

// Test 1: Aproksymacja wielomianem st. 1 funkcji liniowej - dokladna rekonstrukcja
void test_approx_linear() {
    auto c = numlib::polynomial_approximation(f_linear, 0.0L, 2.0L, 1);
    assert(!c.empty());
    // c[0] ~ 1, c[1] ~ 2
    assert(close(c[0], 1.0L, 1e-6L));
    assert(close(c[1], 2.0L, 1e-6L));
    std::cout << "[PASS] test_approx_linear (a0=" << c[0] << " a1=" << c[1] << ")\n";
}

// Test 2: Aproksymacja wielomianem st. 2 funkcji kwadratowej - dokladna rekonstrukcja
void test_approx_quadratic() {
    // x^2 - x + 3 na [0, 4]
    auto c = numlib::polynomial_approximation(f_quadratic, 0.0L, 4.0L, 2);
    assert(!c.empty());
    assert(close(c[0], 3.0L, 1e-5L));
    assert(close(c[1], -1.0L, 1e-5L));
    assert(close(c[2], 1.0L, 1e-5L));
    std::cout << "[PASS] test_approx_quadratic (a0=" << c[0] << " a1=" << c[1] << " a2=" << c[2] << ")\n";
}

// Test 3: Aproksymacja funkcji stalej wielomianem st. 0
void test_approx_constant() {
    auto c = numlib::polynomial_approximation(f_const, 0.0L, 1.0L, 0);
    assert(!c.empty());
    assert(close(c[0], 5.0L, 1e-8L));
    std::cout << "[PASS] test_approx_constant (a0=" << c[0] << ")\n";
}

// Test 4: Aproksymacja sin(x) wielomianem st. 5 — blad max < 0.01 na [0, pi]
void test_approx_sin_degree5() {
    const long double PI = 3.14159265358979323846L;
    auto c = numlib::polynomial_approximation(f_sin, 0.0L, PI, 5);
    assert(!c.empty());
    long double maxErr = 0.0L;
    int pts = 100;
    for (int i = 0; i <= pts; i++) {
        long double x = PI * i / pts;
        long double err = fabsl(sinl(x) - numlib::approx_evaluate(c, x));
        if (err > maxErr) maxErr = err;
    }
    assert(maxErr < 0.01L);
    std::cout << "[PASS] test_approx_sin_degree5 (max blad=" << maxErr << ")\n";
}

// Test 5: Wiekszy stopien => mniejszy blad aproksymacji sin(x)
void test_approx_higher_degree_less_error() {
    const long double PI = 3.14159265358979323846L;
    auto c3 = numlib::polynomial_approximation(f_sin, 0.0L, PI, 3);
    auto c7 = numlib::polynomial_approximation(f_sin, 0.0L, PI, 7);
    int pts = 50;
    long double err3 = 0.0L, err7 = 0.0L;
    for (int i = 0; i <= pts; i++) {
        long double x = PI * i / pts;
        err3 += fabsl(sinl(x) - numlib::approx_evaluate(c3, x));
        err7 += fabsl(sinl(x) - numlib::approx_evaluate(c7, x));
    }
    assert(err7 < err3);
    std::cout << "[PASS] test_approx_higher_degree_less_error" << " (st3=" << err3 << " st7=" << err7 << ")\n";
}

int main() {
    std::cout << "Testy: aproksymacja sredniokwadratowa\n";
    test_approx_linear();
    test_approx_quadratic();
    test_approx_constant();
    test_approx_sin_degree5();
    test_approx_higher_degree_less_error();
    std::cout << "Wszystkie testy aproksymacji: OK\n\n";
    return 0;
}
