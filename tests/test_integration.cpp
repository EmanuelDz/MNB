/**
 * @file test_integration.cpp
 * @brief Testy jednostkowe dla metod calkowania numerycznego.
 */
#include "../include/numlib.h"
#include <iostream>
#include <cmath>
#include <cassert>

static bool close(long double a, long double b, long double tol) {
    return fabsl(a - b) < tol;
}

// f(x) = x^2,  calka od 0 do 1 = 1/3
static long double f_square(long double x) { return x * x; }

// f(x) = sin(x),  calka od 0 do pi = 2
static long double f_sin(long double x) { return sinl(x); }

// f(x) = e^x,  calka od 0 do 1 = e - 1
static long double f_exp(long double x) { return expl(x); }

// f(x) = (1-x)*exp(sin(x)),  calka od -1 do 15 ~ -136.183666144466
static long double f_specific(long double x) { return (1.0L - x) * expl(sinl(x)); }

const long double PI = 3.14159265358979323846L;

// Test 1: Trapezy — calka x^2 od 0 do 1
void test_trapezoidal_square() {
    long double I = numlib::trapezoidal(f_square, 0.0L, 1.0L, 10000);
    assert(close(I, 1.0L / 3.0L, 1e-6L));
    std::cout << "[PASS] test_trapezoidal_square\n";
}

// Test 2: Trapezy — calka sin(x) od 0 do pi
void test_trapezoidal_sin() {
    long double I = numlib::trapezoidal(f_sin, 0.0L, PI, 10000);
    assert(close(I, 2.0L, 1e-6L));
    std::cout << "[PASS] test_trapezoidal_sin\n";
}

// Test 3: Simpson - calka x^2 od 0 do 1 (dokladny dla wielomianow do st. 3)
void test_simpson_square() {
    long double I = numlib::simpson(f_square, 0.0L, 1.0L, 100);
    assert(close(I, 1.0L / 3.0L, 1e-12L));
    std::cout << "[PASS] test_simpson_square\n";
}

// Test 4: Simpson - calka sin(x) od 0 do pi
void test_simpson_sin() {
    long double I = numlib::simpson(f_sin, 0.0L, PI, 1000);
    assert(close(I, 2.0L, 1e-12L));
    std::cout << "[PASS] test_simpson_sin\n";
}

// Test 5: Gauss 4-wezlowy - calka e^x od 0 do 1, dokladny wynik: e-1
void test_gauss4_exp() {
    long double I = numlib::gauss_composite(f_exp, 0.0L, 1.0L, 4, 20);
    long double exact = expl(1.0L) - 1.0L;
    assert(close(I, exact, 1e-12L));
    std::cout << "[PASS] test_gauss4_exp\n";
}

// Test 6: Gauss - funkcja (1-x)*exp(sin(x)) od -1 do 15
void test_gauss_specific_function() {
    long double I = numlib::gauss_composite(f_specific, -1.0L, 15.0L, 4, 50);
    long double exact = -136.183666144466L;
    assert(close(I, exact, 1e-4L));
    std::cout << "[PASS] test_gauss_specific_function\n";
}

// Test 7: Simpson dokladniejszy niz trapezy dla tej samej liczby krokow
void test_simpson_better_than_trapezoidal() {
    long double exact = 2.0L; // calka sin od 0 do pi
    long double errTrap = fabsl(numlib::trapezoidal(f_sin, 0.0L, PI, 100) - exact);
    long double errSimp = fabsl(numlib::simpson(f_sin, 0.0L, PI, 100) - exact);
    assert(errSimp < errTrap);
    std::cout << "[PASS] test_simpson_better_than_trapezoidal\n";
}

int main() {
    std::cout << "Testy: calkowanie numeryczne\n";
    test_trapezoidal_square();
    test_trapezoidal_sin();
    test_simpson_square();
    test_simpson_sin();
    test_gauss4_exp();
    test_gauss_specific_function();
    test_simpson_better_than_trapezoidal();
    std::cout << "Wszystkie testy calkowania: OK\n\n";
    return 0;
}
