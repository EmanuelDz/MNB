/**
 * @file test_interpolation.cpp
 * @brief Testy jednostkowe dla metod interpolacji (Newton, Lagrange, Horner).
 */
#include "../include/numlib.h"
#include <iostream>
#include <cmath>
#include <cassert>

static bool close(long double a, long double b, long double tol = 1e-8L) {
    return fabsl(a - b) < tol;
}

// Test 1: Newton - wartosci w wezlach musza sie zgadzac
//   f(x) = x^2 + 1,  wezly: 0,1,2,3
void test_newton_in_nodes() {
    std::vector<long double> xi  = {0.0L, 1.0L, 2.0L, 3.0L};
    std::vector<long double> fxi = {1.0L, 2.0L, 5.0L, 10.0L};
    auto c = numlib::newton_coefficients(xi, fxi);
    for (int i = 0; i < 4; i++) {
        long double y = numlib::newton_evaluate(xi[i], xi, c);
        assert(close(y, fxi[i]));
    }
    std::cout << "[PASS] test_newton_in_nodes\n";
}

// Test 2: Newton — interpolacja miedzy wezlami
//   f(x) = x^2 + 1  =>  f(1.5) = 3.25
void test_newton_between_nodes() {
    std::vector<long double> xi  = {0.0L, 1.0L, 2.0L};
    std::vector<long double> fxi = {1.0L, 2.0L, 5.0L};
    auto c = numlib::newton_coefficients(xi, fxi);
    assert(close(numlib::newton_evaluate(1.5L, xi, c), 3.25L));
    std::cout << "[PASS] test_newton_between_nodes\n";
}

// Test 3: Lagrange - ten sam wynik co Newton
void test_lagrange_consistent_with_newton() {
    std::vector<long double> xi  = {0.0L, 1.0L, 2.0L};
    std::vector<long double> fxi = {1.0L, 2.0L, 5.0L};
    auto c = numlib::newton_coefficients(xi, fxi);
    long double yN = numlib::newton_evaluate(1.5L, xi, c);
    long double yL = numlib::lagrange_evaluate(1.5L, xi, fxi);
    assert(close(yN, yL, 1e-12L));
    std::cout << "[PASS] test_lagrange_consistent_with_newton\n";
}

// Test 4: Horner - poprawna wartosc wielomianu
//   p(x) = 1 + 2x + 3x^2  w x=2  =>  17
void test_horner_value() {
    std::vector<long double> a = {1.0L, 2.0L, 3.0L};
    assert(close(numlib::horner(a, 2.0L), 17.0L));
    std::cout << "[PASS] test_horner_value\n";
}

// Test 5: Horner vs metoda naturalna - identyczne wyniki
void test_horner_equals_natural() {
    std::vector<long double> a = {-3.0L, 0.0L, 5.0L, -1.0L, 2.0L};
    long double x = 1.7L;
    long double yH = numlib::horner(a, x);
    long double yN = numlib::polynomial_natural(a, x);
    assert(close(yH, yN, 1e-12L));
    std::cout << "[PASS] test_horner_equals_natural\n";
}

// Test 6: Lagrange dla funkcji sinus (5 wezlow)
void test_lagrange_sin() {
    const long double PI = 3.14159265358979323846L;
    std::vector<long double> xi  = {0.0L, PI/4, PI/2, 3*PI/4, PI};
    std::vector<long double> fxi;
    for (auto xi_val : xi) fxi.push_back(sinl(xi_val));
    // Interpolacja w srodku: sin(pi/2) = 1
    long double y = numlib::lagrange_evaluate(PI / 2.0L, xi, fxi);
    assert(close(y, 1.0L, 1e-8L));
    std::cout << "[PASS] test_lagrange_sin\n";
}

int main() {
    std::cout << "Testy: interpolacja\n";
    test_newton_in_nodes();
    test_newton_between_nodes();
    test_lagrange_consistent_with_newton();
    test_horner_value();
    test_horner_equals_natural();
    test_lagrange_sin();
    std::cout << "Wszystkie testy interpolacji: OK\n\n";
    return 0;
}
