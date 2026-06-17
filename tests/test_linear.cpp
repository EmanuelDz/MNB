/**
 * @file test_linear.cpp
 * @brief Testy jednostkowe dla metod rozwiazywania ukladow rownan liniowych (LU).
 */
#include "../include/numlib.h"
#include <iostream>
#include <cmath>
#include <cassert>

static bool close(long double a, long double b, long double tol = 1e-8L) {
    return fabsl(a - b) < tol;
}

// Test 1: Uklad 2x2
//   2x + y = 5
//   4x + 3y = 11   =>  x=2, y=1
void test_lu_2x2() {
    std::vector<long double> A = {2, 1, 4, 3};
    std::vector<long double> b = {5, 11};
    auto x = numlib::lu_solve(A, b, 2);
    assert(!x.empty());
    assert(close(x[0], 2.0L) && close(x[1], 1.0L));
    std::cout << "[PASS] test_lu_2x2\n";
}

// Test 2: Uklad 3x3
//   2x + y - z = 8
//  -3x - y + 2z = -11
//  -2x + y + 2z = -3   =>  x=2, y=3, z=-1
void test_lu_3x3() {
    std::vector<long double> A = { 2, 1,-1, -3,-1, 2, -2, 1, 2};
    std::vector<long double> b = {8, -11, -3};
    auto x = numlib::lu_solve(A, b, 3);
    assert(!x.empty());
    assert(close(x[0], 2.0L) && close(x[1], 3.0L) && close(x[2], -1.0L));
    std::cout << "[PASS] test_lu_3x3\n";
}

// Test 3: Macierz diagonalna  diag(1,2,3)*x = (3,6,9) => x=(3,3,3)
void test_lu_diagonal() {
    std::vector<long double> A = {1,0,0, 0,2,0, 0,0,3};
    std::vector<long double> b = {3, 6, 9};
    auto x = numlib::lu_solve(A, b, 3);
    assert(!x.empty());
    assert(close(x[0], 3.0L) && close(x[1], 3.0L) && close(x[2], 3.0L));
    std::cout << "[PASS] test_lu_diagonal\n";
}

// Test 4: Macierz osobliwa - oczekujemy pustego wektora
void test_lu_singular() {
    std::vector<long double> A = {0, 1, 0, 2};
    std::vector<long double> b = {1, 2};
    auto x = numlib::lu_solve(A, b, 2);
    assert(x.empty());
    std::cout << "[PASS] test_lu_singular\n";
}

// Test 5: Sprawdzenie A*x = b po rozwiazaniu (4x4)
void test_lu_verify_ax_eq_b() {
    std::vector<long double> A = { 4, 3, 2, 1, 1, 5, 3, 2, 2, 1, 6, 4, 3, 2, 1, 7};
    std::vector<long double> b = {10, 11, 13, 13};
    int n = 4;
    auto x = numlib::lu_solve(A, b, n);
    assert(!x.empty());
    // Weryfikacja A*x == b
    for (int i = 0; i < n; i++) {
        long double Ax_i = 0;
        for (int j = 0; j < n; j++) Ax_i += A[i * n + j] * x[j];
        assert(close(Ax_i, b[i], 1e-8L));
    }
    std::cout << "[PASS] test_lu_verify_ax_eq_b\n";
}

int main() {
    std::cout << "Testy: uklady rownan liniowych (LU)\n";
    test_lu_2x2();
    test_lu_3x3();
    test_lu_diagonal();
    test_lu_singular();
    test_lu_verify_ax_eq_b();
    std::cout << "Wszystkie testy LU: OK\n\n";
    return 0;
}
