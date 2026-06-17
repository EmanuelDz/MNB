/**
 * @file test_ode.cpp
 * @brief Testy jednostkowe dla metod rozwiazywania rownan rozniczkowych.
 *
 * Zadanie testowe: y' = -y, y(0) = 1  =>  rozwiazanie dokladne: y(t) = e^{-t}
 */
#include "../include/numlib.h"
#include <iostream>
#include <cmath>
#include <cassert>

static bool close(long double a, long double b, long double tol) {
    return fabsl(a - b) < tol;
}

// Prawa strona: y' = -y
static long double f_neg_y(long double /*t*/, long double y) { return -y; }

// Prawa strona: y' = t + y, y(0) = 0  =>  y(t) = e^t - t - 1
static long double f_t_plus_y(long double t, long double y) { return t + y; }

// Test 1: Euler — y' = -y, y(0) = 1, sprawdz y(1) ~ e^{-1}
void test_euler_decay() {
    auto sol = numlib::euler(f_neg_y, 1.0L, 0.0L, 1.0L, 10000);
    long double y_num   = sol.back().second;
    long double y_exact = expl(-1.0L);
    assert(close(y_num, y_exact, 1e-4L));
    std::cout << "[PASS] test_euler_decay (blad=" << fabsl(y_num - y_exact) << ")\n";
}

// Test 2: Heun — y' = -y, y(0) = 1, sprawdz y(1) ~ e^{-1} z wyzsza dokladnoscia
void test_heun_decay() {
    auto sol = numlib::heun(f_neg_y, 1.0L, 0.0L, 1.0L, 1000);
    long double y_num   = sol.back().second;
    long double y_exact = expl(-1.0L);
    assert(close(y_num, y_exact, 1e-6L));
    std::cout << "[PASS] test_heun_decay (blad=" << fabsl(y_num - y_exact) << ")\n";
}

// Test 3: Punkt srodkowy — y' = -y, y(0) = 1
void test_midpoint_decay() {
    auto sol = numlib::midpoint(f_neg_y, 1.0L, 0.0L, 1.0L, 1000);
    long double y_num   = sol.back().second;
    long double y_exact = expl(-1.0L);
    assert(close(y_num, y_exact, 1e-6L));
    std::cout << "[PASS] test_midpoint_decay (blad=" << fabsl(y_num - y_exact) << ")\n";
}

// Test 4: RK4 — y' = -y, y(0) = 1, bardzo wysoka dokladnosc
void test_rk4_decay() {
    auto sol = numlib::rk4(f_neg_y, 1.0L, 0.0L, 1.0L, 100);
    long double y_num   = sol.back().second;
    long double y_exact = expl(-1.0L);
    assert(close(y_num, y_exact, 1e-10L));
    std::cout << "[PASS] test_rk4_decay (blad=" << fabsl(y_num - y_exact) << ")\n";
}

// Test 5: RK4 bardziej dokladny niz Euler dla tej samej liczby krokow
void test_rk4_better_than_euler() {
    int N = 50;
    long double y_exact = expl(-1.0L);
    auto solE = numlib::euler(f_neg_y, 1.0L, 0.0L, 1.0L, N);
    auto solR = numlib::rk4(f_neg_y, 1.0L, 0.0L, 1.0L, N);
    long double errE = fabsl(solE.back().second - y_exact);
    long double errR = fabsl(solR.back().second - y_exact);
    assert(errR < errE);
    std::cout << "[PASS] test_rk4_better_than_euler" << " (Euler=" << errE << ", RK4=" << errR << ")\n";
}

// Test 6: Liczba punktow wyjsciowych == N+1
void test_output_size() {
    int N = 200;
    assert((int)numlib::euler(f_neg_y, 1.0L, 0.0L, 1.0L, N).size() == N + 1);
    assert((int)numlib::heun(f_neg_y, 1.0L, 0.0L, 1.0L, N).size() == N + 1);
    assert((int)numlib::midpoint(f_neg_y, 1.0L, 0.0L, 1.0L, N).size() == N + 1);
    assert((int)numlib::rk4(f_neg_y, 1.0L, 0.0L, 1.0L, N).size() == N + 1);
    std::cout << "[PASS] test_output_size\n";
}

// Test 7: RK4 — y' = t + y, y(0) = 0,  y(2) = e^2 - 3
void test_rk4_nonautonomous() {
    auto sol = numlib::rk4(f_t_plus_y, 0.0L, 0.0L, 2.0L, 1000);
    long double y_num   = sol.back().second;
    long double y_exact = expl(2.0L) - 3.0L;
    assert(close(y_num, y_exact, 1e-8L));
    std::cout << "[PASS] test_rk4_nonautonomous (blad=" << fabsl(y_num - y_exact) << ")\n";
}

int main() {
    std::cout << "Testy: rownania rozniczkowe (ODE)\n";
    test_euler_decay();
    test_heun_decay();
    test_midpoint_decay();
    test_rk4_decay();
    test_rk4_better_than_euler();
    test_output_size();
    test_rk4_nonautonomous();
    std::cout << "Wszystkie testy ODE: OK\n\n";
    return 0;
}
