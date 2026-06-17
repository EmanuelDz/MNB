/**
 * @file example_ode_interpolation.cpp
 * @brief Przyklad 2: Rownania rozniczkowe, interpolacja i aproksymacja.
 *
 * Demonstracja:
 *   1. Rozwiazanie y' = -2y, y(0) = 1 czterema metodami; wynik dokladny: e^{-2t}
 *   2. Interpolacja Newtona punktow z sinusa; ocena bledu.
 *   3. Aproksymacja sredniokwadratowa f(x) = exp(-x)*sin(3x) wielomianem 5. stopnia.
 */
#include "../include/numlib.h"
#include <iostream>
#include <iomanip>
#include <cmath>

//ODE
static long double ode_rhs(long double /*t*/, long double y) { return -2.0L * y; }

//Interpolacja
const long double PI = 3.14159265358979323846L;

//Aproksymacja
static long double h(long double x) { return expl(-x) * sinl(3.0L * x); }

int main() {
    std::cout << std::fixed << std::setprecision(10);

    // 1. ROWNANIA ROZNICZKOWE: y' = -2y, y(0)=1
    std::cout << "ODE: y' = -2y,  y(0)=1  =>  y(t) = e^(-2t)\n";
    int N = 100;
    long double t_end = 2.0L;

    auto solE  = numlib::euler(ode_rhs, 1.0L, 0.0L, t_end, N);
    auto solH  = numlib::heun(ode_rhs, 1.0L, 0.0L, t_end, N);
    auto solM  = numlib::midpoint(ode_rhs, 1.0L, 0.0L, t_end, N);
    auto solR  = numlib::rk4(ode_rhs, 1.0L, 0.0L, t_end, N);
    long double exact_end = expl(-2.0L * t_end);

    std::cout << "Wynik w t=" << t_end << ", dokladny=" << exact_end << "\n";
    std::cout << std::setw(22) << "Euler"     << " = " << solE.back().second << "  blad=" << fabsl(solE.back().second - exact_end) << "\n";
    std::cout << std::setw(22) << "Heun"      << " = " << solH.back().second << "  blad=" << fabsl(solH.back().second - exact_end) << "\n";
    std::cout << std::setw(22) << "Punkt srod." << " = " << solM.back().second << "  blad=" << fabsl(solM.back().second - exact_end) << "\n";
    std::cout << std::setw(22) << "RK4"       << " = " << solR.back().second << "  blad=" << fabsl(solR.back().second - exact_end) << "\n";

    std::cout << "\nPierwsze 5 krokow RK4 (t, y_num, y_dokl, blad):\n";
    for (int i = 0; i <= 5; i++) {
        long double t  = solR[i].first;
        long double yn = solR[i].second;
        long double ye = expl(-2.0L * t);
        std::cout << "  t=" << t << "  y_num=" << yn << "  y_dokl=" << ye << "  blad=" << fabsl(yn - ye) << "\n";
    }

    // 2. INTERPOLACJA NEWTONA — sin(x) na 8 wezlach
    std::cout << "\n=== Interpolacja Newtona: sin(x) na [0, pi] (8 wezlow) ===\n";
    int Kn = 8;
    std::vector<long double> xi, fxi;
    for (int i = 0; i < Kn; i++) {
        long double x = PI * i / (Kn - 1);
        xi.push_back(x);
        fxi.push_back(sinl(x));
    }
    auto coeffs = numlib::newton_coefficients(xi, fxi);

    std::cout << "Wspolczynniki Newtona:\n";
    for (int i = 0; i < (int)coeffs.size(); i++) std::cout << "  a[" << i << "] = " << coeffs[i] << "\n";

    std::cout << "\nWartosci interpolowane (x nie bedace wezlami):\n";
    std::cout << std::setw(12) << "x" << std::setw(18) << "sin(x)" << std::setw(18) << "interp(x)" << std::setw(16) << "blad\n";
    long double maxErr = 0.0L;
    for (int i = 0; i < 20; i++) {
        long double x  = PI * (i + 0.5L) / 20.0L;
        long double ye = sinl(x);
        long double yi = numlib::newton_evaluate(x, xi, coeffs);
        long double e  = fabsl(ye - yi);
        if (e > maxErr) maxErr = e;
        std::cout << std::setw(12) << x << std::setw(18) << ye << std::setw(18) << yi << std::setw(16) << e << "\n";
    }
    std::cout << "Maksymalny blad interpolacji: " << maxErr << "\n";

    // 3. APROKSYMACJA: h(x) = exp(-x)*sin(3x), stopien 5
    std::cout << "\nAproksymacja: exp(-x)*sin(3x) na [0,2], wielomian st.5 \n";
    auto ac = numlib::polynomial_approximation(h, 0.0L, 2.0L, 5);
    std::cout << "Wspolczynniki wielomianu aproksymujacego:\n";
    for (int i = 0; i < (int)ac.size(); i++) std::cout << "  a[" << i << "] = " << ac[i] << "\n";

    std::cout << "\nPorownanie w wybranych punktach:\n";
    std::cout << std::setw(10) << "x" << std::setw(18) << "f(x)" << std::setw(18) << "P(x)" << std::setw(14) << "blad\n";
    long double maxErrA = 0.0L;
    for (int i = 0; i <= 20; i++) {
        long double x  = 2.0L * i / 20.0L;
        long double fe = h(x);
        long double fa = numlib::approx_evaluate(ac, x);
        long double e  = fabsl(fe - fa);
        if (e > maxErrA) maxErrA = e;
        std::cout << std::setw(10) << x << std::setw(18) << fe << std::setw(18) << fa << std::setw(14) << e << "\n";
    }
    std::cout << "Maksymalny blad aproksymacji: " << maxErrA << "\n";

    // 4. UKLAD ROWNAN LINIOWYCH
    std::cout << "\n Uklad rownan (LU): 2x+y=5,  4x+3y=11 \n";
    std::vector<long double> A = {2, 1, 4, 3};
    std::vector<long double> b = {5, 11};
    auto x = numlib::lu_solve(A, b, 2);
    std::cout << "x[0] = " << x[0] << "  (oczekiwane 2)\n";
    std::cout << "x[1] = " << x[1] << "  (oczekiwane 1)\n";

    return 0;
}
