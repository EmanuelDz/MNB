/**
 * @file example_integration_nonlinear.cpp
 * @brief Przyklad 1: Calkowanie numeryczne i szukanie pierwiastkow.
 *
 * Demonstracja:
 *   1. Calkowanie funkcji (1-x)*exp(sin(x)) od -1 do 15 czterema metodami.
 *   2. Porownanie bledow wzgledem wartosci dokladnej.
 *   3. Szukanie pierwiastka rownania cos(x) = x metodami Newtona i siecznych.
 */
#include "../include/numlib.h"
#include <iostream>
#include <iomanip>
#include <cmath>

static long double f(long double x) {
    return (1.0L - x) * expl(sinl(x));
}
static long double g(long double x)  { return cosl(x) - x; }
static long double dg(long double x) { return -sinl(x) - 1.0L; }

int main() {
    std::cout << std::fixed << std::setprecision(10);

    // 1. CALKOWANIE NUMERYCZNE
    std::cout << "Calkowanie (1-x)*exp(sin(x)) od -1 do 15\n";
    long double exact = -136.183666144466L;
    long double a = -1.0L, b = 15.0L;
    long double I_trap  = numlib::trapezoidal(f, a, b, 50);
    long double I_simp  = numlib::simpson(f, a, b, 50);
    long double I_g4    = numlib::gauss_composite(f, a, b, 4, 50);
    std::cout << "Metoda trapezow (n=50)     : " << I_trap << "  blad=" << fabsl(I_trap  - exact) << "\n";
    std::cout << "Metoda Simpsona (n=50)     : " << I_simp << "  blad=" << fabsl(I_simp  - exact) << "\n";
    std::cout << "Gauss 4-wezl. (m=50 przed.): " << I_g4 << "  blad=" << fabsl(I_g4    - exact) << "\n";

    std::cout << "\nKwadratura Gaussa - wplyw liczby wezlow i przedzialow:\n";
    std::cout << std::setw(8) << "wezly" << std::setw(12) << "przedzialow" << std::setw(22) << "wynik" << std::setw(16) << "blad\n";
    for (int nodes : {2, 3, 4}) {
        for (int m : {10, 20, 50}) {
            long double I = numlib::gauss_composite(f, a, b, nodes, m);
            std::cout << std::setw(8)  << nodes << std::setw(12) << m << std::setw(22) << I << std::setw(16) << fabsl(I - exact) << "\n";
        }
    }

    // 2. ROWNANIA NIELINIOWE: cos(x) = x
    std::cout << "\n=== Szukanie pierwiastka: cos(x) - x = 0 ===\n";
    long double rootNA = numlib::newton_analytic(g, dg, 0.5L, 1e-12L, 100, true);
    std::cout << "Newton (analityczny): x = " << rootNA << "  f(x) = " << g(rootNA) << "\n\n";

    long double rootNN = numlib::newton_numeric(g, 0.5L, 1e-12L, 100, true);
    std::cout << "Newton (numeryczny):  x = " << rootNN << "  f(x) = " << g(rootNN) << "\n\n";

    long double rootSec = numlib::secant(g, 0.0L, 1.0L, 1e-12L, 100, true);
    std::cout << "Metoda siecznych:     x = " << rootSec << "  f(x) = " << g(rootSec) << "\n";

    return 0;
}
