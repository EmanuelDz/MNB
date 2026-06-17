/**
 * @file integration.h
 * @brief Calkowanie numeryczne: metoda trapezow, Simpsona, kwadratura Gaussa-Legendre'a.
 */
#pragma once

namespace numlib {

/**
 * @brief Calkowanie numeryczne metoda trapezow (zlozony wzor trapezow).
 *
 * Podzial przedzialu [a,b] na n rownoodleglych podprzedzialow.
 *
 * @param func  Funkcja podcalkowa f(x)
 * @param a     Lewy kraniec przedzialu
 * @param b     Prawy kraniec przedzialu
 * @param n     Liczba podprzedzialow (im wieksza, tym dokladniejsza)
 * @return      Przyblizenie calki oznaczonej
 *
 * @example
 *   long double I = numlib::trapezoidal([](long double x){ return x*x; }, 0.0L, 1.0L, 1000);
 *   // I ~ 0.333333...
 */
long double trapezoidal(long double (*func)(long double), long double a, long double b, int n);

/**
 * @brief Calkowanie numeryczne wzorem Simpsona (1/3).
 *
 * Podzial [a,b] na 2n podprzedzialow.
 *
 * @param func  Funkcja podcalkowa
 * @param a     Lewy kraniec
 * @param b     Prawy kraniec
 * @param n     Liczba par podprzedzialow (lacznie 2n punktow wewnetrznych)
 * @return      Przyblizenie calki
 *
 * @example
 *   long double I = numlib::simpson([](long double x){ return sinl(x); }, 0.0L, M_PIl, 100);
 *   // I ~ 2.0
 */
long double simpson(long double (*func)(long double), long double a, long double b, int n);

/**
 * @brief Kwadratura Gaussa na jednym przedziale [a,b].
 *
 * Obslugiwane rzedy: nodes = 2, 3 lub 4.
 * Dla wielomianow stopnia < 2*nodes wynik jest dokladny.
 *
 * @param func   Funkcja podcalkowa
 * @param a      Lewy kraniec
 * @param b      Prawy kraniec
 * @param nodes  Liczba wezlow: 2, 3 lub 4
 * @return       Przyblizenie calki (0.0 dla nieobslugiwanych wartosci nodes)
 *
 * @example
 *   long double I = numlib::gauss_quadrature(myFunc, 0.0L, 1.0L, 4);
 */
long double gauss_quadrature(long double (*func)(long double), long double a, long double b, int nodes);

/**
 * @brief Zlozona kwadratura Gaussa - podzial na m podprzedzialow.
 *
 * Podzial [a,b] na m rownych czesci; na kazdej stosowana kwadratura Gaussa
 * z podana liczba wezlow.
 *
 * @param func       Funkcja podcalkowa
 * @param a          Lewy kraniec
 * @param b          Prawy kraniec
 * @param nodes      Liczba wezlow na podprzedzial (2, 3 lub 4)
 * @param intervals  Liczba podprzedzialow
 * @return           Przyblizenie calki
 *
 * @example
 *   long double I = numlib::gauss_composite(f, -1.0L, 15.0L, 4, 50);
 */
long double gauss_composite(long double (*func)(long double), long double a, long double b, int nodes, int intervals);

} // namespace numlib
