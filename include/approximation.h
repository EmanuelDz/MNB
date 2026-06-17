/**
 * @file approximation.h
 * @brief Aproksymacja sredniokwadratowa funkcji wielomianem algebraicznym.
 */
#pragma once
#include <vector>

namespace numlib {

/**
 * @brief Aproksymacja sredniokwadratowa funkcji wielomianem stopnia N.
 *
 * Wyznacza wspolczynniki wielomianu P(x) = a0 + a1*x + ... + aN*x^N
 * minimalizujacego blad sredniokwadratowy ||f - P||^2 na przedziale [a,b].
 *
 * Algorytm:
 *  1. Buduje macierz Grama G[i][j] = integral(x^i * x^j, a, b)
 *     kwadratura Gaussa (4 wezly, podzielona na gauss_intervals czesci).
 *  2. Oblicza wektor prawej strony rhs[i] = integral(f(x)*x^i, a, b).
 *  3. Rozwiazuje uklad normalny G*a = rhs metoda LU.
 *
 * @param func            Aproksymowana funkcja f(x)
 * @param a               Lewy kraniec przedzialu
 * @param b               Prawy kraniec przedzialu
 * @param degree          Stopien wielomianu aproksymujacego (>= 0)
 * @param gauss_intervals Liczba podprzedzialow do calkowania Gaussa (domyslnie 100)
 * @return                Wspolczynniki [a0, a1, ..., aN]; pusty gdy rozklad LU sie nie powiedzie
 *
 * @example
 *   auto coeffs = numlib::polynomial_approximation(sinl, 0.0L, 3.14159L, 4);
 *   long double y = numlib::approx_evaluate(coeffs, 1.5L);
 */
std::vector<long double> polynomial_approximation(long double (*func)(long double),long double a, long double b,int degree,int gauss_intervals = 100);

/**
 * @brief Oblicza wartosc wielomianu aproksymujacego w punkcie x.
 *
 * @param coeffs  Wspolczynniki z polynomial_approximation()
 * @param x       Punkt obliczeniowy
 * @return        Wartosc wielomianu P(x) = sum(coeffs[i] * x^i)
 *
 * @example
 *   auto c = numlib::polynomial_approximation(sinl, 0.0L, 3.14L, 4);
 *   long double y = numlib::approx_evaluate(c, 1.0L); // ~sin(1) = 0.8415
 */
long double approx_evaluate(const std::vector<long double>& coeffs, long double x);

} // namespace numlib
