/**
 * @file interpolation.h
 * @brief Interpolacja wielomianowa: Newton, Lagrange, Horner.
 */
#pragma once
#include <vector>

namespace numlib {

/**
 * @brief Oblicza wspolczynniki wielomianu Newtona (ilorazy roznicowe).
 *
 * @param xi   Wezly interpolacji (musza byc rozne)
 * @param fxi  Wartosci funkcji w wezlach
 * @return     Wspolczynniki [a0, a1, ..., a_{n-1}] wielomianu Newtona
 *
 * @example
 *   auto c = numlib::newton_coefficients({0,1,2}, {1,2,5});
 *   // c[0]=1, c[1]=1, c[2]=1  =>  p(x) = 1 + 1*(x-0) + 1*(x-0)*(x-1)
 */
std::vector<long double> newton_coefficients(const std::vector<long double>& xi, const std::vector<long double>& fxi);

/**
 * @brief Oblicza wartosc wielomianu Newtona w punkcie x (schemat Hornera).
 *
 * @param x       Punkt obliczeniowy
 * @param xi      Wezly interpolacji (identyczne jak przy wywolaniu newton_coefficients)
 * @param coeffs  Wspolczynniki z newton_coefficients()
 * @return        Wartosc wielomianu w punkcie x
 *
 * @example
 *   auto c = numlib::newton_coefficients({0,1,2}, {1,2,5});
 *   long double y = numlib::newton_evaluate(1.5L, {0,1,2}, c); // 3.25
 */
long double newton_evaluate(long double x, const std::vector<long double>& xi, const std::vector<long double>& coeffs);

/**
 * @brief Oblicza wartosc wielomianu interpolacyjnego Lagrange'a w punkcie x.
 *
 * @param x    Punkt obliczeniowy
 * @param xi   Wezly interpolacji
 * @param fxi  Wartosci funkcji w wezlach
 * @return     Wartosc wielomianu Lagrange'a
 *
 * @example
 *   long double y = numlib::lagrange_evaluate(1.5L, {0,1,2}, {1,2,5}); // 3.25
 */
long double lagrange_evaluate(long double x, const std::vector<long double>& xi, const std::vector<long double>& fxi);

/**
 * @brief Wyznacza wartosc wielomianu metoda Hornera (optymalnie obliczeniowo).
 *
 * Wielomian: p(x) = a[0] + a[1]*x + a[2]*x^2 + ... + a[n-1]*x^{n-1}
 *
 * @param a  Wspolczynniki [a0, a1, ..., a_{n-1}] (od najnizszej potegi)
 * @param x  Punkt obliczeniowy
 * @return   Wartosc wielomianu
 *
 * @example
 *   // 1 + 2x + 3x^2 w x=2: 1+4+12=17
 *   long double y = numlib::horner({1.0L, 2.0L, 3.0L}, 2.0L); // 17
 */
long double horner(const std::vector<long double>& a, long double x);

/**
 * @brief Wyznacza wartosc wielomianu metoda naturalna (suma a[i]*x^i).
 *
 * @param a  Wspolczynniki [a0, a1, ..., a_{n-1}]
 * @param x  Punkt obliczeniowy
 * @return   Wartosc wielomianu
 *
 * @example
 *   long double y = numlib::polynomial_natural({1.0L, 2.0L, 3.0L}, 2.0L); // 17
 */
long double polynomial_natural(const std::vector<long double>& a, long double x);

} // namespace numlib
