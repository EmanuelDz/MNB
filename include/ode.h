/**
 * @file ode.h
 * @brief Rozwiazywanie rownan rozniczkowych zwyczajnych (RRZ):
 *        Euler, Heun, punkt srodkowy, Runge-Kutta 4. rzedu.
 */
#pragma once
#include <vector>
#include <utility>

namespace numlib {

/// Typ wskaznika na prawą stronę równania: y'(t) = f(t, y)
using ODEFunc = long double (*)(long double t, long double y);

/// Wynik: wektor par (t_i, y_i)
using ODEResult = std::vector<std::pair<long double, long double>>;

/**
 * @brief Metoda Eulera (1. rzad) — rozwiazuje y'(t) = f(t,y), y(t0) = y0.
 *
 * Schemat: y_{i+1} = y_i + h * f(t_i, y_i)
 *
 * @param f      Prawa strona rownania f(t, y)
 * @param y0     Warunek poczatkowy y(t0) = y0
 * @param t0     Poczatek przedzialu
 * @param t_end  Koniec przedzialu
 * @param N      Liczba krokow czasowych
 * @return       Wektor par (t_i, y_i) dla i = 0, 1, ..., N
 *
 * @example
 *   // y' = -y, y(0) = 1  ->  y(t) = e^(-t)
 *   auto f = [](long double t, long double y){ return -y; };
 *   auto sol = numlib::euler(f, 1.0L, 0.0L, 5.0L, 10000);
 */
ODEResult euler(ODEFunc f, long double y0, long double t0, long double t_end, int N);

/**
 * @brief Metoda Heuna (2. rzad).
 *
 * Schemat predyktor-korektor:
 *   k1 = f(t_i, y_i)
 *   k2 = f(t_i + h, y_i + h*k1)
 *   y_{i+1} = y_i + h/2 * (k1 + k2)
 *
 * @param f      Prawa strona rownania
 * @param y0     Warunek poczatkowy
 * @param t0     Poczatek przedzialu
 * @param t_end  Koniec przedzialu
 * @param N      Liczba krokow
 * @return       Wektor par (t_i, y_i)
 *
 * @example
 *   auto sol = numlib::heun(f, 1.0L, 0.0L, 5.0L, 1000);
 */
ODEResult heun(ODEFunc f, long double y0, long double t0, long double t_end, int N);

/**
 * @brief Metoda punktu srodkowego (2. rzad).
 *
 * Schemat:
 *   k1 = f(t_i, y_i)
 *   k2 = f(t_i + h/2, y_i + h/2 * k1)
 *   y_{i+1} = y_i + h * k2
 *
 * @param f      Prawa strona rownania
 * @param y0     Warunek poczatkowy
 * @param t0     Poczatek przedzialu
 * @param t_end  Koniec przedzialu
 * @param N      Liczba krokow
 * @return       Wektor par (t_i, y_i)
 *
 * @example
 *   auto sol = numlib::midpoint(f, 1.0L, 0.0L, 5.0L, 1000);
 */
ODEResult midpoint(ODEFunc f, long double y0, long double t0, long double t_end, int N);

/**
 * @brief Metoda Rungego-Kutty 4. rzedu (RK4).
 *
 * Schemat:
 *   k1 = f(t_i, y_i)
 *   k2 = f(t_i + h/2, y_i + h/2*k1)
 *   k3 = f(t_i + h/2, y_i + h/2*k2)
 *   k4 = f(t_i + h,   y_i + h*k3)
 *   y_{i+1} = y_i + h/6 * (k1 + 2*k2 + 2*k3 + k4)
 *
 *
 * @param f      Prawa strona rownania
 * @param y0     Warunek poczatkowy
 * @param t0     Poczatek przedzialu
 * @param t_end  Koniec przedzialu
 * @param N      Liczba krokow
 * @return       Wektor par (t_i, y_i)
 *
 * @example
 *   auto sol = numlib::rk4(f, 1.0L, 0.0L, 5.0L, 100);
 */
ODEResult rk4(ODEFunc f, long double y0, long double t0, long double t_end, int N);

} // namespace numlib
