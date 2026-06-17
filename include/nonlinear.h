/**
 * @file nonlinear.h
 * @brief Rozwiazywanie rownan nieliniowych f(x) = 0:
 *        metoda Newtona (analityczna i numeryczna) oraz metoda siecznych.
 */
#pragma once

namespace numlib {

/**
 * @brief Numeryczne wyznaczanie pochodnej metoda roznic centralnych.
 *
 * f'(x) ≈ (f(x + eps) - f(x - eps)) / (2 * eps)
 *
 * @param f    Funkcja
 * @param x    Punkt, w ktorym liczymy pochodna
 * @param eps  Krok (domyslnie 1e-7)
 * @return     Przyblizenie f'(x)
 *
 * @example
 *   long double d = numlib::numerical_derivative([](long double x){ return x*x; }, 3.0L);
 *   // d ~ 6.0
 */
long double numerical_derivative(long double (*f)(long double),long double x,long double eps = 1e-7L);

/**
 * @brief Metoda Newtona z analityczna pochodna - szuka pierwiastka f(x) = 0.
 *
 * Schemat: x_{n+1} = x_n - f(x_n) / f'(x_n)
 * Zbieznosc kwadratowa w poblizu prostego pierwiastka.
 *
 * @param f        Funkcja, ktorej pierwiastek szukamy
 * @param df       Analityczna pochodna funkcji f
 * @param x0       Przyblizenie poczatkowe
 * @param tol      Tolerancja zbieznosci |f(x)| < tol (domyslnie 1e-12)
 * @param maxIter  Maksymalna liczba iteracji (domyslnie 100)
 * @param verbose  Czy wypisywac kolejne iteracje na stdout (domyslnie false)
 * @return         Przyblizenie pierwiastka; NaN gdy brak zbieznosci
 *
 * @example
 *   auto f  = [](long double x){ return x*x - 2; };
 *   auto df = [](long double x){ return 2*x; };
 *   long double root = numlib::newton_analytic(f, df, 1.5L); // ~ sqrt(2)
 */
long double newton_analytic(long double (*f)(long double),long double (*df)(long double), long double x0, long double tol = 1e-12L, int maxIter = 100, bool verbose = false);
/**
 * @brief Metoda Newtona z numeryczna pochodna (rozniczkowanie centralne).
 *
 * Nie wymaga podania analitycznej pochodnej.
 *
 * @param f        Funkcja, ktorej pierwiastek szukamy
 * @param x0       Przyblizenie poczatkowe
 * @param tol      Tolerancja zbieznosci
 * @param maxIter  Maksymalna liczba iteracji
 * @param verbose  Wypisywanie iteracji
 * @return         Przyblizenie pierwiastka; NaN gdy brak zbieznosci
 *
 * @example
 *   long double root = numlib::newton_numeric([](long double x){ return x*x - 2; }, 1.5L);
 */
long double newton_numeric(long double (*f)(long double),long double x0, long double tol = 1e-12L, int maxIter = 100, bool verbose = false);

/**
 * @brief Metoda siecznych - szuka pierwiastka f(x) = 0.
 *
 * Schemat: x_{n+1} = x_n - f(x_n) * (x_n - x_{n-1}) / (f(x_n) - f(x_{n-1}))
 * Nie wymaga pochodnej.
 *
 * @param f		   Funkcja
 * @param x0       Pierwsze przyblizenie poczatkowe
 * @param x1       Drugie przyblizenie poczatkowe (rozne od x0)
 * @param tol      Tolerancja zbieznosci
 * @param maxIter  Maksymalna liczba iteracji
 * @param verbose  Wypisywanie iteracji
 * @return         Przyblizenie pierwiastka; NaN przy dzieleniu przez zero
 *
 * @example
 *   long double root = numlib::secant([](long double x){ return x*x - 2; }, 1.0L, 2.0L);
 */
long double secant(long double (*f)(long double),long double x0, long double x1, long double tol = 1e-12L, int maxIter = 100, bool verbose = false);

} // namespace numlib
