/**
 * @file numlib.h
 * @brief Glowny plik naglowkowy biblioteki numerycznej numlib.
 *
 * Biblioteka numlib implementuje nastepujace kategorie metod numerycznych:
 *  - Rozwiazywanie ukladow rownan liniowych (rozklad LU)
 *  - Interpolacja (Newton, Lagrange, Horner)
 *  - Aproksymacja sredniokwadratowa wielomianem
 *  - Calkowanie numeryczne (trapezy, Simpson, kwadratura Gaussa-Legendre'a)
 *  - Rozwiazywanie rownan rozniczkowych (Euler, Heun, punkt srodkowy, RK4)
 *  - Rozwiazywanie rownan nieliniowych (Newton analityczny/numeryczny, sieczne)
 *
 * Uzycie:
 * @code
 *   #include "numlib.h"
 *   long double I = numlib::trapezoidal(myFunc, 0.0L, 1.0L, 1000);
 * @endcode
 */
#pragma once

#include "linear_systems.h"
#include "interpolation.h"
#include "approximation.h"
#include "integration.h"
#include "ode.h"
#include "nonlinear.h"
