/**
 * @file linear_systems.h
 * @brief Rozwiazywanie ukladow rownan liniowych metoda rozkladu LU.
 */
#pragma once
#include <vector>

namespace numlib {

/**
 * @brief Rozklad LU macierzy kwadratowej (bez przestawiania wierszy).
 *
 * Rozkłada macierz A = L*U, gdzie L jest macierza dolna trojkatna z jedynkami
 * na przekatnej, a U jest macierza gorna trojkatna.
 *
 * @param A   Macierz wejsciowa zapisana wierszami w wektorze 1D (rozmiar n*n)
 * @param n   Rozmiar macierzy
 * @param L   [out] Macierz dolna trojkatna (n x n)
 * @param U   [out] Macierz gorna trojkatna (n x n)
 * @return    true jezeli rozklad sie powiodl, false gdy na przekatnej U pojawi sie zero
 *
 * @example
 *   vector<long double> A = {2,1,-1, -3,-1,2, -2,1,2};
 *   vector<vector<long double>> L, U;
 *   bool ok = numlib::lu_decompose(A, 3, L, U);
 */
bool lu_decompose(const std::vector<long double>& A, int n, std::vector<std::vector<long double>>& L, std::vector<std::vector<long double>>& U);

/**
 * @brief Rozwiazuje uklad rownan liniowych Ax = b metoda rozkladu LU.
 *
 * Wewnetrznie wykonuje rozklad LU, nastepnie podstawianie w przod (Lz=b)
 * i podstawianie wstecz (Ux=z).
 *
 * @param A   Macierz wspolczynnikow (wierszami w wektorze 1D), rozmiar n*n
 * @param b   Wektor prawej strony, rozmiar n
 * @param n   Rozmiar ukladu
 * @return    Wektor rozwiazania x; pusty jezeli rozklad jest niemozliwy
 *
 * @example
 *   // 2x + y = 5,  4x + 3y = 11  ->  x=2, y=1
 *   vector<long double> A = {2,1, 4,3};
 *   vector<long double> b = {5, 11};
 *   auto x = numlib::lu_solve(A, b, 2);
 */
std::vector<long double> lu_solve(const std::vector<long double>& A, const std::vector<long double>& b, int n);

} // namespace numlib
