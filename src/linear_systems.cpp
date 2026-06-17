#include "../include/linear_systems.h"
#include <iostream>
#include <cmath>

namespace numlib {

bool lu_decompose(const std::vector<long double>& A, int n, std::vector<std::vector<long double>>& L, std::vector<std::vector<long double>>& U) {
    L.assign(n, std::vector<long double>(n, 0.0L));
    U.assign(n, std::vector<long double>(n, 0.0L));
    for (int i = 0; i < n; i++) L[i][i] = 1.0L;
    for (int k = 0; k < n; k++) {
        // Wyznacz wiersz k macierzy U
        for (int j = k; j < n; j++) {
            long double sum = 0;
            for (int s = 0; s < k; s++) sum += L[k][s] * U[s][j];
            U[k][j] = A[k * n + j] - sum;
        }
        // Wyznacz kolumne k macierzy L
        for (int i = k + 1; i < n; i++) {
            if (fabsl(U[k][k]) < 1e-15L) {
                std::cerr << "Blad LU: zero na przekatnej U[" << k << "][" << k << "]. Rozklad niemozliwy.\n";
                return false;
            }
            long double sum = 0;
            for (int s = 0; s < k; s++) sum += L[i][s] * U[s][k];
            L[i][k] = (A[i * n + k] - sum) / U[k][k];
        }
    }
    return true;
}

std::vector<long double> lu_solve(const std::vector<long double>& A, const std::vector<long double>& b, int n) {
    std::vector<std::vector<long double>> L, U;
    if (!lu_decompose(A, n, L, U)) return {};

    // Podstawianie w przod: Lz = b
    std::vector<long double> z(n, 0.0L);
    for (int i = 0; i < n; i++) {
        long double sum = 0;
        for (int j = 0; j < i; j++) sum += L[i][j] * z[j];
        z[i] = b[i] - sum;
    }

    // Podstawianie wstecz: Ux = z
    std::vector<long double> x(n, 0.0L);
    for (int i = n - 1; i >= 0; i--) {
        if (fabsl(U[i][i]) < 1e-15L) {
            std::cerr << "Blad LU: zero na przekatnej U[" << i << "][" << i << "]. Rozwiazanie niemozliwe.\n";
            return {};
        }
        long double sum = 0;
        for (int j = i + 1; j < n; j++) sum += U[i][j] * x[j];
        x[i] = (z[i] - sum) / U[i][i];
    }
    return x;
}

} // namespace numlib
