#include "../include/approximation.h"
#include "../include/linear_systems.h"
#include <cmath>
#include <vector>

namespace numlib {

std::vector<long double> polynomial_approximation(long double (*func)(long double), long double a, long double b, int degree, int gauss_intervals) {
    int m = degree + 1; // liczba wspolczynnikow
    // Wezly i wagi kwadrury Gaussa (4-wezlowej)
    long double gx1 = sqrtl(3.0L / 7.0L - 2.0L / 7.0L * sqrtl(6.0L / 5.0L));
    long double gx2 = sqrtl(3.0L / 7.0L + 2.0L / 7.0L * sqrtl(6.0L / 5.0L));
    long double gw1 = (18.0L + sqrtl(30.0L)) / 36.0L;
    long double gw2 = (18.0L - sqrtl(30.0L)) / 36.0L;
    std::vector<long double> gx = {-gx2, -gx1, gx1, gx2};
    std::vector<long double> gw = {gw2, gw1, gw1, gw2};

    long double hSub = (b - a) / gauss_intervals;

    // Macierz Grama (zapisana wierszami w wektorze 1D) i wektor prawej strony
    std::vector<long double> G(m * m, 0.0L);
    std::vector<long double> rhs(m, 0.0L);

    for (int i = 0; i < m; i++) {
        // Elementy macierzy G[i][j] = integral( x^i * x^j, a, b )
        for (int j = 0; j < m; j++) {
            long double val = 0.0L;
            for (int s = 0; s < gauss_intervals; s++) {
                long double sA   = a + s * hSub;
                long double sB   = sA + hSub;
                long double mid  = (sA + sB) / 2.0L;
                long double half = (sB - sA) / 2.0L;
                long double sub  = 0.0L;
                for (int g = 0; g < 4; g++) sub += gw[g] * powl(mid + half * gx[g], (long double)(i + j));
                val += sub * half;
            }
            G[i * m + j] = val;
        }

        // Elementy prawej strony rhs[i] = integral( f(x) * x^i, a, b )
        long double val = 0.0L;
        for (int s = 0; s < gauss_intervals; s++) {
            long double sA   = a + s * hSub;
            long double sB   = sA + hSub;
            long double mid  = (sA + sB) / 2.0L;
            long double half = (sB - sA) / 2.0L;
            long double sub  = 0.0L;
            for (int g = 0; g < 4; g++) sub += gw[g] * func(mid + half * gx[g]) * powl(mid + half * gx[g], (long double)i);
            val += sub * half;
        }
        rhs[i] = val;
    }

    // Rozwiaz uklad normalny G * coeffs = rhs metoda LU
    return lu_solve(G, rhs, m);
}

long double approx_evaluate(const std::vector<long double>& coeffs, long double x) {
    long double result = 0.0L;
    for (int i = 0; i < (int)coeffs.size(); i++) {
        result += coeffs[i] * powl(x, (long double)i);
    }
    return result;
}

} // namespace numlib
