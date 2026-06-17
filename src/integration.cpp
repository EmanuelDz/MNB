#include "../include/integration.h"
#include <vector>
#include <cmath>

namespace numlib {

long double trapezoidal(long double (*func)(long double), long double a, long double b, int n) {
    long double h = (b - a) / n;
    long double sum = 0.5L * (func(a) + func(b));
    for (int i = 1; i < n; i++) {
        sum += func(a + i * h);
    }
    return sum * h;
}

long double simpson(long double (*func)(long double), long double a, long double b, int n) {
    long double h = (b - a) / (2.0L * n);
    long double sum = func(a) + func(b);
    for (int i = 1; i <= n; i++) {
        sum += 4.0L * func(a + (2 * i - 1) * h);
    }
    for (int i = 1; i < n; i++) {
        sum += 2.0L * func(a + (2 * i) * h);
    }
    return (h / 3.0L) * sum;
}

long double gauss_quadrature(long double (*func)(long double), long double a, long double b, int nodes) {
    std::vector<long double> x, w;

    if (nodes == 2) {
        x = {-1.0L / sqrtl(3.0L), 1.0L / sqrtl(3.0L)};
        w = {1.0L, 1.0L};
    } else if (nodes == 3) {
        long double x1 = sqrtl(3.0L / 5.0L);
        x = {-x1, 0.0L, x1};
        w = {5.0L / 9.0L, 8.0L / 9.0L, 5.0L / 9.0L};
    } else if (nodes == 4) {
        long double x1 = sqrtl(3.0L / 7.0L - 2.0L / 7.0L * sqrtl(6.0L / 5.0L));
        long double x2 = sqrtl(3.0L / 7.0L + 2.0L / 7.0L * sqrtl(6.0L / 5.0L));
        long double w1 = (18.0L + sqrtl(30.0L)) / 36.0L;
        long double w2 = (18.0L - sqrtl(30.0L)) / 36.0L;
        x = {-x2, -x1, x1, x2};
        w = {w2, w1, w1, w2};
    } else {
        return 0.0L; // nieobslugiwana liczba wezlow
    }

    long double mid  = (b + a) / 2.0L;
    long double half = (b - a) / 2.0L;
    long double sum  = 0.0L;
    for (int i = 0; i < nodes; i++) {
        sum += w[i] * func(mid + half * x[i]);
    }
    return half * sum;
}

long double gauss_composite(long double (*func)(long double), long double a, long double b, int nodes, int intervals) {
    long double h = (b - a) / intervals;
    long double sum = 0.0L;
    for (int j = 0; j < intervals; j++) {
        long double subA = a + j * h;
        long double subB = subA + h;
        sum += gauss_quadrature(func, subA, subB, nodes);
    }
    return sum;
}

} // namespace numlib
