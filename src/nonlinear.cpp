#include "../include/nonlinear.h"
#include <cmath>
#include <iostream>
#include <limits>

namespace numlib {

long double numerical_derivative(long double (*f)(long double), long double x, long double eps) {
    return (f(x + eps) - f(x - eps)) / (2.0L * eps);
}

long double newton_analytic(long double (*f)(long double), long double (*df)(long double), long double x0, long double tol, int maxIter, bool verbose) {
    long double x = x0;
    for (int i = 0; i < maxIter; i++) {
        long double fx  = f(x);
        long double dfx = df(x);
        if (fabsl(dfx) < 1e-15L) {
            if (verbose) std::cout << "[Newton] Pochodna zerowa w iter " << i << "\n";
            return std::numeric_limits<long double>::quiet_NaN();
        }
        x -= fx / dfx;
        if (verbose) std::cout << "[Newton] iter " << i + 1 << ":  x = " << x << "  f(x) = " << f(x) << "\n";
        if (fabsl(f(x)) < tol) break;
    }
    return x;
}

long double newton_numeric(long double (*f)(long double), long double x0, long double tol, int maxIter, bool verbose) {
    long double x = x0;
    for (int i = 0; i < maxIter; i++) {
        long double fx  = f(x);
        long double dfx = numerical_derivative(f, x);
        if (fabsl(dfx) < 1e-15L) {
            if (verbose) std::cout << "[Newton-num] Pochodna zerowa w iter " << i << "\n";
            return std::numeric_limits<long double>::quiet_NaN();
        }
        x -= fx / dfx;
        if (verbose) std::cout << "[Newton-num] iter " << i + 1 << ":  x = " << x << "  f(x) = " << f(x) << "\n";
        if (fabsl(f(x)) < tol) break;
    }
    return x;
}

long double secant(long double (*f)(long double), long double x0, long double x1, long double tol, int maxIter, bool verbose) {
    long double xa = x0, xb = x1;
    for (int i = 0; i < maxIter; i++) {
        long double fa = f(xa), fb = f(xb);
        if (fabsl(fb - fa) < 1e-15L) {
            if (verbose) std::cout << "[Sieczne] Dzielenie przez zero w iter " << i << "\n";
            return std::numeric_limits<long double>::quiet_NaN();
        }
        long double xc = xb - fb * (xb - xa) / (fb - fa);
        xa = xb;
        xb = xc;
        if (verbose) std::cout << "[Sieczne] iter " << i + 1 << ":  x = " << xb << "  f(x) = " << f(xb) << "\n";
        if (fabsl(f(xb)) < tol) break;
    }
    return xb;
}

} // namespace numlib
