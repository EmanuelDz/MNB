#include "../include/interpolation.h"
#include <cmath>

namespace numlib {

std::vector<long double> newton_coefficients(const std::vector<long double>& xi, const std::vector<long double>& fxi) {
    int n = (int)xi.size();
    std::vector<long double> a = fxi; // kopia - beda modyfikowane in-place
    for (int i = 1; i < n; i++) {
        for (int j = n - 1; j >= i; j--) {
            a[j] = (a[j] - a[j - 1]) / (xi[j] - xi[j - i]);
        }
    }
    return a;
}

long double newton_evaluate(long double x, const std::vector<long double>& xi, const std::vector<long double>& coeffs) {
    int n = (int)coeffs.size();
    if (n == 0) return 0.0L;
    // Schemat Hornera dla postaci Newtona
    long double result = coeffs[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        result = result * (x - xi[i]) + coeffs[i];
    }
    return result;
}

long double lagrange_evaluate(long double x, const std::vector<long double>& xi, const std::vector<long double>& fxi) {
    int n = (int)xi.size();
    long double result = 0.0L;
    for (int i = 0; i < n; i++) {
        long double li = 1.0L;
        for (int j = 0; j < n; j++) {
            if (i != j && xi[i] != xi[j]) {
                li *= (x - xi[j]) / (xi[i] - xi[j]);
            }
        }
        result += fxi[i] * li;
    }
    return result;
}

long double horner(const std::vector<long double>& a, long double x) {
    if (a.empty()) return 0.0L;
    int n = (int)a.size();
    long double result = a[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        result = result * x + a[i];
    }
    return result;
}

long double polynomial_natural(const std::vector<long double>& a, long double x) {
    long double result = 0.0L;
    for (int i = 0; i < (int)a.size(); i++) {
        result += a[i] * powl(x, (long double)i);
    }
    return result;
}

} // namespace numlib
