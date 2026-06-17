#include "../include/ode.h"

namespace numlib {

ODEResult euler(ODEFunc f, long double y0, long double t0, long double t_end, int N) {
    long double h = (t_end - t0) / N;
    ODEResult result;
    result.reserve(N + 1);
    long double t = t0, y = y0;
    for (int i = 0; i <= N; i++) {
        result.push_back({t, y});
        long double yn = y + h * f(t, y);
        t = t0 + (i + 1) * h;
        y = yn;
    }
    return result;
}

ODEResult heun(ODEFunc f, long double y0, long double t0, long double t_end, int N) {
    long double h = (t_end - t0) / N;
    ODEResult result;
    result.reserve(N + 1);
    long double t = t0, y = y0;
    for (int i = 0; i <= N; i++) {
        result.push_back({t, y});
        long double k1 = f(t, y);
        long double k2 = f(t + h, y + h * k1);
        long double yn = y + h * 0.5L * (k1 + k2);
        t = t0 + (i + 1) * h;
        y = yn;
    }
    return result;
}

ODEResult midpoint(ODEFunc f, long double y0, long double t0, long double t_end, int N) {
    long double h = (t_end - t0) / N;
    ODEResult result;
    result.reserve(N + 1);
    long double t = t0, y = y0;
    for (int i = 0; i <= N; i++) {
        result.push_back({t, y});
        long double k1 = f(t, y);
        long double k2 = f(t + h * 0.5L, y + h * 0.5L * k1);
        long double yn = y + h * k2;
        t = t0 + (i + 1) * h;
        y = yn;
    }
    return result;
}

ODEResult rk4(ODEFunc f, long double y0, long double t0, long double t_end, int N) {
    long double h = (t_end - t0) / N;
    ODEResult result;
    result.reserve(N + 1);
    long double t = t0, y = y0;
    for (int i = 0; i <= N; i++) {
        result.push_back({t, y});
        long double k1 = f(t,            y);
        long double k2 = f(t + h / 2.0L, y + h / 2.0L * k1);
        long double k3 = f(t + h / 2.0L, y + h / 2.0L * k2);
        long double k4 = f(t + h,        y + h * k3);
        long double yn = y + (h / 6.0L) * (k1 + 2.0L * k2 + 2.0L * k3 + k4);
        t = t0 + (i + 1) * h;
        y = yn;
    }
    return result;
}

} // namespace numlib
