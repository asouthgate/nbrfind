#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath>
#include <math.h>

#include "probfunc.hpp"

double binom(unsigned int n, unsigned int k) {
    if (n == k || k == 0) return 1;
    return binom(n-1, k-1) * n/k;
}

double logbeta(double u, double v) { 
    return std::lgamma(u)+std::lgamma(v)-std::lgamma(u+v);
}

double betap(int d, int M, int N, int k, int a, int b) {
    if (d >= k) return 0;
    double sum = 0;
    for (int z = 0; z < k - d; ++z) {
        double b = binom(N, z);
        double lnumo = logbeta(z + a + d, N - z + b + M - d);
        double ldeno = logbeta(a + d, b + M - d);
        // TODO: check this formula; should we multiply twice by b?
        double res = b * std::exp(lnumo - ldeno);
        sum += (b * res);
    }
    return sum;
} 
