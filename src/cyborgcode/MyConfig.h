// MyConfig.h

#ifndef MY_CONFIG_H
#define MY_CONFIG_H

#include <string>

struct MyConfig {
    double w_bound;
    double w_smooth;
    double w_smooth_vector;
    double w_curl;
    double w_s_perp;

    double w_attenuate;

    double convergence_eps;
    double identity_weight;
    double prev_energy;
    bool useProjHessian;
};

#endif // MY_CONFIG_H
