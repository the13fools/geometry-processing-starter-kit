// MyConfig.cpp

#include "MyConfig.h"
#include <iostream>

// Define default values for MyConfig members
MyConfig::MyConfig()
    : w_bound(1.0),
      w_smooth(0.1),
      w_smooth_vector(0.1),
      w_curl(0.1),
      w_s_perp(0.1),
      w_attenuate(0.1),
      convergence_eps(1e-6),
      identity_weight(0.1),
      prev_energy(0.0),
      useProjHessian(true) {
}

// Implement any other MyConfig-related functionality here if needed
