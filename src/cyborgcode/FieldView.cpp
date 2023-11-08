#include "FieldView.h"

// Utility function to convert Field_View enums to strings
std::string fieldViewToString(Field_View view) {
    switch (view) {
        case vec_norms: return "Vector Norms";
        case delta_norms: return "Delta Norms";
        case vec_dirch: return "Vec Dirichlet";
        case moment_dirch: return "Moment Dirichlet";
        case primal_curl_residual: return "Primal Curl Residual";
        case sym_curl_residual: return "Symmetric Curl Residual";
        case gui_free: return "Free";
        default: return "Unknown";
    }
}
