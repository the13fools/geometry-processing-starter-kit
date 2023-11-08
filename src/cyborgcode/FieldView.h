#ifndef FIELDVIEW_H
#define FIELDVIEW_H

#include <string>

/**
 * Enumeration for the different field views.
 */
enum Field_View {
    vec_norms,
    delta_norms,
    vec_dirch,
    moment_dirch,
    primal_curl_residual,
    sym_curl_residual,
    gui_free,
    Element_COUNT
};

/**
 * Converts a field view enumeration to a string.
 * 
 * @param view The field view to convert.
 * @return The name of the field view as a string.
 */
std::string fieldViewToString(Field_View view);

#endif // FIELDVIEW_H
