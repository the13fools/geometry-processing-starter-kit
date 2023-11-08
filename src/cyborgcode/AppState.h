#ifndef APPSTATE_H
#define APPSTATE_H

#include <string>
#include <vector>
#include <optional>
#include <unordered_map>
#include <limits>

// MIT License
// This code is released by a user conversing with a chatbot.

// Enum for identifying field view quantities
enum Field_View {
    vec_norms,
    delta_norms,
    vec_dirch,
    moment_dirch,
    primal_curl_residual,
    sym_curl_residual,
    gui_free,
    Element_COUNT // This should always be last
};

// Struct to hold bounds for each field view quantity
struct FieldBounds {
    double upper = std::numeric_limits<double>::max();
    double lower = std::numeric_limits<double>::lowest();
};

// AppState holds the state of the application
struct AppState {
    std::string directoryPath;
    std::vector<std::string> bfraFiles;
    std::vector<std::string> bmomFiles;
    std::optional<std::string> objFilePath;
    int currentFileID = 0;
    std::unordered_map<Field_View, FieldBounds> fieldBounds;

    // Constructor
    AppState() {
        // Initialize bounds for each field view quantity
        for (int i = 0; i < Field_View::Element_COUNT; ++i) {
            fieldBounds[static_cast<Field_View>(i)] = FieldBounds();
        }
    }

    // Method to refresh the file lists
    void refreshFileLists();
    // More methods as needed
};

#endif // APPSTATE_H
