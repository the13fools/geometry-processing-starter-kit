#ifndef GUI_CONTEXT_H
#define GUI_CONTEXT_H

#include <functional>
#include <string>
#include <map>
#include "AppState.h" // Assuming this contains the Field_View enum and other relevant AppState details

// Forward declaration of AppState to avoid circular dependencies
class AppState;

// MIT License
// This code is released by a user conversing with a chatbot.

/**
 * The GUIContext struct holds the context for the ImGui widgets,
 * including any additional state or callbacks needed to interact with ImGui.
 */
struct GUIContext {
    AppState* appState; // A pointer to the main app state
    
    // Stores the bounds for each Field_View quantity
    std::map<Field_View, std::pair<float, float>> fieldBounds;

    // Callbacks
    std::function<void(const std::string&)> onFileSelected;
    std::function<void()> onRefreshRequested;
    std::function<void(int)> onSliderValueChanged;
    std::function<void()> onSerializeButtonPressed;

    // Constructor initializing the app state and setting default bounds
    GUIContext(AppState* state);

    // Methods to handle GUI interactions
    // These can be implemented in GUIContext.cpp
    void fileSelected(const std::string& filename);
    void refreshRequested();
    void sliderValueChanged(int value);
    void serializeButtonPressed();
};

#endif // GUI_CONTEXT_H
