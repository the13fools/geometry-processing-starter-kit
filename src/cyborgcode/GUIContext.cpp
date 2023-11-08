#include "GUIContext.h"
#include "ImGui/imgui.h"

// Constructor initializing the app state and setting default bounds
GUIContext::GUIContext(AppState* state) : appState(state) {
    // Initialize default bounds for each Field_View quantity
    for (int i = 0; i < Field_View::Element_COUNT; ++i) {
        fieldBounds[static_cast<Field_View>(i)] = std::make_pair(0.0f, 1.0f);
    }
}

// Method to handle file selection in GUI
void GUIContext::fileSelected(const std::string& filename) {
    if (onFileSelected) {
        onFileSelected(filename);
    }
    // AppState should have a method to handle file selection
    appState->selectFile(filename);
}

// Method to handle GUI refresh button press
void GUIContext::refreshRequested() {
    if (onRefreshRequested) {
        onRefreshRequested();
    }
    // AppState should have a method to refresh the view or data
    appState->refreshData();
}

// Method to handle slider value change in GUI
void GUIContext::sliderValueChanged(int value) {
    if (onSliderValueChanged) {
        onSliderValueChanged(value);
    }
    // AppState should have a method to update based on the slider's value
    appState->updateSliderValue(value);
}

// Method to handle serialization button press in GUI
void GUIContext::serializeButtonPressed() {
    if (onSerializeButtonPressed) {
        onSerializeButtonPressed();
    }
    // AppState should have a method to handle serialization
    appState->serializeData();
}
