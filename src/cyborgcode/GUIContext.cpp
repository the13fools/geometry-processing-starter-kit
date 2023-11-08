#include "GUIContext.h"

// Define the GUI interaction methods
void GUIContext::fileSelected(const std::string& filename) {
    if (onFileSelected) {
        onFileSelected(filename);
    }
    // AppState should have a method to handle file selection
    appState->selectFile(filename);
}

void GUIContext::refreshRequested() {
    if (onRefreshRequested) {
        onRefreshRequested();
    }
    // AppState should have a method to refresh the view or data
    appState->refreshData();
}

void GUIContext::sliderValueChanged(int value) {
    if (onSliderValueChanged) {
        onSliderValueChanged(value);
    }
    // AppState should have a method to update based on the slider's value
    appState->updateSliderValue(value);
}

void GUIContext::serializeButtonPressed() {
    if (onSerializeButtonPressed) {
        onSerializeButtonPressed();
    }
    // AppState should have a method to handle serialization
    appState->serializeData();
}

// Other GUIContext methods can be implemented below
