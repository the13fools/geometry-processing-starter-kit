#include "AppState.h"
#include "FileParser.h"
#include "Serialization.h"
#include <iostream>

// Implement the refreshFileLists to populate the lists of files
void AppState::refreshFileLists() {
    FileParser fileParser(directoryPath);

    // Refresh bfraFiles and bmomFiles
    fileParser.parseBfraFiles(bfraFiles);
    fileParser.parseBmomFiles(bmomFiles);

    // Refresh objFilePath
    fileParser.parseObjFile(objFilePath);

    // Additional logic to handle other types of files if needed
}

// Method to log current variables (if logging is enabled)
void AppState::logCurrentVariables() {
    if (logVariables) {
        // Implement serialization functions to log variables to a file
        std::string logFilePath = directoryPath + "/log.json";

        // Serialize variables to log file
        if (Serialization::serializeAppState(*this, logFilePath)) {
            std::cout << "Logged variables to: " << logFilePath << std::endl;
        } else {
            std::cerr << "Failed to log variables." << std::endl;
        }
    }
}

// Method to update the current frame and moment indices
void AppState::updateIndices(int frameIndex, int momentIndex) {
    currentFrameIndex = frameIndex;
    currentMomentIndex = momentIndex;
    // Additional logic to handle index updates
}

// Method to serialize the current app state to a file
bool AppState::serializeData() {
    std::string dataFilePath = directoryPath + "/appState.json";

    // Serialize app state to a JSON file
    if (Serialization::serializeAppState(*this, dataFilePath)) {
        std::cout << "Serialized app state to: " << dataFilePath << std::endl;
        return true;
    } else {
        std::cerr << "Failed to serialize app state." << std::endl;
        return false;
    }
}

// Method to deserialize app state from a file
bool AppState::deserializeData() {
    std::string dataFilePath = directoryPath + "/appState.json";

    // Deserialize app state from a JSON file
    if (Serialization::deserializeAppState(*this, dataFilePath)) {
        std::cout << "Deserialized app state from: " << dataFilePath << std::endl;
        return true;
    } else {
        std::cerr << "Failed to deserialize app state." << std::endl;
        return false;
    }
}

// Additional AppState-related functionality can be added here as needed.
