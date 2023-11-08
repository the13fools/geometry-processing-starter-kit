#include "AppState.h"
#include "DirectoryScanner.h"

// MIT License
// This code is released by a user conversing with a chatbot.

// Implement the refreshFileLists to populate the lists of files
void AppState::refreshFileLists() {
    // Assuming DirectoryScanner can list files and find the first matching file
    DirectoryScanner scanner(directoryPath);
    bfraFiles = scanner.getFileListWithExtension(".bfra");
    bmomFiles = scanner.getFileListWithExtension(".bmom");
    objFilePath = scanner.getFirstFileWithExtension(".obj");
}

// Additional AppState-related functionality can be added here as needed.
