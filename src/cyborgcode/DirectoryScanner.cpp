#include "DirectoryScanner.h"
#include <filesystem>
#include <algorithm>

namespace fs = std::filesystem;

DirectoryScanner::DirectoryScanner(const std::string& path)
    : directoryPath(path) {}

std::vector<std::string> DirectoryScanner::getFileListWithExtension(const std::string& extension) const {
    std::vector<std::string> fileList;
    for (const auto& entry : fs::directory_iterator(directoryPath)) {
        if (entry.is_regular_file() && entry.path().extension() == extension) {
            fileList.push_back(entry.path().filename().string());
        }
    }
    // Sort the file list to ensure consistent order
    std::sort(fileList.begin(), fileList.end());
    return fileList;
}

std::optional<std::string> DirectoryScanner::getFirstFileWithExtension(const std::string& extension) const {
    for (const auto& entry : fs::directory_iterator(directoryPath)) {
        if (entry.is_regular_file() && entry.path().extension() == extension) {
            return entry.path().filename().string();
        }
    }
    return {}; // Return empty optional if no file is found
}
