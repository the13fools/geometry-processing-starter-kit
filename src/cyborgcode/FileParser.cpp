#include "FileParser.h"
#include <filesystem>
#include <algorithm>
#include <regex>

namespace fs = std::filesystem;

FileParser::FileParser(const std::string& directoryPath)
    : directoryPath(directoryPath) {
    scanDirectory();
}

void FileParser::scanDirectory() {
    for (const auto& entry : fs::directory_iterator(directoryPath)) {
        if (entry.is_regular_file()) {
            std::string filename = entry.path().filename().string();
            if (std::regex_match(filename, std::regex("primal_\\d{5}\\.bfra$"))) {
                bfraFiles.push_back(entry.path().string());
            } else if (std::regex_match(filename, std::regex("optvars_\\d{5}\\.bmom$"))) {
                bmomFiles.push_back(entry.path().string());
            } else if (filename.ends_with(".obj") && objFilePath.empty()) {
                objFilePath = entry.path().string(); // Assuming only one .obj file
            }
        }
    }

    findLargestIDFile();
}

void FileParser::findLargestIDFile() {
    auto fileIdComparator = [](const std::string& file1, const std::string& file2) {
        int id1 = std::stoi(file1.substr(file1.find_last_of('_') + 1, 5));
        int id2 = std::stoi(file2.substr(file2.find_last_of('_') + 1, 5));
        return id1 > id2;
    };

    if (!bfraFiles.empty()) {
        std::sort(bfraFiles.begin(), bfraFiles.end(), fileIdComparator);
        largestBfraFile = bfraFiles.back();
    }
    if (!bmomFiles.empty()) {
        std::sort(bmomFiles.begin(), bmomFiles.end(), fileIdComparator);
        largestBmomFile = bmomFiles.back();
    }
}

// Rest of FileParser implementation...
