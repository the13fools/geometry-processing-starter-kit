#ifndef FILE_PARSER_H
#define FILE_PARSER_H

#include <string>
#include <vector>
#include <optional>

/**
 * The FileParser class is responsible for parsing files of different types
 * within a given directory and retrieving relevant data.
 */
class FileParser {
public:
    // Constructor with the directory to scan
    explicit FileParser(const std::string& directoryPath);

    // Parse a file with a given ID
    bool parseFileWithID(int fileId);

    // Retrieve the largest file's data
    bool parseLargestFile();

    // Load additional .obj file if found
    bool parseObjFile();

    // Get the path to the largest bfra file
    const std::string& getLargestBfraFile() const { return largestBfraFile; }

    // Get the path to the largest bmom file
    const std::string& getLargestBmomFile() const { return largestBmomFile; }

private:
    std::string directoryPath;
    std::vector<std::string> bfraFiles;
    std::vector<std::string> bmomFiles;
    std::string objFilePath; // Path to the .obj file if found
    std::string largestBfraFile; // Path to the largest bfra file
    std::string largestBmomFile; // Path to the largest bmom file

    // Helper functions to find and sort files, populate the file lists
    void scanDirectory();
    void findLargestIDFile();
};

#endif // FILE_PARSER_H
