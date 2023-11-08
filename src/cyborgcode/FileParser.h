#ifndef FILE_PARSER_H
#define FILE_PARSER_H

#include <string>
#include <vector>
#include <optional>
#include <Eigen/Dense>
#include "Serialization.h" // Include the Serialization header

// Define the FileType enum
enum class FileType {
    BFRA,
    BMOM,
    OBJ,
    FRA,
    MOM,
};

/**
 * The FileParser class is responsible for parsing files of different types
 * within a given directory and retrieving relevant data.
 */
class FileParser {
public:
    // Constructor with the directory to scan
    explicit FileParser(const std::string& directoryPath);

    // Parse a file with a given ID
    bool parseFileWithID(Eigen::VectorXd& data, FileType fileType, int fileId);

    // Retrieve the largest file's data
    bool parseLargestFile(Eigen::VectorXd& data, FileType fileType);

    // Load an .obj file if found
    bool parseObjFile(Eigen::MatrixXd& V, Eigen::MatrixXi& F);

    // Update the current directory path
    void setDirectoryPath(const std::string& directoryPath);

private:
    std::string directoryPath;
    std::vector<std::string> bfraFiles;
    std::vector<std::string> bmomFiles;
    std::string objFilePath;
    std::string largestBfraFile;
    std::string largestBmomFile;

    // Helper functions to scan and sort files
    void scanDirectory();
    void findLargestIDFile();
};

#endif // FILE_PARSER_H
