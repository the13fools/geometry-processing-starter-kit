#ifndef DIRECTORYSCANNER_H
#define DIRECTORYSCANNER_H

#include <string>
#include <vector>
#include <optional>

// MIT License
// This code is released by a user conversing with a chatbot.

class DirectoryScanner {
public:
    DirectoryScanner(const std::string& path);

    std::vector<std::string> getFileListWithExtension(const std::string& extension) const;
    std::optional<std::string> getFirstFileWithExtension(const std::string& extension) const;

private:
    std::string directoryPath;
};

#endif // DIRECTORYSCANNER_H
