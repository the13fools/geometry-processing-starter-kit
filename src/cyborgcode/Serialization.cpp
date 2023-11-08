#include "Serialization.h"
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

// Serialize Eigen vector to a binary file
bool serializeVector(const Eigen::VectorXd& vec, const std::string& filepath) {
    try {
        std::ofstream outFile(filepath, std::ios::binary);
        if (!outFile.is_open()) {
            std::cerr << "Error: Unable to open file for writing: " << filepath << std::endl;
            return false;
        }
        
        // Write vector size
        int size = static_cast<int>(vec.size());
        outFile.write(reinterpret_cast<char*>(&size), sizeof(int));

        // Write vector data
        outFile.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(double));
        outFile.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to serialize vector: " << e.what() << std::endl;
        return false;
    }
}

// Deserialize Eigen vector from a binary file
bool deserializeVector(Eigen::VectorXd& vec, const std::string& filepath) {
    try {
        std::ifstream inFile(filepath, std::ios::binary);
        if (!inFile.is_open()) {
            std::cerr << "Error: Unable to open file for reading: " << filepath << std::endl;
            return false;
        }

        // Read vector size
        int size;
        inFile.read(reinterpret_cast<char*>(&size), sizeof(int));

        // Read vector data
        vec.resize(size);
        inFile.read(reinterpret_cast<char*>(vec.data()), size * sizeof(double));
        inFile.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to deserialize vector: " << e.what() << std::endl;
        return false;
    }
}

// Serialize Eigen matrix to a binary file
bool serializeMatrix(const Eigen::MatrixXd& mat, const std::string& filepath) {
    try {
        std::ofstream outFile(filepath, std::ios::binary);
        if (!outFile.is_open()) {
            std::cerr << "Error: Unable to open file for writing: " << filepath << std::endl;
            return false;
        }
        
        // Write matrix rows and cols
        int rows = static_cast<int>(mat.rows());
        int cols = static_cast<int>(mat.cols());
        outFile.write(reinterpret_cast<char*>(&rows), sizeof(int));
        outFile.write(reinterpret_cast<char*>(&cols), sizeof(int));

        // Write matrix data
        outFile.write(reinterpret_cast<const char*>(mat.data()), rows * cols * sizeof(double));
        outFile.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to serialize matrix: " << e.what() << std::endl;
        return false;
    }
}

// Deserialize Eigen matrix from a binary file
bool deserializeMatrix(Eigen::MatrixXd& mat, const std::string& filepath) {
    try {
        std::ifstream inFile(filepath, std::ios::binary);
        if (!inFile.is_open()) {
            std::cerr << "Error: Unable to open file for reading: " << filepath << std::endl;
            return false;
        }

        // Read matrix rows and cols
        int rows, cols;
        inFile.read(reinterpret_cast<char*>(&rows), sizeof(int));
        inFile.read(reinterpret_cast<char*>(&cols), sizeof(int));

        // Read matrix data
        mat.resize(rows, cols);
        inFile.read(reinterpret_cast<char*>(mat.data()), rows * cols * sizeof(double));
        inFile.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to deserialize matrix: " << e.what() << std::endl;
        return false;
    }
}

// Serialize MyConfig struct to a JSON file
bool serializeConfig(const MyConfig& config, const std::string& filepath) {
    try {
        json j;
        j["w_bound"] = config.w_bound;
        j["w_smooth"] = config.w_smooth;
        j["w_smooth_vector"] = config.w_smooth_vector;
        j["w_curl"] = config.w_curl;
        j["w_s_perp"] = config.w_s_perp;
        j["w_attenuate"] = config.w_attenuate;
        j["convergence_eps"] = config.convergence_eps;
        j["identity_weight"] = config.identity_weight;
        j["prev_energy"] = config.prev_energy;
        j["useProjHessian"] = config.useProjHessian;

        std::ofstream outFile(filepath);
        if (!outFile.is_open()) {
            std::cerr << "Error: Unable to open file for writing JSON: " << filepath << std::endl;
            return false;
        }

        outFile << j.dump(4); // Pretty-print with 4 spaces of indentation
        outFile.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to serialize config to JSON: " << e.what() << std::endl;
        return false;
    }
}

// Deserialize MyConfig struct from a JSON file
bool deserializeConfig(MyConfig& config, const std::string& filepath) {
    try {
        std::ifstream inFile(filepath);
        if (!inFile.is_open()) {
            std::cerr << "Error: Unable to open file for deserialization: " << filepath << std::endl;
            return false;
        }

        json j;
        inFile >> j;

        config.w_bound = j.at("w_bound");
        config.w_smooth = j.at("w_smooth");
        config.w_smooth_vector = j.at("w_smooth_vector");
        config.w_curl = j.at("w_curl");
        config.w_s_perp = j.at("w_s_perp");
        config.w_attenuate = j.at("w_attenuate");
        config.convergence_eps = j.at("convergence_eps");
        config.identity_weight = j.at("identity_weight");
        config.prev_energy = j.at("prev_energy");
        config.useProjHessian = j.at("useProjHessian");

        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to deserialize config from JSON: " << e.what() << std::endl;
        return false;
    }
}

// Serialize a vector of doubles to a CSV file
bool writeCSV(const std::vector<double>& data, const std::string& filepath) {
    try {
        std::ofstream outFile(filepath);
        if (!outFile.is_open()) {
            std::cerr << "Error: Unable to open file for writing CSV: " << filepath << std::endl;
            return false;
        }

        for (const double& value : data) {
            outFile << value << "\n";
        }

        outFile.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to write data to CSV: " << e.what() << std::endl;
        return false;
    }
}

// Deserialize a vector of doubles from a CSV file
bool readCSV(std::vector<double>& data, const std::string& filepath) {
    try {
        std::ifstream inFile(filepath);
        if (!inFile.is_open()) {
            std::cerr << "Error: Unable to open file for reading CSV: " << filepath << std::endl;
            return false;
        }

        data.clear(); // Clear the existing data

        double value;
        while (inFile >> value) {
            data.push_back(value);
        }

        inFile.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to read data from CSV: " << e.what() << std::endl;
        return false;
    }
}

// Serialize JSON data to a JSON file (using nlohmann/json)
bool serializeJSON(const nlohmann::json& jsonData, const std::string& filepath) {
    try {
        std::ofstream outFile(filepath);
        if (!outFile.is_open()) {
            std::cerr << "Error: Unable to open file for writing JSON: " << filepath << std::endl;
            return false;
        }

        outFile << jsonData.dump(4); // Pretty-print with 4 spaces of indentation
        outFile.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to serialize JSON data to file: " << e.what() << std::endl;
        return false;
    }
}

// Deserialize JSON data from a JSON file (using nlohmann/json)
bool deserializeJSON(nlohmann::json& jsonData, const std::string& filepath) {
    try {
        std::ifstream inFile(filepath);
        if (!inFile.is_open()) {
            std::cerr << "Error: Unable to open file for reading JSON: " << filepath << std::endl;
            return false;
        }

        inFile >> jsonData;
        inFile.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to deserialize JSON data from file: " << e.what() << std::endl;
        return false;
    }
}
