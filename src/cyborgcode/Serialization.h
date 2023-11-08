#ifndef SERIALIZATION_H
#define SERIALIZATION_H

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "MyConfig.h"

// Eigen Serialization
bool serializeVector(const Eigen::VectorXd& vec, const std::string& filepath);
bool deserializeVector(Eigen::VectorXd& vec, const std::string& filepath);
bool serializeMatrix(const Eigen::MatrixXd& mat, const std::string& filepath);
bool deserializeMatrix(Eigen::MatrixXd& mat, const std::string& filepath);

// MyConfig Serialization
bool serializeConfig(const MyConfig& config, const std::string& filepath);
bool deserializeConfig(MyConfig& config, const std::string& filepath);

// CSV Serialization
bool writeCSV(const std::vector<double>& data, const std::string& filepath);
bool readCSV(std::vector<double>& data, const std::string& filepath);

#endif // SERIALIZATION_H
