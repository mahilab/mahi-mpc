#pragma once

#include <Eigen/Dense>
#include <vector>

class LinearModel
{
public:
    virtual Eigen::MatrixXd get_A(const std::vector<double>& x, const std::vector<double>& control) = 0;
    virtual Eigen::MatrixXd get_B(const std::vector<double>& x, const std::vector<double>& control) = 0;
    virtual Eigen::MatrixXd get_Xdot(const std::vector<double>& x, const std::vector<double>& control) = 0;
    static std::vector<double> flatten_matrix(Eigen::MatrixXd matrix);
    LinearModel(/* args */);
};