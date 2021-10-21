#pragma once

#include <Mahi/Casadi/LinearModel.hpp>

LinearModel::LinearModel(/* args */){

}

std::vector<double> LinearModel::flatten_matrix(Eigen::MatrixXd matrix){
    std::vector<double> flattened_matrix;
    flattened_matrix.reserve(matrix.size());
    for (auto i = 0; i < matrix.cols(); i++){
        for (size_t j = 0; j < matrix.rows(); j++){
            flattened_matrix.push_back(matrix(j,i));
        }
    }
    return flattened_matrix;
}