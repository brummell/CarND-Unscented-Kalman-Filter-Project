#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    /**
    TODO:
      * Calculate the RMSE here.
    */
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;

    // should consider throwing error here and for div by 0 in jacobian
    if (estimations.size() != ground_truth.size()
        || estimations.size() == 0) {
        std::cout << "Invalid estimation or ground_truth data" << std::endl;
        return rmse;
    }

    // convert to modern syntax
    for (unsigned int i = 0; i < estimations.size(); ++i) {
        VectorXd residual = estimations[i] - ground_truth[i];
        residual = residual.array() * residual.array(); //coefficient-wise multiplication
        rmse += residual;
    }

    //calculate the mean
    rmse = rmse / estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;
}

double Tools::NormalizeAngle(const double in_angle) {
    // TODO: This is not an ideal implementation (atan2), but it'll do for now
    auto angle = in_angle;
    while (angle > M_PI) angle -= 2. * M_PI;
    while (angle < -M_PI) angle += 2. * M_PI;
    return angle;
}