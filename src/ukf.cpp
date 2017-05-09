#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

#pragma clang diagnostic push
#pragma ide diagnostic ignored "IncompatibleTypes"
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    use_laser_ = true; // if this is false, laser measurements will be ignored (except during init)
    use_radar_ = true; // if this is false, radar measurements will be ignored (except during init)
    x_ = VectorXd(5); // initial state vector
    P_ = MatrixXd(5, 5); // initial covariance matrix
    // TODO: These are way off (this was when they started at 30.
    std_a_ = 0.2; // Process noise standard deviation longitudinal acceleration in m/s^2
    std_yawdd_ = 0.2; // Process noise standard deviation yaw acceleration in rad/s^2
    std_laspx_ = 0.15; // Laser measurement noise standard deviation position1 in m
    std_laspy_ = 0.15; // Laser measurement noise standard deviation position2 in m
    std_radr_ = 0.3; // Radar measurement noise standard deviation radius in m
    std_radphi_ = 0.03; // Radar measurement noise standard deviation angle in rad
    std_radrd_ = 0.3; // Radar measurement noise standard deviation radius change in m/s
    VectorXd weights_;     ///* Weights of sigma points
    n_x_ = 5;     ///* State dimension
    n_aug_ = 7;     ///* Augmented state dimension
    lambda_ = 3 - n_aug_;     ///* Sigma point spreading parameter
    NIS_radar_ = 0;    ///* the current NIS for radar
    NIS_laser_ = 0;     ///* the current NIS for laser
    previous_timestamp_ = 0;
    H_ = MatrixXd(2, 4);
    H_laser_ = MatrixXd(2, 4);
    H_laser_ << 1, 0, 0, 0,
            0, 1, 0, 0;
    Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);



}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} measurement_packa The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
    /*****************************************************************************
    *  Initialization
    ****************************************************************************/
    if (!is_initialized_) {
        x_ = VectorXd(5);
        x_ << 1, 1, 1, 1, 1;

        P_ = MatrixXd(5, 5);
        // TODO: this may need to be changed
        P_ << 1, 0, 0, 0, 0,
                0, 1, 0, 0, 0,
                0, 0, 1, 0, 0,
                0, 0, 0, 1, 0,
                0, 0, 0, 0, 1;

        previous_timestamp_ = measurement_pack.timestamp_;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            // Convert radar from polar to cartesian coordinates and initialize state.
            auto rho = measurement_pack.raw_measurements_(0);
            auto phi = measurement_pack.raw_measurements_(1);
            auto rho_dot = measurement_pack.raw_measurements_(2);

            x_ << rho * cos(phi),
                    rho * sin(phi),
                    rho_dot, // seems a slight bit better initialization than 0s at basically no cost.
                    phi,
                    0;

        } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            // Initialize state.
            x_ << measurement_pack.raw_measurements_(0),
                    measurement_pack.raw_measurements_(1),
                    0,
                    0,
                    0;
        }
        is_initialized_ = true;
        return; // TODO: rebuild flow with else
    }

    /*****************************************************************************
    *  Prediction
    ****************************************************************************/
    //compute the time elapsed between the current and previous measurements
    auto dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;    //dt - expressed in seconds
    previous_timestamp_ = measurement_pack.timestamp_;

    Prediction(dt);

    /*****************************************************************************
     *  Update
     ****************************************************************************/
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
//        H_ = Hj_;
//        R_ = R_radar_;
        UpdateRadar(measurement_pack.raw_measurements_);
        cout << "RADAR" << endl;
    } else {
        H_ = H_laser_;
//        ekf_.R_ = R_laser_;
        UpdateLidar(measurement_pack.raw_measurements_);
        cout << "LIDAR" << endl;


    }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    /**
    TODO:

    Complete this function! Estimate the object's location. Modify the state
    vector, x_. Predict sigma points, the state, and the state covariance matrix.
    */
    // TODO: inject process covariance terms into mean estimate vector
    // TODO: use formulae to create set of sigma points, X, from current distribution (d=7x1)
    // TODO: propagate each point through the process model (CRTV)
    // TODO: update augmented covariance matrix (top left plus process noise in bottom right)

    VectorXd x_aug = VectorXd(7);    //create augmented mean vector
    MatrixXd P_aug = MatrixXd(7, 7);    //create augmented state covariance
    MatrixXd X_sigma_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);     //create sigma point matrix

    //create augmented mean state
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    //create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5, 5) = P_;
    P_aug(5, 5) = std_a_ * std_a_;
    P_aug(6, 6) = std_yawdd_ * std_yawdd_;

    MatrixXd L = P_aug.llt().matrixL();     // Use Cholesky decomposition to find sqrt of process covariance matri

    //create augmented sigma points
    X_sigma_aug.col(0) = x_aug;
    for (int i = 0; i < n_aug_; i++) {
        X_sigma_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
        X_sigma_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
    }

//    MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

    //predict sigma points
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        //extract values for better readability
        double p_x = X_sigma_aug(0, i);
        double p_y = X_sigma_aug(1, i);
        double v = X_sigma_aug(2, i);
        double yaw = X_sigma_aug(3, i);
        double yawd = X_sigma_aug(4, i);
        double nu_a = X_sigma_aug(5, i);
        double nu_yawdd = X_sigma_aug(6, i);

        //predicted state values
        double px_p, py_p;

        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
            py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
        } else {
            px_p = p_x + v * delta_t * cos(yaw);
            py_p = p_y + v * delta_t * sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd * delta_t;
        double yawd_p = yawd;

        //add noise
        px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
        py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
        v_p = v_p + nu_a * delta_t;

        yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
        yawd_p = yawd_p + nu_yawdd * delta_t;

        //write predicted sigma point into right column
        Xsig_pred(0, i) = px_p;
        Xsig_pred(1, i) = py_p;
        Xsig_pred(2, i) = v_p;
        Xsig_pred(3, i) = yaw_p;
        Xsig_pred(4, i) = yawd_p;
    }

    VectorXd weights = VectorXd(2 * n_aug_ + 1); //create vector for weights
    VectorXd x = VectorXd(n_x_); //create vector for predicted state
    MatrixXd P = MatrixXd(n_x_, n_x_); //create covariance matrix for prediction

    // set weights
    double weight_0 = lambda_ / (lambda_ + n_aug_);
    weights(0) = weight_0;
    for (int i = 1; i < 2 * n_aug_ + 1; i++) {  //2n+1 weights
        double weight = 0.5 / (n_aug_ + lambda_);
        weights(i) = weight;
    }

    //predicted state mean
    // TODO: make this work for in-place case
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x_ = x_ + weights(i) * Xsig_pred.col(i);
    }

    //predicted state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        VectorXd x_diff = Xsig_pred.col(i) - x;         // state difference
        while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;         //angle normalization
        while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;
        P_ = P_ + weights(i) * x_diff * x_diff.transpose();
    }
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const VectorXd &z) {
    /**
    TODO:

    Complete this function! Use lidar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the lidar NIS.
    */
    // TODO: update occurs via regular KF methods (I believe?), can probably copy from EKFs
    // TODO: run NIS calculation
//    MatrixXd Ht_(4,2);
//    Ht_ = H_.transpose();
//    VectorXd y = z - H_ * x_;
//    MatrixXd S_(2,2);
//    S_ = H_ * P_ * Ht_;
////    S_ = H_ * P_ * Ht_ + R_;
//    MatrixXd S_i(2,2);
//    S_i = S_.inverse();
//    MatrixXd K_(4,2);
//    K_ = P_ * Ht_ * S_i;
//    MatrixXd I_(4, 4);
//    I_.setIdentity(4,4);
//    x_ = x_ + K_ * y;
//    P_ = (I_ - K_ * H_) * P_;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const VectorXd &z) { // TODO: This is definitely fucked up!!!
    /**
    TODO:

    Complete this function! Use radar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the radar NIS.
    */
    // TODO: reuse sigma points from prediction step partial output, run them through non-linear measurement functions
    // TODO: use modified Kalman gain, continue using same old stuff from EKF
    // TODO: update covariance with additive noise R
    // TODO: run NIS calculation

    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;

    //set vector for weights
    VectorXd weights = VectorXd(2 * n_aug_ + 1);
    double weight_0 = lambda_ / (lambda_ + n_aug_);
    weights(0) = weight_0;
    for (int i = 1; i < 2 * n_aug_ + 1; i++) {  //2n+1 weights
        double weight = 0.5 / (n_aug_ + lambda_);
        weights(i) = weight;
    }

    //create example matrix with sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1); // d = 3x15
    VectorXd z_pred = VectorXd(n_z); // d = 3x1     //create example vector for mean predicted measurement
    MatrixXd S = MatrixXd(n_z, n_z);  // d = 3x3   //create example matrix for predicted measurement covariance

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        VectorXd z_diff = Zsig.col(i) - z_pred;         //residual


        // TODO: THESE HAVE GOTTA GO!!!! ATAN2
        //angle normalization
        while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
        while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

        VectorXd x_diff = Xsig_pred.col(i) - x_;         // state difference

        //angle normalization
        while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
        while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

        Tc = Tc + weights(i) * x_diff * z_diff.transpose();
    }

    MatrixXd K = Tc * S.inverse();     //Kalman gain K;
    VectorXd z_diff = z - z_pred;     //residual

    //angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();
}


#pragma clang diagnostic pop