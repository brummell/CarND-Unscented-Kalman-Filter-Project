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


// TODO: REFACTOR THIS!!!!! THIS MAKES SENSE TO ACTUALLY HAVE A GOOD IMPLEMENTATION WITH THEORETICAL AND PRACTICAL DOCS
// TODO: COULD EVEN EVENTUALLY BE A LIBRARY


UKF::UKF() {
    use_laser_ = true; // if this is false, laser measurements will be ignored (except during init)
    use_radar_ = true; // if this is false, radar measurements will be ignored (except during init)
    x_ = VectorXd(5); // initial state vector
    P_ = MatrixXd(5, 5); // initial covariance matrix
    std_a_ = 0.2; // Process noise standard deviation longitudinal acceleration in m/s^2
    std_yawdd_ = 0.2; // Process noise standard deviation yaw acceleration in rad/s^2
    std_laspx_ = 0.15; // Laser measurement noise standard deviation position1 in m
    std_laspy_ = 0.15; // Laser measurement noise standard deviation position2 in m
    std_radr_ = 0.3; // Radar measurement noise standard deviation radius in m
    std_radphi_ = 0.03; // Radar measurement noise standard deviation angle in rad
    std_radrd_ = 0.3; // Radar measurement noise standard deviation radius change in m/s
    n_x_ = 5;     ///* State dimension
    // TODO: Give dims to sensor modalities or just extract them.
    n_aug_ = 7;     ///* Augmented state dimension
    lambda_ = 3 - n_aug_;     ///* Sigma point spreading parameter
    NIS_radar_ = 0;    ///* the current NIS for radar
    NIS_laser_ = 0;     ///* the current NIS for laser
    previous_timestamp_ = 0;
    Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
    weights_ = VectorXd(2 * n_aug_ + 1); //create vector for weights;     ///* Weights of sigma points

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} measurement_pack The latest measurement data of
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
        P_ <<   1, 0, 0, 0, 0,
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
        UpdateRadar(measurement_pack.raw_measurements_);
        cout << "RADAR" << endl;
    } else {
        UpdateLidar(measurement_pack.raw_measurements_);
        cout << "LIDAR" << endl;
    }
}

void UKF::SetWeights() {
    double weight_0 = lambda_ / (lambda_ + n_aug_);
    weights_(0) = weight_0;
    for (int i = 1; i < 2 * n_aug_ + 1; i++) {
        double weight = 0.5 / (n_aug_ + lambda_);
        weights_(i) = weight;
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    // rather than generating a second set of sigma points to be able to propogate though the noise
    //   part of the process model, include them as dims in original sig point calc, model will squeeze back to 5d
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

    //predict sigma points
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
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

    SetWeights();

    //predicted state mean
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x_ = x_ + weights_(i) * Xsig_pred.col(i);
    }

    //predicted state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        VectorXd x_diff = Xsig_pred.col(i) - x_;         // state difference
        x_diff(3) = tools.NormalizeAngle(x_diff(3));  //angle normalization
        P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
    }
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const VectorXd &z) {
    // TODO: run NIS calculation

    //set measurement dimension, lidar measures (x,y) pairs
    int n_z = 2;

    SetWeights();

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    /*******************************************************************************
     * Convert predictions to measurement space
     ******************************************************************************/
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        double p_x = Xsig_pred(0, i);
        double p_y = Xsig_pred(1, i);

        // measurement model
        Zsig(0, i) = p_x;
        Zsig(1, i) = p_y;
    }

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i); // weights are for this
    }

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z, n_z);
    R << std_laspx_ * std_laspx_, 0,
            0, std_laspy_ * std_laspy_;
    S = S + R;

    /*******************************************************************************
     * Measurement Update
     ******************************************************************************/

    MatrixXd Tc = MatrixXd(n_x_, n_z);     //create matrix for cross correlation Tc

    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        VectorXd z_diff = Zsig.col(i) - z_pred; //residual
        VectorXd x_diff = Xsig_pred.col(i) - x_; // state difference
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    MatrixXd K = Tc * S.inverse(); // Kalman gain K;
    VectorXd z_diff = z - z_pred;     //residual

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();
}


/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const VectorXd &z) {
    // TODO: run NIS calculation

    int n_z = 3;  //  set measurement dimension, radar can measure r, phi, and r_dot

    // TODO: think this is sensor modality agnostic, extract
    SetWeights();

    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1); //create matrix for sigma points in measurement space


    /*******************************************************************************
     * Convert predictions to measurement space
     ******************************************************************************/
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        // extract values for better readibility
        double p_x = Xsig_pred(0, i);
        double p_y = Xsig_pred(1, i);
        double v = Xsig_pred(2, i);
        double yaw = Xsig_pred(3, i);

        double v1 = cos(yaw) * v;
        double v2 = sin(yaw) * v;

        // measurement model
        Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);                        //r
        Zsig(1, i) = atan2(p_y, p_x);                                 //phi
        Zsig(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y);   //r_dot
    }

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        // using weights instead of simple mean to potentially take advantage of distro knowledge
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        VectorXd z_diff = Zsig.col(i) - z_pred;   //residual
        z_diff(1) = tools.NormalizeAngle(z_diff(1)); //angle normalization
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z, n_z);
    R << std_radr_ * std_radr_, 0, 0,
            0, std_radphi_ * std_radphi_, 0,
            0, 0, std_radrd_ * std_radrd_;
    S = S + R;

   /*******************************************************************************
    * Measurement Update
    ******************************************************************************/
    //calculate cross correlation matrix Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        VectorXd z_diff = Zsig.col(i) - z_pred; //residual
        z_diff(1) = tools.NormalizeAngle(z_diff(1)); //angle normalization
        VectorXd x_diff = Xsig_pred.col(i) - x_;  // state difference
        x_diff(3) = tools.NormalizeAngle(x_diff(3)); //angle normalization;
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    MatrixXd K = Tc * S.inverse(); //Kalman gain K;
    VectorXd z_diff = z - z_pred; //residual

    z_diff(1) = tools.NormalizeAngle(z_diff(1)); //angle normalization

    // update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();
}


#pragma clang diagnostic pop