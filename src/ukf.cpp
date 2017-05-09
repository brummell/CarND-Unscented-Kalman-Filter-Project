#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);

    // TODO: These are way off (this was when they started at 30.)
    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 0.2;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.2;

    // Laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15;

    // Laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;

    // Radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;

    // Radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.03;

    // Radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.3;

    ///* Weights of sigma points
    VectorXd weights_;

    ///* State dimension
    n_x_ = 5;

    ///* Augmented state dimension
    n_aug_ = 7;

    // TODO: this needs work
    ///* Sigma point spreading parameter
    lambda_ = 3 - n_aug_;

    ///* the current NIS for radar
    NIS_radar_ = 0;

    ///* the current NIS for laser
    NIS_laser_ = 0;

    previous_timestamp_ = 0;


    /**
    Hint: one or more values initialized above might be wildly off...
    */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
    /**
    TODO:

    Complete this function! Make sure you switch between lidar and radar
    measurements.
    */
    // TODO: probably can be more or less taken from EKF project
    /*****************************************************************************
 *  Initialization
 ****************************************************************************/
    if (!is_initialized_) {
        // first measurement
        x_ = VectorXd(5);
        x_ << 1, 1, 1, 1, 1;

        P_ = MatrixXd(5, 5);
        // TODO: this may need to be changed
        P_ << 10000, 0, 0, 0, 0,
                0, 10000, 0, 0, 0,
                0, 0, 10000, 0, 0,
                0, 0, 0, 10000, 0,
                0, 0, 0, 0, 10000;

        Q_ = MatrixXd(5, 5);

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
        return;
    }

    /*****************************************************************************
    *  Prediction
    ****************************************************************************/
    //compute the time elapsed between the current and previous measurements
    auto dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;    //dt - expressed in seconds
    previous_timestamp_ = measurement_pack.timestamp_;

    //Modify the F matrix so that the time is integrated
    F_(0, 2) = dt;
    F_(1, 3) = dt;

    //set the process covariance matrix Q
    Q_ << pow(dt, 4.0) / 4 * noise_ax, 0, pow(dt, 3.0) / 2 * noise_ax, 0,
            0, pow(dt, 4.0) / 4 * noise_ay, 0, pow(dt, 3.0) / 2 * noise_ay,
            pow(dt, 3.0) / 2 * noise_ax, 0, pow(dt, 2.0) * noise_ax, 0,
            0, pow(dt, 3.0) / 2 * noise_ay, 0, pow(dt, 2.0) * noise_ay;


    Prediction(dt);

    /*****************************************************************************
     *  Update
     ****************************************************************************/
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        H_ = Hj_;
        R_ = R_radar_;
        UpdateRadar(measurement_pack.raw_measurements_);
        cout << "RADAR" << endl;
    } else {
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;
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


    MatrixXd X_sigma = MatrixXd(n_x_, 2 * n_x_ + 1);    // Create the sigma point matrix


    MatrixXd A = P_.llt().matrixL();    // Use Cholesky decomposition to find sqrt of process covariance matrix


    // Set mean sigma point, then loop, setting the rest
    X_sigma.col(0) = x_;
    for (int i = 0; i < n_x_; i++) {
        X_sigma.col(i + 1) = x_ + sqrt(lambda + n_x) * A.col(i);
        X_sigma.col(i + 1 + n_x_) = x_ - sqrt(lambda + n_x) * A.col(i);
    }

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

    MatrixXd L = P_aug.llt().matrixL();     //create square root matrix

    //create augmented sigma points
    X_sigma_aug.col(0) = x_aug;
    for (int i = 0; i < n_aug_; i++) {
        X_sigma_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
        X_sigma_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
    }

    MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

    //predict sigma points
    for (int i = 0; i< 2*n_aug_+1; i++)
    {
        //extract values for better readability
        double p_x = X_sigma_aug(0,i);
        double p_y = X_sigma_aug(1,i);
        double v = X_sigma_aug(2,i);
        double yaw = X_sigma_aug(3,i);
        double yawd = X_sigma_aug(4,i);
        double nu_a = X_sigma_aug(5,i);
        double nu_yawdd = X_sigma_aug(6,i);

        //predicted state values
        double px_p, py_p;

        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;

        //add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;

        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;

        //write predicted sigma point into right column
        Xsig_pred(0,i) = px_p;
        Xsig_pred(1,i) = py_p;
        Xsig_pred(2,i) = v_p;
        Xsig_pred(3,i) = yaw_p;
        Xsig_pred(4,i) = yawd_p;
    }

    VectorXd weights = VectorXd(2*n_aug_+1); //create vector for weights
    VectorXd x = VectorXd(n_x_); //create vector for predicted state
    MatrixXd P = MatrixXd(n_x_, n_x_); //create covariance matrix for prediction

    // set weights
    double weight_0 = lambda_/(lambda_ + n_aug_);
    weights(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
        double weight = 0.5/(n_aug_+lambda_);
        weights(i) = weight;
    }

    //predicted state mean
    // TODO: make this work for in-place case
    x.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x = x+ weights(i) * Xsig_pred.col(i);
    }

    //predicted state covariance matrix
    P.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

        VectorXd x_diff = Xsig_pred.col(i) - x;         // state difference
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;         //angle normalization
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        P = P + weights(i) * x_diff * x_diff.transpose() ;
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Use lidar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the lidar NIS.
    */
    // TODO: update occurs via regular KF methods (I believe?), can probably copy from EKFs
    // TODO: run NIS calculation
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
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
    void UKF::UpdateState(VectorXd *x_out, MatrixXd *P_out) {

        //set state dimension
        int n_x = 5;

        //set augmented dimension
        int n_aug = 7;

        //set measurement dimension, radar can measure r, phi, and r_dot
        int n_z = 3;

        //define spreading parameter
        double lambda = 3 - n_aug;

        //set vector for weights
        VectorXd weights = VectorXd(2 * n_aug + 1);
        double weight_0 = lambda / (lambda + n_aug);
        weights(0) = weight_0;
        for (int i = 1; i < 2 * n_aug + 1; i++) {  //2n+1 weights
            double weight = 0.5 / (n_aug + lambda);
            weights(i) = weight;
        }

        //create example matrix with predicted sigma points
        MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
        Xsig_pred <<
                  5.9374, 6.0640, 5.925, 5.9436, 5.9266, 5.9374, 5.9389, 5.9374, 5.8106, 5.9457, 5.9310, 5.9465, 5.9374, 5.9359, 5.93744,
                1.48, 1.4436, 1.660, 1.4934, 1.5036, 1.48, 1.4868, 1.48, 1.5271, 1.3104, 1.4787, 1.4674, 1.48, 1.4851, 1.486,
                2.204, 2.2841, 2.2455, 2.2958, 2.204, 2.204, 2.2395, 2.204, 2.1256, 2.1642, 2.1139, 2.204, 2.204, 2.1702, 2.2049,
                0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337, 0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188, 0.5367, 0.535048,
                0.352, 0.29997, 0.46212, 0.37633, 0.4841, 0.41872, 0.352, 0.38744, 0.40562, 0.24347, 0.32926, 0.2214, 0.28687, 0.352, 0.318159;

        //create example vector for predicted state mean
        VectorXd x = VectorXd(n_x);
        x <<
          5.93637,
                1.49035,
                2.20528,
                0.536853,
                0.353577;

        //create example matrix for predicted state covariance
        MatrixXd P = MatrixXd(n_x, n_x);
        P <<
          0.0054342, -0.002405, 0.0034157, -0.0034819, -0.00299378,
                -0.002405, 0.01084, 0.001492, 0.0098018, 0.00791091,
                0.0034157, 0.001492, 0.0058012, 0.00077863, 0.000792973,
                -0.0034819, 0.0098018, 0.00077863, 0.011923, 0.0112491,
                -0.0029937, 0.0079109, 0.00079297, 0.011249, 0.0126972;

        //create example matrix with sigma points in measurement space
        MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
        Zsig <<
             6.1190, 6.2334, 6.1531, 6.1283, 6.1143, 6.1190, 6.1221, 6.1190, 6.0079, 6.0883, 6.1125, 6.1248, 6.1190, 6.1188, 6.12057,
                0.24428, 0.2337, 0.27316, 0.24616, 0.24846, 0.24428, 0.24530, 0.24428, 0.25700, 0.21692, 0.24433, 0.24193, 0.24428, 0.24515, 0.245239,
                2.1104, 2.2188, 2.0639, 2.187, 2.0341, 2.1061, 2.1450, 2.1092, 2.0016, 2.129, 2.0346, 2.1651, 2.1145, 2.0786, 2.11295;

        //create example vector for mean predicted measurement
        VectorXd z_pred = VectorXd(n_z);
        z_pred <<
               6.12155,
                0.245993,
                2.10313;

        //create example matrix for predicted measurement covariance
        MatrixXd S = MatrixXd(n_z, n_z);
        S <<
          0.0946171, -0.000139448, 0.00407016,
                -0.000139448, 0.000617548, -0.000770652,
                0.00407016, -0.000770652, 0.0180917;

        //create example vector for incoming radar measurement
        VectorXd z = VectorXd(n_z);
        z <<
          5.9214,
                0.2187,
                2.0062;

        //create matrix for cross correlation Tc
        MatrixXd Tc = MatrixXd(n_x, n_z);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

        //calculate cross correlation matrix
        Tc.fill(0.0);
        for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

            //residual
            VectorXd z_diff = Zsig.col(i) - z_pred;
            //angle normalization
            while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
            while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

            // state difference
            VectorXd x_diff = Xsig_pred.col(i) - x;
            //angle normalization
            while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
            while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

            Tc = Tc + weights(i) * x_diff * z_diff.transpose();
        }

        //Kalman gain K;
        MatrixXd K = Tc * S.inverse();

        //residual
        VectorXd z_diff = z - z_pred;

        //angle normalization
        while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
        while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

        //update state mean and covariance matrix
        x = x + K * z_diff;
        P = P - K * S * K.transpose();
    }
