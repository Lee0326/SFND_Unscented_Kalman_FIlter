#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  is_initialized_ = false;

  use_laser_ = false;

  use_radar_ = false;

  Xsig_pred_ = MatrixXd(5,11);

  n_x_ = 5;

  n_aug_ = 7;
  
  lambda_ = 3 - n_x_;
       
  P_<< 0.15, 0, 0, 0, 0,
      0, 0.15, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  // Generating Sigma Points
  if (!is_initialized_)
  {
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      x_[0] = meas_package.raw_measurements_[0];
      x_[1] = meas_package.raw_measurements_[1];
      x_[2] = 0;
      x_[3] = 0;
      x_[4] = 0; 

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rhodot = meas_package.raw_measurements_[2];
      x_[0] = rho*cos(phi);
      x_[1] = rho*sin(phi);
      x_[2] = rhodot;
    }
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return; 
  }
  // Predict
  double delta_t = meas_package.timestamp_-time_us_;
  Prediction(delta_t);
  //Update
  if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  // Generating the Sigma Points
  MatrixXd Xsig = MatrixXd(n_x_, 2*n_x_+1);
  MatrixXd A = P_.llt().matrixL();

  for (int i=0; i<2*n_x_+1; i++)
  {
    Xsig.col(i) = x_;
  }
  Xsig.block<5,5>(0,1) += sqrt(lambda_+n_x_)*A;
  Xsig.block<5,5>(0,n_x_+1) -= sqrt(lambda_+n_x_)*A;
  // Augmentation
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_,2*n_aug_+1);
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
  
  x_aug.head(5) = x_;
  x_aug[5] = 0;
  x_aug[6] = 0;

  P_aug.fill(0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //Create Augmented Sigma Points
  MatrixXd L = P_aug.llt().matrixL();
  for (int i=0; i<2*n_aug_+1; i++)
  {
    Xsig_aug.col(i) = x_aug;
  }
  Xsig_aug.block<7,7>(0,1) += sqrt(lambda_+n_aug_)*L;
  Xsig_aug.block<7,7>(0,n_aug_+1) += sqrt(lambda_+n_aug_)*L;
  // Sigma Point Prediction
  

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}