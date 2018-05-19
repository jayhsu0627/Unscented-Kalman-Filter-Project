#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/4.0;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  ///* State dimension
  n_x_ = 5;
  ///* Augmented state dimension
  n_aug_ = 7;
  augmented_columns = 2*n_aug_ + 1;
  Xsig_pred_ = MatrixXd::Zero(n_x_, augmented_columns);
  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;
  ///* Weights of sigma points
  weights_ = VectorXd(augmented_columns);
  weights_.fill(1.0/(2.0*(lambda_+n_aug_)));
  weights_[0] = lambda_/(lambda_+n_aug_);
  P_ = MatrixXd::Identity(n_x_, n_x_);
  Q_ = MatrixXd::Zero(2,2);
  Q_(0,0) = std_a_*std_a_;
  Q_(1,1) = std_yawdd_*std_yawdd_;
  R_lidar = MatrixXd::Zero(2, 2);
  R_lidar(0,0) = std_laspx_*std_laspx_;
  R_lidar(1,1) = std_laspy_*std_laspy_;
  R_radar = MatrixXd::Zero(3,3);
  R_radar(0,0) = std_radr_*std_radr_;
  R_radar(1,1) = std_radphi_*std_radphi_;
  R_radar(2,2) = std_radrd_*std_radrd_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_)
  {
    time_us_ = meas_package.timestamp_;
    x_.fill(0.0);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      float phi = meas_package.raw_measurements_[0];
      float rho = meas_package.raw_measurements_[1];
      float px = cos(rho)*phi;
      float py = sin(rho)*phi;
      x_(0) = px;
      x_(1) = py;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
    }

    is_initialized_ = true;
  }

  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;


  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
  else //LASER
  {
    UpdateLidar(meas_package);
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
  MatrixXd augmentedSigmaPoints = GenerateAugmentedSigmaPoints();
  PredictAugmentedSigmaPoints(augmentedSigmaPoints, delta_t);
  PredictMeanAndCovariance();
}

MatrixXd UKF::GenerateAugmentedSigmaPoints() {

  //create augmented mean vector
  VectorXd x_aug = VectorXd::Zero(n_aug_);
  x_aug.head(n_x_) = x_;

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2,2) = Q_;

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, augmented_columns);

  MatrixXd A = P_aug.llt().matrixL();
  Xsig_aug.col(0)  = x_aug;
  MatrixXd offset = sqrt(lambda_+n_aug_)*A;

  for (int i = 0; i < n_aug_; i++)
  {
      Xsig_aug.col(i+1) = x_aug + offset.col(i);
      Xsig_aug.col(i+1+n_aug_) = x_aug - offset.col(i);
  }

  return Xsig_aug;
}


void UKF::PredictAugmentedSigmaPoints(MatrixXd augmentedSigmaPoints, float dt)
{
  Xsig_pred_.fill(0.0);

  for(int i = 0; i < augmented_columns; i++)
  {
      Xsig_pred_.col(i) = PredictAugmentedSigmaPoint(augmentedSigmaPoints.col(i), dt);
  }
}

VectorXd UKF::PredictAugmentedSigmaPoint(VectorXd xsig_pred_col, float dt)
{
    double px = xsig_pred_col[0];
    double py = xsig_pred_col[1];
    double v = xsig_pred_col[2];
    double yaw = xsig_pred_col[3];
    double yaw_rate = xsig_pred_col[4];
    double acceleration_noise = xsig_pred_col[5];
    double yaw_acceleration_noise = xsig_pred_col[6];

    bool isYawRateZero = fabs(yaw_rate) < 0.01;

    double px_predicted = 0;
    double py_predicted = 0;

    if(isYawRateZero)
    {
        px_predicted = v*cos(yaw)*dt;
        py_predicted = v*sin(yaw)*dt;
    }
    else
    {
        px_predicted = (v/yaw_rate)*(sin(yaw+yaw_rate*dt)-sin(yaw));
        py_predicted = (v/yaw_rate)*(-cos(yaw+yaw_rate*dt)+cos(yaw));
    }

    double v_predicted = 0;
    double yaw_predicted = yaw_rate*dt;
    double yaw_rate_predicted = 0;

    double dt_squared = dt*dt;
    double px_noise = 0.5*dt_squared*cos(yaw)*acceleration_noise;
    double py_noise = 0.5*dt_squared*sin(yaw)*acceleration_noise;
    double v_noise = dt*acceleration_noise;
    double yaw_noise = 0.5*dt_squared*yaw_acceleration_noise;
    double yaw_rate_noise = dt*yaw_acceleration_noise;

    VectorXd x_predicted = VectorXd(5);
    x_predicted(0) = px + px_predicted + px_noise;
    x_predicted(1) = py + py_predicted + py_noise;
    x_predicted(2) = v + v_predicted + v_noise;
    x_predicted(3) = yaw + yaw_predicted + yaw_noise;
    x_predicted(4) = yaw_rate + yaw_rate_predicted + yaw_rate_noise;


    return x_predicted;
}

void UKF::PredictMeanAndCovariance(){

    //predict state mean
    x_.fill(0.0);
    for (int i = 0; i < augmented_columns; i++)
    {
        x_ = x_ + weights_[i]*Xsig_pred_.col(i);
    }

    //predict state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < augmented_columns; i++)
    {
      VectorXd aux = SubstractAndKeepAngleNormalizedIfNecessary(Xsig_pred_.col(i), x_, 3);
      P_ = P_ + weights_[i]*aux*aux.transpose();
    }

}
VectorXd UKF::SubstractAndKeepAngleNormalizedIfNecessary(VectorXd a, VectorXd b, int indexOfAngle){
  VectorXd aux = a - b;
  if (indexOfAngle != -1)
  {
    while (aux(indexOfAngle)> M_PI) aux(indexOfAngle)-=2.*M_PI;
    while (aux(indexOfAngle)<-M_PI) aux(indexOfAngle)+=2.*M_PI;
  }

  return aux;
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
  MatrixXd Zsig = PredictLidar();

  int n_z = 2;
  VectorXd z = VectorXd::Zero(n_z);
  float px = meas_package.raw_measurements_[0];
  float py = meas_package.raw_measurements_[1];
  z << px, py;

  nis_lidar = UpdateState(z, Zsig, R_lidar, -1);
  //cout << nis_lidar << endl;
}
MatrixXd UKF::PredictLidar()
{
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z, augmented_columns);
  for (int i = 0; i < augmented_columns; i++)
  {
    double px = Xsig_pred_.col(i)[0];
    double py = Xsig_pred_.col(i)[1];


    Zsig(0,i) = px;
    Zsig(1,i) = py;
  }

  return Zsig;
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
  MatrixXd Zsig = PredictRadar();

  int n_z = 3;
  VectorXd z = VectorXd::Zero(n_z);
  z(0) = meas_package.raw_measurements_[0];
  z(1) = meas_package.raw_measurements_[1];
  z(2) = meas_package.raw_measurements_[2];

  nis_radar = UpdateState(z, Zsig, R_radar, 1);
  //cout << nis_radar << endl;
}
  MatrixXd UKF::PredictRadar()
{
  int n_z = 3;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z, augmented_columns);

  //mean predicted measurement
  for (int i = 0; i < augmented_columns; i++)
  {
    double px = Xsig_pred_.col(i)[0];
    double py = Xsig_pred_.col(i)[1];
    double v = Xsig_pred_.col(i)[2];
    double yaw = Xsig_pred_.col(i)[3];

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    double rho = sqrt(px*px + py*py);
    double phi = atan2(py,px);

    double rho_rate = 0.0;
    if (fabs(rho) > 0.1)
      rho_rate = (px*cos(yaw)*v + py*sin(yaw)*v)/rho;


    Zsig(0,i) = rho;
    Zsig(1,i) = phi;
    Zsig(2,i) = rho_rate;
  }

  return Zsig;
}

float UKF::UpdateState(VectorXd z, MatrixXd Zsig, MatrixXd R, int indexToNormalizeInZ)
{
  int n_z = z.size();
  //mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);
  for(int i = 0; i < augmented_columns; i++)
  {
    z_pred = z_pred + weights_[i]*Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z, n_z);
  for(int i = 0; i < augmented_columns; i++)
  {
    //VectorXd aux = Zsig.col(i) - z_pred;
    VectorXd aux = SubstractAndKeepAngleNormalizedIfNecessary(Zsig.col(i), z_pred, indexToNormalizeInZ);
    S = S + weights_[i]*aux*aux.transpose();
  }

  S = S + R;

  // Update State
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

  for(int i = 0; i < augmented_columns; i++)
  {

    //VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd z_diff = SubstractAndKeepAngleNormalizedIfNecessary(Zsig.col(i), z_pred, indexToNormalizeInZ);

    VectorXd x_diff = SubstractAndKeepAngleNormalizedIfNecessary(Xsig_pred_.col(i), x_, 3);

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();


  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;

  x_ = x_ + K*z_diff;
  P_ = P_ - K*S*K.transpose();

  return z_diff.transpose()*S.inverse()*z_diff;
}
