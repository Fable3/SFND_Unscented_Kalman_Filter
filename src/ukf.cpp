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
  //P_ = MatrixXd(5, 5);
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  // Bicycles can slowly accelerating, but they can break efficiently
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  // it would take 1 second from going straight (yawd=0) to turn into a circle which takes 16 seconds to complete
  std_yawdd_ = 2*M_PI/16;
  
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
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  // create vector for weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  // set weights
  int i;
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (i = 1; i < 2 * n_aug_ + 1; i++) weights_(i) = 1 / (lambda_ + n_aug_) / 2;

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
	if (!is_initialized_)
	{
		//x: [pos1 pos2 vel_abs yaw_angle yaw_rate]
		x_ = VectorXd::Zero(5);
		if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER)
		{
			x_(0) = meas_package.raw_measurements_(0);
			x_(1) = meas_package.raw_measurements_(1);
		}
		else
		{
			double p = meas_package.raw_measurements_(0);
			double fi = meas_package.raw_measurements_(1);
			x_(0) = cos(fi)*p;
			x_(1) = sin(fi)*p;
		}
		time_us_ = meas_package.timestamp_;
		is_initialized_ = true;
	}
	else
	{
		if (!use_laser_ &&  meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER ||
			!use_radar_ &&  meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR)
		{
			return;
		}
		double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
		time_us_ = meas_package.timestamp_;
		Prediction(dt);
		if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER)
		{
			UpdateLidar(meas_package);
		}
		else
		{
			UpdateRadar(meas_package);
		}

	}
}


void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
	// create augmented mean vector
	VectorXd x_aug = VectorXd::Zero(n_aug_);

	// create augmented state covariance
	MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);

	// create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	// create augmented mean state
	x_aug.head(n_x_) = x_;

	// create augmented covariance matrix
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug(5, 5) = std_a_ * std_a_;
	P_aug(6, 6) = std_yawdd_ * std_yawdd_;

	// create square root matrix
	MatrixXd A = P_aug.llt().matrixL();

	// create augmented sigma points
	// set first column of sigma point matrix
	Xsig_aug.col(0) = x_aug;

	// set remaining sigma points
	for (int i = 0; i < n_aug_; ++i) {
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
	}

	// create matrix with predicted sigma points as columns
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

	// predict sigma points
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		double px = Xsig_aug(0, i);
		double py = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);
		// avoid division by zero
		if (fabs(yawd) > 0.0001)
		{
			Xsig_pred_(0, i) = px + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw)) + 0.5*delta_t*delta_t*cos(yaw)*nu_a;
			Xsig_pred_(1, i) = py + v / yawd * (-cos(yaw + yawd * delta_t) + cos(yaw)) + 0.5*delta_t*delta_t*sin(yaw)*nu_a;
			Xsig_pred_(2, i) = v + 0 + delta_t * nu_a;
			Xsig_pred_(3, i) = yaw + yawd * delta_t + 0.5*delta_t*delta_t*nu_yawdd;
			Xsig_pred_(4, i) = yawd + 0 + delta_t * nu_yawdd;
		}
		else
		{
			Xsig_pred_(0, i) = px + v * cos(yaw)*delta_t + 0.5*delta_t*delta_t*cos(yaw)*nu_a;
			Xsig_pred_(1, i) = py + v * sin(yaw)*delta_t + 0.5*delta_t*delta_t*sin(yaw)*nu_a;
			Xsig_pred_(2, i) = v + 0 + delta_t * nu_a;
			Xsig_pred_(3, i) = yaw + yawd * delta_t + 0.5*delta_t*delta_t*nu_yawdd;
			Xsig_pred_(4, i) = yawd + 0 + delta_t * nu_yawdd;
		}
	}

	// predict state mean
	x_ = VectorXd::Zero(n_x_);
	int i;
	for (i = 0; i < 2 * n_aug_ + 1; i++)
		x_ += weights_(i)*Xsig_pred_.col(i);

	// predict state covariance matrix
	P_ = MatrixXd::Zero(n_x_, n_x_);
	for (i = 0; i < 2 * n_aug_ + 1; i++)
	{
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		x_diff(3) = fmod(x_diff(3) + 5*M_PI, 2 * M_PI) - M_PI;
		P_ += weights_(i)*(x_diff*x_diff.transpose());
	}
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
	int n_z = 2;
	// create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	// mean predicted measurement
	VectorXd z_pred = VectorXd::Zero(n_z);

	// measurement covariance matrix S
	MatrixXd S = MatrixXd::Zero(n_z, n_z);

	// transform sigma points into measurement space
	int i;
	for (i = 0; i < 2 * n_aug_ + 1; ++i)
	{
		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);
		Zsig(0, i) = px;
		Zsig(1, i) = py;
	}

	// calculate mean predicted measurement
	for (i = 0; i < 2 * n_aug_ + 1; ++i)
	{
		z_pred += weights_(i)*Zsig.col(i);
	}

	// calculate innovation covariance matrix S
	for (i = 0; i < 2 * n_aug_ + 1; ++i)
	{
		VectorXd z_diff = Zsig.col(i) - z_pred;
		S += weights_(i)*(z_diff*z_diff.transpose());
	}
	S(0, 0) += std_laspx_ * std_laspx_;
	S(1, 1) += std_laspy_ * std_laspy_;

	// create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
	// calculate cross correlation matrix
	for (i = 0; i < 2 * n_aug_ + 1; ++i)
	{
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		x_diff(3) = fmod(x_diff(3) + M_PI, 2 * M_PI) - M_PI;
		VectorXd z_diff = Zsig.col(i) - z_pred;
		Tc += weights_(i)*x_diff*z_diff.transpose();
	}

	// calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	// update state mean and covariance matrix
	VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
	x_ = x_ + K * z_diff;
	P_ = P_ - K * S*K.transpose();
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
   
	// set measurement dimension, radar can measure r, phi, and r_dot
	int n_z = 3;

	// create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	// mean predicted measurement
	VectorXd z_pred = VectorXd::Zero(n_z);

	// measurement covariance matrix S
	MatrixXd S = MatrixXd::Zero(n_z, n_z);

	// transform sigma points into measurement space
	int i;
	for (i = 0; i < 2 * n_aug_ + 1; ++i)
	{
		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);
		double yawd = Xsig_pred_(4, i);
		double dist = sqrt(px*px + py * py);
		Zsig(0, i) = dist;
		Zsig(1, i) = atan2(py, px);
		if (dist < 0.0001)
		{
			Zsig(2, i) = 0;
		}
		else
		{
			Zsig(2, i) = (px*cos(yaw)*v + py * sin(yaw)*v) / dist;
		}
	}

	// calculate mean predicted measurement
	for (i = 0; i < 2 * n_aug_ + 1; ++i)
	{
		z_pred += weights_(i)*Zsig.col(i);
	}

	// calculate innovation covariance matrix S
	for (i = 0; i < 2 * n_aug_ + 1; ++i)
	{
		VectorXd z_diff = Zsig.col(i) - z_pred;
		z_diff(1) = fmod(z_diff(1) + 5 * M_PI, 2 * M_PI) - M_PI;
		S += weights_(i)*(z_diff*z_diff.transpose());
	}
	S(0, 0) += std_radr_ * std_radr_;
	S(1, 1) += std_radphi_ * std_radphi_;
	S(2, 2) += std_radrd_ * std_radrd_;

	// create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
	// calculate cross correlation matrix
	for (i = 0; i < 2 * n_aug_ + 1; ++i)
	{
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		x_diff(3) = fmod(x_diff(3) + M_PI, 2 * M_PI) - M_PI;
		VectorXd z_diff = Zsig.col(i) - z_pred;
		z_diff(1) = fmod(z_diff(1) + 5*M_PI, 2 * M_PI) - M_PI;
		Tc += weights_(i)*x_diff*z_diff.transpose();
	}

	// calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	// update state mean and covariance matrix
	VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
	z_diff(1) = fmod(z_diff(1) + 5*M_PI, 2 * M_PI) - M_PI;
	x_ = x_ + K * z_diff;
	P_ = P_ - K * S*K.transpose();
}
