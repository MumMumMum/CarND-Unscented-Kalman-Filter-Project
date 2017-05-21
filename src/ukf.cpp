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

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2; //0.2
  //[.09, .10, .40, .30].
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .7;

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = 5;

  n_aug_ =  7;
  n_aug_sigma = 2*n_aug_+1;

  lambda_ =  3 - n_aug_;

  NIS_radar_ = 7;

  NIS_laser_= 6;

  weights_ = VectorXd(n_aug_sigma);

  //InitWeights();
  weights_.fill(0.5 /(lambda_ + n_aug_));
  weights_(0) = lambda_/(lambda_ + n_aug_);
  Xsig_pred_ = MatrixXd(n_x_, n_aug_sigma);
  Xsig_aug_ = MatrixXd(n_aug_, n_aug_sigma);

  //add measurement noise covariance matrix
  R_lidar_ = MatrixXd(2,2);
  R_lidar_ <<    std_a_*std_a_, 0,
                0, std_a_*std_a_;

  //add measurement noise covariance matrix
  R_radar_ = MatrixXd(3,3);
  R_radar_ <<    std_radr_*std_radr_, 0, 0,
                 0, std_radphi_*std_radphi_, 0,
                 0, 0,std_radrd_*std_radrd_;

  is_initialized_ = false;
}
/*
// As per suggestion in review removed for loop to initalize weights
 //Thus init weights function got eliminated
void UKF:: InitWeights(void)
{
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<n_aug_sigma; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

}
*/


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
  //Initalization
  if(!is_initialized_){

      x_ << 1, 1, .1, .1,.1;
      if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

        /**
           Convert radar from polar to cartesian coordinates and initialize state.
           */
        float radialDistance = meas_package.raw_measurements_[0];
        float anglePhi = meas_package.raw_measurements_[1];
        float radialValocity = meas_package.raw_measurements_[2];
        x_(0) = radialDistance * cos(anglePhi);
        x_(1) = radialDistance * sin(anglePhi);
        x_(2) = radialValocity/100 ;// this is objects radial velocity and not mesaured
        x_(3) = 0 ;
        x_(4) = 0;


      }else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
       x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0,0;
      }
    P_ <<   0.5, 0, 0, 0, 0, //px
            0, 0.5, 0, 0, 0,   //py
            0, 0, 0.1, 0, 0, //v
            0, 0, 0, 0.01, 0,     //turning angle psi;
            0, 0, 0, 0, 0.01;// rate of turning angle psi_dt
      previous_timestamp_ = meas_package.timestamp_;
      is_initialized_ = true;
      return;
    }
    float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;    //dt - expressed in seconds
    previous_timestamp_ = meas_package.timestamp_;
    //Predict the state vector and covariance matrix.
    //cout<<"sensor type "<<meas_package.sensor_type_<<endl;

    Prediction( dt);
  /*****************************************************************************
 *  Update
 ****************************************************************************/

  /**
    * Use the sensor type to perform the the measurement prediction
   */

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {


    UpdateRadar(meas_package);

  } else {

    UpdateLidar(meas_package);
  }
  }


void UKF::Augmentation_X_(void){
//create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  Xsig_aug_ = MatrixXd(n_aug_, n_aug_sigma);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //create augmented mean state
  x_aug.head(5) = x_;
  //Longitudnal noise for acc
  x_aug(5) = 0;
  //Longitudnal noise for psi
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug_.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
}

// to find sigma points we try to find out a matrix whose covariance and mean is  same
//as our give vector x and covariance p
void UKF::SigmaPointPrediction(double delta_t ){


  for (int i = 0; i< n_aug_sigma; i++)
  {
    //extract values for better readability
    double p_x           = Xsig_aug_(0,i);
    double p_y           = Xsig_aug_(1,i);
    double v             = Xsig_aug_(2,i);
    double yaw           = Xsig_aug_(3,i);
    double yawd          = Xsig_aug_(4,i);
    double nu_a          = Xsig_aug_(5,i);
    double nu_yawdd      = Xsig_aug_(6,i);

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
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
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
  Augmentation_X_();
  SigmaPointPrediction(delta_t);
  /*x_.fill(0.0);
  for (int i = 0; i < n_aug_sigma; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }*/
  // This suggestion of array mathematics was given in review.
  x_ = Xsig_pred_ * weights_;


  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_aug_sigma; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization

    /*while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;*/

    tools_.Normalize_angle(x_diff(3));
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
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
  int n_z = 2;
  /*MatrixXd Zsig = MatrixXd(n_z, n_aug_sigma);
  //transform sigma points into measurement space
  for (int i = 0; i < n_aug_sigma; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    Zsig(0,i) = p_x;  //px
    Zsig(1,i) = p_y; //py

  }*/
  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z,  n_aug_sigma);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < n_aug_sigma; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < n_aug_sigma; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  S = S + R_lidar_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_aug_sigma; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization


    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;


    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z = VectorXd(n_z);
  //z = meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
  z = meas_package.raw_measurements_;
  //z = meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
  VectorXd z_diff = z - z_pred;

  NIS_laser_  =  z_diff.transpose()*S.inverse()*z_diff;
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
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

  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, n_aug_sigma);
  //transform sigma points into measurement space
  for (int i = 0; i < n_aug_sigma; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
                          //r
    double eps = 0.001;
    double c1 = std::max(eps,sqrt(p_x*p_x + p_y*p_y) );
    Zsig(0,i) = c1;
    if (p_y == 0 & p_x == 0)
      Zsig(1,i) = 0.001;
    else
      Zsig(1,i) = atan2(p_y,p_x);
    //Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / c1;   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < n_aug_sigma; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < n_aug_sigma; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
  /*  while (z_diff(1)>= M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<= -M_PI) z_diff(1)+=2.*M_PI;*/

    tools_.Normalize_angle(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }


  S = S + R_radar_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_aug_sigma; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
   /* while (z_diff(1)>= M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<= -M_PI) z_diff(1)+=2.*M_PI;*/
    tools_.Normalize_angle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    /*while (x_diff(3)>= M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<=-M_PI) x_diff(3)+=2.*M_PI;*/
    tools_.Normalize_angle(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z = VectorXd(n_z);
  //z = meas_package.raw_measurements_[0], meas_package.raw_measurements_[1],meas_package.raw_measurements_[2];
  z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  /*//angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;*/

  tools_.Normalize_angle(z_diff(1));

  NIS_radar_  =  z_diff.transpose()*S.inverse()*z_diff;
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

}

/*

void UKF::UpdateLidar(MeasurementPackage meas_package) {}
void UKF::UpdateRadar(MeasurementPackage meas_package) {}*/
