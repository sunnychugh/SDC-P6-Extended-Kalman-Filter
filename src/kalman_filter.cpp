#include "kalman_filter.h"
#include "tools.h"
#include <iostream>
#include <math.h>

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict()
{
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z)
{
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  /*
  For lidar measurements, the error equation is y = z - H * x'. 
  For radar measurements, the functions that map the x vector [px, py, vx, vy] to polar coordinates are 
  non-linear. Instead of using H to calculate y = z - H * x', for radar measurements 
  you'll have to use the equations that map from cartesian to polar coordinates: y = z - h(x').
  */
  // VectorXd z_pred = H_ * x_; // avoid this way , hence update later

  double px = x_[0];
  double py = x_[1];
  double vx = x_[2];
  double vy = x_[3];

  // to avoid nan values if starting readings of radar data is zero
  double epsilon = 0.0001;
  if ((fabs(px) < epsilon && fabs(py) < epsilon) || (fabs((px * px) + (py * py)) < epsilon))
  {
    px = epsilon;
    py = epsilon;
  }

  double rho = sqrt(px * px + py * py);
  double phi = atan2(py, px);
  double rho_dot = (px * vx + py * vy) / rho;

  VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot;
  // measurement - prediction in polar coordinates
  VectorXd y = z - z_pred;

  //cout << "z_pred[1]: " << z_pred[1] << endl;
  //cout << "y[1]: " << y[1] << endl;

  // The Kalman filter is expecting small angle values between the range -pi and pi.
  // HINT: when working in radians, you can add 2π or subtract 2π until the angle is within the desired range.
  while (y[1] < -M_PI)
    y[1] += 2 * M_PI;
  while (y[1] > M_PI)
    y[1] -= 2 * M_PI;

  //cout << "y[1]: " << y[1] << endl;

  Tools t;
  MatrixXd Hj_ = t.CalculateJacobian(x_);

  MatrixXd Hjt = Hj_.transpose();
  MatrixXd S = Hj_ * P_ * Hjt + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHjt = P_ * Hjt;
  MatrixXd K = PHjt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj_) * P_;
}
