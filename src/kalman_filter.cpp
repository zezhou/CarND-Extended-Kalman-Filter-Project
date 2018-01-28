#include "kalman_filter.h"
#include <assert.h>
#include <math.h>
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
	P_ = F_ *  P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - H_ * x_;
  MeasurementUpdate(y, H_, R_);

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  
  VectorXd y = z - Convert(x_, z);
  MatrixXd Hj = CalJacobian(x_);
  MeasurementUpdate(y, Hj, R_);
}


VectorXd KalmanFilter::Convert(const VectorXd& cartesian, const VectorXd& polar_base) {
    float px = cartesian[0];
    float py = cartesian[1];
    float vx = cartesian[2];
    float vy = cartesian[3];

    float rho = sqrt(px * px + py * py);

    VectorXd polar(3);
    polar << rho,
             Normalize(atan2(py, px), polar_base[1]),
             (px * vx + py * vy) / rho;

    return polar;
}

void KalmanFilter::MeasurementUpdate(const VectorXd &y, const MatrixXd &H, const MatrixXd &R) {
    MatrixXd Ht = H.transpose();
    MatrixXd S = H * P_ * Ht + R;
    MatrixXd K = P_ * Ht * S.inverse();
    // New estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H) * P_;
}

MatrixXd KalmanFilter::CalJacobian(const VectorXd& x_state) {
    MatrixXd Hj(3,4);
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    float c1 = px * px + py * py;
    float c2 = sqrt(c1);
    float c3 = (c1 * c2);
    Hj << (px / c2), (py / c2), 0, 0,
          -(py / c1), (px / c1), 0, 0,
          py * (vx * py - vy * px) / c3, px * (vy * px - vx * py) / c3, px / c2, py / c2;

    return Hj;
}

float KalmanFilter::Normalize(const float& angle, const float& base_angle) {
    float norm_angle = angle;
    while (norm_angle - base_angle > M_PI) {
        norm_angle -= 2 * M_PI;
    }
    while (base_angle - norm_angle > M_PI) {
        norm_angle += 2 * M_PI;
    }
    return norm_angle;
}