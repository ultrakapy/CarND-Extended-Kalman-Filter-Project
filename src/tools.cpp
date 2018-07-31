#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::invalid_argument;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() > 0 && ground_truth.size() > 0 
      && estimations.size() == ground_truth.size()) {
    //accumulate squared residuals
    VectorXd diff;
    for(unsigned int i=0; i < estimations.size(); ++i){
      diff = estimations[i] - ground_truth[i];
      rmse = rmse.array() + diff.array()*diff.array();
    }
    
    //calculate the mean
    rmse = rmse / estimations.size();
    
    //calculate the squared root
    rmse = rmse.array().sqrt();
  }
  
  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float posSquaredSum = px*px + py*py;
  float sqrtPos = sqrt(posSquaredSum);
  float posPow32 = pow(posSquaredSum, 3/2);
  
  //check division by zero
  if ((px == 0 && py == 0) || fabs(posSquaredSum) < 0.0001) {
    std::cout << "HERE" << "\n";
    throw invalid_argument("CalculateJacobian() - Division by Zero");
  } else {
    //compute the Jacobian matrix
        
    Hj << px/sqrtPos, py/sqrtPos, 0, 0,
      -(py/posSquaredSum), px/posSquaredSum, 0, 0,
      (py*(vx*py-vy*px))/posPow32, (px*(vy*px-vx*py))/posPow32, px/sqrtPos, py/sqrtPos;
  }
  
  return Hj;
}
