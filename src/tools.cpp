#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
 VectorXd rmse(4);
 rmse << 0,0,0,0;

 // check the validity of inputs
 if (estimations.size() == 0 || estimations.size() != ground_truth.size()){
   cout << "Not Validation Error!!!";
   return rmse; 
}
 
 // accumulate squared residuals
 for (int i=0; i < estimations.size(); ++i){
   VectorXd residual = estimations[i]-ground_truth[i];
   residual = residual.array() * residual.array();
   rmse += residual;  
 }

 rmse /= estimations.size();
 rmse = rmse.array().sqrt();

 return rmse;
}
