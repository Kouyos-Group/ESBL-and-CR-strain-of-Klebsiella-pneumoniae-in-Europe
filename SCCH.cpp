#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector SCCH(double SC, double CC, double SH, double CH,
                   double beta_community,
                   double beta_hospital,
                   double clearance_rate,
                   double hospitalization_rate,
                   double discharge_rate) {
  int N = 1000000;
  NumericMatrix out(4,N);
  NumericVector out1(4);
  out(0,0) = SC;
  out(1,0) = CC;
  out(2,0) = SH;
  out(3,0) = CH;
  double t_step = 0.01;
  for (int i = 1; i <= (N-1); i += 1) {
    out(0,i)= out(0,i-1) + t_step*(- out(0,i-1) * hospitalization_rate + out(2,i-1) * discharge_rate - beta_community * out(0,i-1) * out(1,i-1)/(out(0,i-1)+out(1,i-1)) + out(1,i-1)  * clearance_rate);
    
    out(1,i)= out(1,i-1) + t_step*(- out(1,i-1) * hospitalization_rate + out(3,i-1) * discharge_rate + beta_community * out(0,i-1) * out(1,i-1)/(out(0,i-1)+out(1,i-1)) - out(1,i-1)  * clearance_rate);
    
    out(2,i)= out(2,i-1) + t_step*(+ out(0,i-1) * hospitalization_rate - out(2,i-1) * discharge_rate - beta_hospital * out(2,i-1) * out(3,i-1)/(out(2,i-1)+out(3,i-1)) + out(3,i-1)  * clearance_rate);
    
    out(3,i)= out(3,i-1) + t_step*(+ out(1,i-1) * hospitalization_rate - out(3,i-1) * discharge_rate + beta_hospital * out(2,i-1) * out(3,i-1)/(out(2,i-1)+out(3,i-1)) - out(3,i-1)  * clearance_rate);
    
  }
  out1(0) = out(0,(N-1));
  out1(1) = out(1,(N-1));
  out1(2) = out(2,(N-1));
  out1(3) = out(3,(N-1));
  return(out1);
}



