#include <stdlib.h>
#include <Rcpp.h>
using namespace Rcpp;

//' @title rwMetropolis using Rcpp
//' @description rwMetropolis using Rcpp
//' @param n the number of samples
//' @param sigma sigma of normal distribution
//' @param x0 initial value
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' rwMetropolis(1,50,10)
//' }
//' @export
// [[Rcpp::export]]
List rwMetropolis(float sigma,float x0,int n) {
  List out(n);
  out[0]=0;
  out[1]=x0;
  int k=0;
  float y=0;
  for(int i = 2; i < n; i++) {
   y=rnorm(1, out[i-1], sigma)[0];
   if (runif(1)[0] <= (exp(fabs(out[i-1])-fabs(y))))
    out[i] = y; 
	else {
    out[i] = out[i-1];
    k++;
	}    
}
  out[0]=k;
  return out;
}
