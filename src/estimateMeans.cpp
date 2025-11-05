#include<Rcpp.h>

//' Print the version number
//' 
//' Print the version number and date.
//' @export
// [[Rcpp::export]]
void version(){
  Rcpp::Rcout << "TargetedShrinkage version 1.0, 3rd April 2023." << std::endl;
}