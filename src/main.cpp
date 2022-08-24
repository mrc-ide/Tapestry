
#include "main.h"

using namespace Rcpp;

//------------------------------------------------
// square a vector of values
List square_cpp(List args) {
  
  // print message to console
  Rcout << "running C++ square_cpp function\n";
  
  // get inputs from Rcpp format to base C++ format
  std::vector<double> x = as<std::vector<double>>(args("x"));
  
  // square values
  for (int i = 0; i < int(x.size()); ++i) {
    x[i] *= x[i];
  }
  
  // return as Rcpp list
  List ret = List::create(Named("x_squared") = x);
  return ret;
}
