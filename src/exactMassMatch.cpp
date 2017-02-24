#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List exactMassMatchCpp(List Res, NumericVector unknowns, NumericVector metMasses, 
                       CharacterVector metMassesNames, double ppm, double mill){
  int n = unknowns.size();
  for(int i = 0; i < n; ++i){
    Res[i] = metMassesNames[(metMasses < unknowns[i] + (unknowns[i] / mill) * ppm) & (metMasses > unknowns[i] - (unknowns[i] / mill) * ppm)];
  }
  return Res;
}