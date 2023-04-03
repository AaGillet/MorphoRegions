#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

int prod(IntegerVector x) {
  int p = x[0];

  if (x.size() > 1) {
    for (int i = 1; i < x.size(); i++) {
      p *= x[i];
    }
  }

  return p;
}

// [[Rcpp::export]]
IntegerMatrix combosC(IntegerVector y, int m, int nbp) {
  int n = y.size();
  int p = n - m - m * nbp + 1;
  int ncombos = prod(p - 1 + seq_len(nbp))/prod(seq_len(nbp));

  IntegerVector mins = m * seq_len(nbp);
  IntegerVector maxes = rep(n, nbp) - rev(mins);
  IntegerVector counters = clone(mins);

  IntegerMatrix combos(ncombos, nbp);

  LogicalVector reset = rep(false, nbp);

  int i;
  for (int k = 0; k < ncombos; k++) {
    combos(k, _) = counters;
    counters[nbp - 1]++;
    for (i = nbp - 1; i > 0; i--) {
      if (counters[i] > maxes[i]) {
        reset[i] = true;
        counters[i - 1]++;
      }
    }
    if (any(reset).is_true()) {
      for (i = 1; i < nbp; i++) {
        if (reset[i]) {
          counters[i] = counters[i - 1] + m;
          reset[i] = false;
        }
      }
    }
  }

  return combos;
}
