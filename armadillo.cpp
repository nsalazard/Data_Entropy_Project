#include <armadillo>
#include <iostream>

int main() {
  arma::mat A(14, 15, arma::fill::randu);
  arma::mat B(14, 15, arma::fill::randu);

  std::cout << arma::det(A * B.t()) << std::endl;

  return 0;
}
