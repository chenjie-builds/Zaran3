#include"LastSquare.h"
using namespace zaran;
DVector LastSquare::Solver()
{
  if (solver_type_ == 0)
    SolverNormalMatrix();
  else if (solver_type_ == 1)
    SolverSVD();
  else
    SolverColPivHouseholderQR();
  return x_;
}
void LastSquare::SolverSVD()
{
  x_ = A_.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_).transpose();
}

void LastSquare::SolverColPivHouseholderQR()
{
  x_ = A_.colPivHouseholderQr().solve(b_).transpose();
}

void LastSquare::SolverNormalMatrix()
{
  x_ = (A_.transpose() * A_).ldlt().solve(A_.transpose() * b_).transpose();
}
