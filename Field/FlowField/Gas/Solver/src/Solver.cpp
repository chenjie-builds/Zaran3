#include "solver.h"
using namespace zaran;
Solver::Solver()
{
	index_ = -1;
	name_ = "notSet";
	gridList_ = nullptr;
	gridIndex_ = -1;
	para_ = nullptr;
}

Solver::~Solver()
{

}

SolverParaPtr& Solver::GetPara()
{
	return para_;
}
