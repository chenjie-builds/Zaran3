#include "solver.h"
#include"Log.h"
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

SolverPara* Solver::GetPara()
{
	Log::info("using Solver::GetPara()");
	return para_;
}
