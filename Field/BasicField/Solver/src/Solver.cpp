#include "solver.h"
#include"Log.h"
using namespace zaran;
Solver::Solver()
{
	index_ = -1;
	name_ = "notSet";
	m_grid = nullptr;
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
