#include "solver.h"
#include"Log.h"
using namespace zaran;

zaran::Solver::Solver(int index, string name, SolverPara* para, GridBase* grid):
	index_(index),name_(name),para_(para),m_grid(grid)
{
}

Solver::~Solver()
{
	delete[] para_;

}

SolverPara* Solver::GetPara()
{
	Log::info("using Solver::GetPara()");
	return para_;
}

void zaran::Solver::Init()
{
	para_->Init();
}
