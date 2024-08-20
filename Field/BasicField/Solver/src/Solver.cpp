#include "Solver.h"
#include"Log.h"
using namespace zaran;

zaran::Solver::Solver(int index, string name, SolverPara* para, GridBase* grid):
	index_(index),name_(name),m_para(para),m_grid(grid)
{
}

Solver::~Solver()
{
	delete[] m_para;

}

SolverPara* Solver::GetPara()
{
	Log::info("using Solver::GetPara()");
	return m_para;
}

void zaran::Solver::Init()
{
	m_para->Init();
}
