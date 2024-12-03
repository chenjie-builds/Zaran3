#include "Solver.h"
#include"Log.h"
namespace zaran
{


Solver::Solver(int index, string name, SolverParam* para, GridBase* grid):
	index_(index),name_(name),m_para(para),m_grid(grid)
{
}

Solver::~Solver()
{
	delete[] m_para;
}

SolverParam* Solver::GetPara()
{
	Log::info("using Solver::GetPara()");
	return m_para;
}

void Solver::Init()
{
	m_para->Init();
}

}