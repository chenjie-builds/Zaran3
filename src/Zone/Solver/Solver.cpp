#include "Solver.h"
#include"Log.h"
namespace zaran
{


Solver::Solver(Id index, string name, std::shared_ptr<SolverParam>para, std::shared_ptr<GridBase>grid) :
	index_(index),name_(name),m_para(para),m_grid(grid)
{
}

Solver::~Solver()
{
}


void Solver::Init()
{
	m_para->Init();
}

}