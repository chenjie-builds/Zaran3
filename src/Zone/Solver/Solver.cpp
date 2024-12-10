#include "Solver.h"
#include"Log.h"
namespace zaran
{


Solver::Solver(index_type index, string name, shared_ptr<SolverParam>para, shared_ptr<GridBase>grid) :
	m_id(index),m_name(name),m_para(para),m_grid(grid)
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