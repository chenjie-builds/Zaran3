#include "SolverVec.h"
#include "log.h"
using namespace zaran;
SolverVec::SolverVec()
{

}

SolverVec::~SolverVec()
{

}

void SolverVec::AddSolver(Ptr<Solver>& solver)
{
	solverVec_.emplace_back(solver);
}

void SolverVec::CheckSolver()
{
	for (int iSolver = 0; iSolver < solverVec_.size(); ++iSolver)
	{
		if (solverVec_[iSolver]->GetIndex() != iSolver)
		{
			ZaranLog::warn("the index of solver are not match to solver set!");
			ZaranLog::warn("solver index:{}, solver set index:{}", solverVec_[iSolver]->GetIndex(), iSolver);
		}
	}
}

Ptr<Solver>& SolverVec::GetSolverPtr(const int index)
{
	return solverVec_[index];
}
