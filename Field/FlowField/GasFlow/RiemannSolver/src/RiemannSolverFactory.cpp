#include "RiemannSolverFactory.h"
#include "Vanleer.h"
#include "Ausmpw.h"
#include "Roe.h"
#include "HLLC.h"
#include "StegerWarming.h"
#include "log.h"
using namespace zaran;


RiemannSolver* RiemannSolverBuilder::Create(RiemannSolverType type)
{
	if (type == RiemannSolverType::VanLeer)
	{
		return new Vanleer();
	}
	else if (type == RiemannSolverType::Ausmpw)
	{
		return new Ausmpw();
	}
	else if (type == RiemannSolverType::Roe)
	{
		return new Roe();
	}
	else if (type == RiemannSolverType::HLLC)
	{
		return new HLLC();
	}
	else if (type == RiemannSolverType::StegerWarming)
	{
		return new StegerWarming();
	}
	else
	{
		Log::warn("RiemannSolverBuilder::Create: unknown RiemannSolverType");
		Log::warn("RiemannSolverBuilder::Create: Using default RiemannSolverType: VanLeer");
		return new Vanleer();
	}
}
