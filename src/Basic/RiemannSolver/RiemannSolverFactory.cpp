#include "RiemannSolverFactory.h"
#include "Vanleer.h"
#include "Ausmpw.h"
#include "Roe.h"
#include "HLLC.h"
#include "StegerWarming.h"
#include "Log.h"
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

std::unique_ptr<zaran::RiemannSolver> RiemannSolverBuilder::CreateUnique(RiemannSolverType type)
{
	if (type == RiemannSolverType::VanLeer)
	{
		return std::make_unique<Vanleer>();
	}
	else if (type == RiemannSolverType::Ausmpw)
	{
		return std::make_unique<Ausmpw>();
	}
	else if (type == RiemannSolverType::Roe)
	{
		return std::make_unique<Roe>();
	}
	else if (type == RiemannSolverType::HLLC)
	{
		return std::make_unique<HLLC>();
	}
	else if (type == RiemannSolverType::StegerWarming)
	{
		return std::make_unique<StegerWarming>();
	}
	else
	{
		Log::warn("RiemannSolverBuilder::Create: unknown RiemannSolverType");
		Log::warn("RiemannSolverBuilder::Create: Using default RiemannSolverType: VanLeer");
		return std::make_unique<Vanleer>();
	}

}
