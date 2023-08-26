#include "RiemannSolverFactory.h"
#include "Vanleer.h"
#include "Ausmpw.h"
#include "Roe.h"
#include "HLLC.h"
#include "StegerWarming.h"
#include "log.h"
using namespace zaran;

void RiemannSolverFactory::Create(std::shared_ptr<RiemannSolver>& riemannSolver, string& riemannSolverName)
{
	if (riemannSolverName == "VanLeer")
		CreateVanLeer(riemannSolver);
	else if (riemannSolverName == "Roe")
		CreateRoe(riemannSolver);
	else if (riemannSolverName == "Ausmpw")
		CreateAusmpw(riemannSolver);
	else if (riemannSolverName == "HLLC")
		CreateHLLC(riemannSolver);
	else if (riemannSolverName == "StegerWarming")
		CreateStegerWarming(riemannSolver);
	else
	{
		ZaranLog::warn("unsupportted riemann solver name: {}", riemannSolverName);
	}

}

void RiemannSolverFactory::CreateVanLeer(std::shared_ptr<RiemannSolver>& riemannSolver)
{
	riemannSolver = std::make_shared<Vanleer>();
}

void RiemannSolverFactory::CreateHLLC(std::shared_ptr<RiemannSolver>& riemannSolver)
{
	riemannSolver = std::make_shared<HLLC>();
}

void RiemannSolverFactory::CreateRoe(std::shared_ptr<RiemannSolver>& riemannSolver)
{
	riemannSolver = std::make_shared<Roe>();
}

void RiemannSolverFactory::CreateStegerWarming(std::shared_ptr<RiemannSolver>& riemannSolver)
{
	riemannSolver = std::make_shared<StegerWarming>();
}

void RiemannSolverFactory::CreateAusmpw(std::shared_ptr<RiemannSolver>& riemannSolver)
{
	riemannSolver = std::make_shared<Ausmpw>();
}
