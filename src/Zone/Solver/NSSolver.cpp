#include "NSSolver.h"
#include "File.h"
#include <fstream>
#include "Log.h"
#include "PerfectGas.h"
#include "RiemannSolverFactory.h"
namespace zaran
{

	NSSolver::NSSolver(index_type index, string name, shared_ptr<FlowSolverParam> para, shared_ptr < GridBase> grid, shared_ptr < DataManagerNS>data_manager) :
		FlowFieldSolver(index, name, para, grid, data_manager)
	{
		m_data_manager = data_manager;
	}

	NSSolver::~NSSolver()
	{
	}

	void NSSolver::Init()
	{
		InitSolver();
		InitField();
		CalcCoordTransCoef();
	}
	void NSSolver::InitField()
	{
		auto para = GetPara();
		auto init_type = para->GetInitFieldType();
		if (init_type == InitFieldType::FarFlowNoVelocity)
		{
			InitFieldFarFieldZeroVel();
		}
		else if (init_type == InitFieldType::FarFlow)
		{
			InitFieldFarfield();
		}
		else if (init_type == InitFieldType::Backup)
		{
			InitFieldBackup();
		}
		else if (init_type == InitFieldType::Explosion)
		{
			InitFieldExplosion();
		}
		else
		{
			Log::warn("Initialize Failed!");
			Log::warn("Wrong Flow field init parameter:{}", GlobalData::GetString("initFieldType"));
			exit(0);
		}
		Prim2Cons();
		Log::info("Flow Field Initialize Finished!");
	}
	DataManagerNS* NSSolver::GetDataManager()
	{
		return m_data_manager.get();
	}
	void NSSolver::InitSolver()
	{
		Log::info("Initialize NS Solver!");
		FlowFieldSolver::InitSolver();
		RiemannSolverBuilder riemann_solver_builder;
		m_riemann_solver = riemann_solver_builder.CreateUnique(GetPara()->GetRiemannSolverType());
		auto& ref_value = GetPara()->GetDimensionless();
		m_gas = make_shared<PerfectGas>(ref_value.GetRefMw(), ref_value.GetRefGamma(), ref_value);
		Log::info("NS Solver Initialize Finished!");
	}

	void NSSolver::Preprocess()
	{
		auto para = GetPara();
		para->SetCurrentStep(para->GetCurrentStep() + 1);
		BoundaryCondition();
	}

	void NSSolver::Postprocess()
	{
		UpdateField();
		CheckPrimtive();
		FixPrimtive();
		CalcForce();
	}

	void NSSolver::Solve()
	{
		TimeAdvance();
	}

	void NSSolver::TimeAdvance()
	{
		RungeKutta();
	}

	void NSSolver::CalcResidual()
	{
		ZeroResidual();
		CalcConvectionResidual();
		CalcViscousResidual();
		CalcSourceResidual();
	}

	void NSSolver::UpdateField()
	{
		Cons2Prim();
	}

} // namespace zaran