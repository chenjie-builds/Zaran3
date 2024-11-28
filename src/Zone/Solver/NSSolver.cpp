#include "NSSolver.h"
#include "File.h"
#include <fstream>
#include "Log.h"
#include "PerfectGas.h"
#include "RiemannSolverFactory.h"
namespace zaran
{

	NSSolver::NSSolver(int index, string name, FlowSolverPara *para, GridBase *grid, DataManagerNS *data_manager)
		: FlowFieldSolver(index, name, para, grid, data_manager)
	{
	}

	NSSolver::~NSSolver()
	{
		delete[] m_riemann_solver;
		delete[] m_gas;
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
			InitFieldFarFieldNoVelocity();
		}
		else if (init_type == InitFieldType::FarFlow)
		{
			InitFieldFarFlow();
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
	DataManagerNS *NSSolver::GetDataManager()
	{
		return static_cast<DataManagerNS *>(FieldSolver::GetDataManager());
	}
	void NSSolver::InitSolver()
	{
		Log::info("Initialize NS Solver!");
		FlowFieldSolver::InitSolver();
		RiemannSolverBuilder riemann_solver_builder;
		m_riemann_solver = riemann_solver_builder.Create(GetPara()->GetRiemannSolverType());
		auto ref_value = GetPara()->GetDimensionless();
		m_gas = new PerfectGas(ref_value.GetRefMw(), ref_value.GetRefGamma(), ref_value);
		Log::info("NS Solver Initialize Finished!");
	}

	void NSSolver::Preprocess()
	{
		auto para = GetPara();
		para->SetCurrentStep(para->GetCurrentStep() + 1);
		CalcTimeStep();
	}

	void NSSolver::Postprocess()
	{
		UpdateField();
		BoundaryCondition();
		CheckPrimtive();
		FixPrimtive();
		CalcForce();
	}

	void NSSolver::Solve()
	{
		TimeAdvance();
	}

	void NSSolver::CalcTimeStep()
	{
		CalcTimeStepLocal();
		double dt = GlobalData::GetDouble("dt");
		double current_time = GlobalData::GetDouble("currentTime");
		double end_time = GlobalData::GetDouble("endTime");
		if (current_time + dt > end_time)
		{
			dt = end_time - current_time;
			current_time = end_time;
		}
		else
			current_time += dt;
		GlobalData::Update("currentTime", current_time);
		int isSteady = GlobalData::GetInt("isSteady");
		if (isSteady == 0)
		{
			ReduceTimeStep(dt);
		}
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