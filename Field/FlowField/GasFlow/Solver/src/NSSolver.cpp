#include "NSSolver.h"
#include "File.h"
#include <corecrt_math.h>
#include <fstream>
namespace zaran
{

	NSSolver::NSSolver(int index, string name, FlowSolverPara* para, GridBase* grid, FieldData* fieldData)
		:FlowSolver(index, name, para, grid, fieldData)
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
		CalcMetric();
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
		else
		{
			Log::warn("Initialize Failed!");
			Log::warn("Wrong Flow field init parameter:{}", GlobalData::GetString("initFieldType"));
			exit(0);
		}
		Prim2Cons();
		Log::info("Flow Field Initialize Finished!");
	}
	void NSSolver::InitSolver()
	{
		Log::info("Initialize NS Solver!");
		FlowSolver::InitSolver();
		RiemannSolverBuilder riemannSolverFactory;
		m_riemann_solver = riemannSolverFactory.Create(GetPara()->GetRiemannSolverType());
		auto ref_value = GetPara()->GetDimensionless();
		m_gas=new PerfectGas(ref_value.GetRefMw(), ref_value.GetRefGamma(), ref_value);
		Log::info("NS Solver Initialize Finished!");
	}

	void NSSolver::Preprocess()
	{
		CalcTimeStep();
		CalcPrimGrad();
		CalcLimiter();
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
		BoundaryCondition();
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
			SnycTimeStepWithGlobal(dt);
		}
	}
	void NSSolver::TimeAdvance()
	{
		RungeKutta();
	}

	void NSSolver::CalcPrimGrad()
	{
		auto para = GetPara();
		if (para->GetGradScheme() == GradScheme::wls)
		{
			CalcGradWLS();
		}
		else if (para->GetGradScheme() == GradScheme::ufdm)
		{
			CalcGradUFDM();
		}
		else if (para->GetGradScheme() == GradScheme::noGrad)
		{
			NoGradient();
		}
		else
		{
			Log::warn("Unsupported Gradiend Scheme!");
		}
		CalcPrimGradBound();
	}

	void NSSolver::CalcResidual()
	{
		ZeroResidual();
		ConvectiveResidual();
		ViscousResidual();
		SourceTermResidual();
	}

	void NSSolver::CalcGradUFDM()
	{
		Log::warn("TO DO Gradient Function UFDM!");
	}

	void NSSolver::UpdateField()
	{
		Cons2Prim();
	}



	void NSSolver::NoGradient()
	{
		// do nothing
	}

} // namespace zaran