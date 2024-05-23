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

void NSSolver::Init()
	{
		InitSolver();
		InitField();
		CalcMetric();
	}
	void NSSolver::InitField()
	{
		FlowSolverPara* para = GetPara();
		const InitFieldType& init_type = para->GetInitFieldType();
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
		SetEquNum(5);
		FlowSolver::InitSolver();
		std::string riemannSolverType = GlobalData::GetString("riemannSolver");
		RiemannSolverFactory riemannSolverFactory;
		riemannSolverFactory.Create(m_riemann_solver, riemannSolverType);
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
		double endTime = GlobalData::GetDouble("endTime");
		if (current_time + dt > endTime)
		{
			dt = endTime - current_time;
			current_time = endTime;
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
		FlowSolverPara* para = GetPara();
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

	void NSSolver::Prim2Cons(const double* prim, double* cons)
	{
		FlowSolverPara* para = GetPara();
		Gas* gas = para->GetGas();
		double gamma = gas->GetGamma();
		double v2 = prim[1] * prim[1] + prim[2] * prim[2] + prim[3] * prim[3];
		cons[0] = prim[0];
		cons[1] = prim[0] * prim[1];
		cons[2] = prim[0] * prim[2];
		cons[3] = prim[0] * prim[3];
		cons[4] = 0.5 * prim[0] * v2 + prim[4] / (gamma - 1);
	}

	void NSSolver::Cons2Prim(const double* cons, double* prim)
	{
		FlowSolverPara* para = GetPara();
		Gas* gas = para->GetGas();
		double gamma = gas->GetGamma();
		prim[0] = cons[0];
		prim[1] = cons[1] / cons[0];
		prim[2] = cons[2] / cons[0];
		prim[3] = cons[3] / cons[0];
		double v2 = prim[1] * prim[1] + prim[2] * prim[2] + prim[3] * prim[3];
		prim[4] = (cons[4] - 0.5 * cons[0] * v2) * (gamma - 1);
	}

} // namespace zaran