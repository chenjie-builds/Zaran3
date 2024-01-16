#include "vanleer.h"
#include "Log.h"
using namespace zaran;
void Vanleer::Solver(RiemannSolverPara& para)
{
	para.gamma_left = para.gamma_right = 1.4;
	double& rhoL = para.prim_left(0);
	double& uL = para.prim_left(1);
	double& vL = para.prim_left(2);
	double& wL = para.prim_left(3);
	double& pL = para.prim_left(4);
	double& rhoR = para.prim_right(0);
	double& uR = para.prim_right(1);
	double& vR = para.prim_right(2);
	double& wR = para.prim_right(3);
	double& pR = para.prim_right(4);
	double deltaNorm = para.norm.norm();
	//FL+
	para.c_left = sqrt(para.gamma_left * pL / rhoL);
	para.v2_left = uL * uL + vL * vL + wL * wL;
	para.vn_left = uL * para.norm(0) + vL * para.norm(1) + wL * para.norm(2) + para.nt;
	para.mach_left = para.vn_left / (para.c_left * deltaNorm);
	para.energy_left = 0.5 * para.v2_left + pL / (rhoL * (para.gamma_left - 1));
	para.flux.setZero();
	if (para.mach_left <= -1.0)
	{
		//do nothing
	}
	else if (para.mach_left >= 1.0)
	{
		para.flux[0] += rhoL * para.vn_left;
		para.flux[1] += rhoL * para.vn_left * uL + pL * para.norm(0);
		para.flux[2] += rhoL * para.vn_left * vL + pL * para.norm(1);
		para.flux[3] += rhoL * para.vn_left * wL + pL * para.norm(2);
		para.flux[4] += para.vn_left * (rhoL * para.energy_left + pL) - pL * para.nt;
	}
	else
	{
		double fluxMass = 0.25 * rhoL * para.c_left * (para.mach_left + 1) * (para.mach_left + 1);
		para.flux[0] += fluxMass * deltaNorm;
		para.flux[1] += fluxMass * deltaNorm * (para.norm(0) / deltaNorm * (-para.vn_left / deltaNorm + 2 * para.c_left) / para.gamma_left + uL);
		para.flux[2] += fluxMass * deltaNorm * (para.norm(1) / deltaNorm * (-para.vn_left / deltaNorm + 2 * para.c_left) / para.gamma_left + vL);
		para.flux[3] += fluxMass * deltaNorm * (para.norm(2) / deltaNorm * (-para.vn_left / deltaNorm + 2 * para.c_left) / para.gamma_left + wL);
		para.flux[4] += fluxMass * deltaNorm * (para.vn_left / deltaNorm * (-para.vn_left / deltaNorm + 2 * para.c_left) / (para.gamma_left + 1)
			+ 2 * para.c_left * para.c_left / (para.gamma_left * para.gamma_left - 1)
			+ 0.5 * para.v2_left - para.nt * (-para.vn_left / deltaNorm + 2 * para.c_left) / para.gamma_left);
	}
	//FR-
	para.c_right = sqrt(para.gamma_right * pR / rhoR);
	para.v2_right = uR * uR + vR * vR + wR * wR;
	para.vn_right = uR * para.norm(0) + vR * para.norm(1) + wR * para.norm(2) + para.nt;
	para.mach_right = para.vn_right / (para.c_right * deltaNorm);
	para.energy_right = 0.5 * para.v2_right + pR / (rhoR * (para.gamma_right - 1));
	if (para.mach_right <= -1.0)
	{
		para.flux[0] += rhoR * para.vn_right;
		para.flux[1] += rhoR * para.vn_right * uR + pR * para.norm(0);
		para.flux[2] += rhoR * para.vn_right * vR + pR * para.norm(1);
		para.flux[3] += rhoR * para.vn_right * wR + pR * para.norm(2);
		para.flux[4] += para.vn_right * (rhoR * para.energy_right + pR) - pR * para.nt;
	}
	else if (para.mach_right >= 1.0)
	{
		//do nothing
	}
	else
	{
		double fluxMass = -0.25 * rhoR * para.c_right * (para.mach_right - 1) * (para.mach_right - 1);
		para.flux[0] += fluxMass * deltaNorm;
		para.flux[1] += fluxMass * deltaNorm * (para.norm(0) / deltaNorm * (-para.vn_right / deltaNorm - 2 * para.c_right) / para.gamma_right + uR);
		para.flux[2] += fluxMass * deltaNorm * (para.norm(1) / deltaNorm * (-para.vn_right / deltaNorm - 2 * para.c_right) / para.gamma_right + vR);
		para.flux[3] += fluxMass * deltaNorm * (para.norm(2) / deltaNorm * (-para.vn_right / deltaNorm - 2 * para.c_right) / para.gamma_right + wR);
		para.flux[4] += fluxMass * deltaNorm * (para.vn_right / deltaNorm * (-para.vn_right / deltaNorm - 2 * para.c_right) / (para.gamma_right + 1)
			+ 2 * para.c_right * para.c_right / (para.gamma_right * para.gamma_right - 1)
			+ 0.5 * para.v2_right - para.nt * (-para.vn_right / deltaNorm - 2 * para.c_right) / para.gamma_right);
	}

}