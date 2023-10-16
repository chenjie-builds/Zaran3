#include "vanleer.h"
using namespace zaran;
void Vanleer::Solver(Ptr<RiemannSolverPara>& para)
{
	double& rhoL = para->primL(0);
	double& uL = para->primL(1);
	double& vL = para->primL(2);
	double& wL = para->primL(3);
	double& pL = para->primL(4);
	double& rhoR = para->primR(0);
	double& uR = para->primR(1);
	double& vR = para->primR(2);
	double& wR = para->primR(3);
	double& pR = para->primR(4);
	double deltaNorm = para->norm.norm();
	//¼ÆËãFL+
	para->cL = sqrt(para->gammaL * pL / rhoL);
	para->v2L = uL * uL + vL * vL + wL * wL;
	para->vnL = uL * para->norm(0) + vL * para->norm(1) + wL * para->norm(2) + para->nt;
	para->machL = para->vnL / (para->cL * deltaNorm);
	para->eL = 0.5 * para->v2L + pL / (rhoL * (para->gammaL - 1));
	para->flux.setZero();
	if (para->machL <= -1.0)
	{
		//do nothing
	}
	else if (para->machL >= 1.0)
	{
		para->flux[0] += rhoL * para->vnL;
		para->flux[1] += rhoL * para->vnL * uL + pL * para->norm(0);
		para->flux[2] += rhoL * para->vnL * vL + pL * para->norm(1);
		para->flux[3] += rhoL * para->vnL * wL + pL * para->norm(2);
		para->flux[4] += para->vnL * (rhoL * para->eL + pL) - pL * para->nt;
	}
	else
	{
		double fluxMass = 0.25 * rhoL * para->cL * (para->machL + 1) * (para->machL + 1);
		para->flux[0] += fluxMass * deltaNorm;
		para->flux[1] += fluxMass * deltaNorm * (para->norm(0) / deltaNorm * (-para->vnL / deltaNorm + 2 * para->cL) / para->gammaL + uL);
		para->flux[2] += fluxMass * deltaNorm * (para->norm(1) / deltaNorm * (-para->vnL / deltaNorm + 2 * para->cL) / para->gammaL + vL);
		para->flux[3] += fluxMass * deltaNorm * (para->norm(2) / deltaNorm * (-para->vnL / deltaNorm + 2 * para->cL) / para->gammaL + wL);
		para->flux[4] += fluxMass * deltaNorm * (para->vnL / deltaNorm * (-para->vnL / deltaNorm + 2 * para->cL) / (para->gammaL + 1)
			+ 2 * para->cL * para->cL / (para->gammaL * para->gammaL - 1)
			+ 0.5 * para->v2L - para->nt * (-para->vnL / deltaNorm + 2 * para->cL) / para->gammaL);
	}
	//¼ÆËãFR-
	para->cR = sqrt(para->gammaR * pR / rhoR);
	para->v2R = uR * uR + vR * vR + wR * wR;
	para->vnR = uR * para->norm(0) + vR * para->norm(1) + wR * para->norm(2) + para->nt;
	para->machR = para->vnR / (para->cR * deltaNorm);
	para->eR = 0.5 * para->v2R + pR / (rhoR * (para->gammaR - 1));
	if (para->machR <= -1.0)
	{
		para->flux[0] += rhoR * para->vnR;
		para->flux[1] += rhoR * para->vnR * uR + pR * para->norm(0);
		para->flux[2] += rhoR * para->vnR * vR + pR * para->norm(1);
		para->flux[3] += rhoR * para->vnR * wR + pR * para->norm(2);
		para->flux[4] += para->vnR * (rhoR * para->eR + pR) - pR * para->nt;
	}
	else if (para->machR >= 1.0)
	{
		//do nothing
	}
	else
	{
		double fluxMass = -0.25 * rhoR * para->cR * (para->machR - 1) * (para->machR - 1);
		para->flux[0] += fluxMass * deltaNorm;
		para->flux[1] += fluxMass * deltaNorm * (para->norm(0) / deltaNorm * (-para->vnR / deltaNorm - 2 * para->cR) / para->gammaR + uR);
		para->flux[2] += fluxMass * deltaNorm * (para->norm(1) / deltaNorm * (-para->vnR / deltaNorm - 2 * para->cR) / para->gammaR + vR);
		para->flux[3] += fluxMass * deltaNorm * (para->norm(2) / deltaNorm * (-para->vnR / deltaNorm - 2 * para->cR) / para->gammaR + wR);
		para->flux[4] += fluxMass * deltaNorm * (para->vnR / deltaNorm * (-para->vnR / deltaNorm - 2 * para->cR) / (para->gammaR + 1)
			+ 2 * para->cR * para->cR / (para->gammaR * para->gammaR - 1)
			+ 0.5 * para->v2R - para->nt * (-para->vnR / deltaNorm - 2 * para->cR) / para->gammaR);
	}
}