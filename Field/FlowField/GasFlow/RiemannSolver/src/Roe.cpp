#include "Roe.h"
using namespace zaran;
void Roe::Solver(RiemannSolverPara& para)
{
    double& rhoL = para.prim_left(0);
    double& uL = para.prim_left(1);
    double& vL = para.prim_left(2);
    double& wL = para.prim_left(3);
    double& pL = para.prim_left(4);
    double deltaNorm = para.norm.norm();
    para.norm /= deltaNorm;
    para.nt /= deltaNorm;
    para.c_left = sqrt(para.gamma_left * pL / rhoL);
    para.v2_left = uL * uL + vL * vL + wL * wL;
    para.vn_left = uL * para.norm(0) + vL * para.norm(1) + wL * para.norm(2) + para.nt;
    para.energy_left = 0.5 * para.v2_left + pL / (rhoL * (para.gamma_left - 1));
    para.enthalpy_left = para.energy_left + pL / rhoL;
    //FL
    para.flux[0] = rhoL * para.vn_left;
    para.flux[1] = rhoL * para.vn_left * uL + pL * para.norm(0);
    para.flux[2] = rhoL * para.vn_left * vL + pL * para.norm(1);
    para.flux[3] = rhoL * para.vn_left * wL + pL * para.norm(2);
    para.flux[4] = para.vn_left * (rhoL * para.energy_left + pL) - pL * para.nt;

    double& rhoR = para.prim_right(0);
    double& uR = para.prim_right(1);
    double& vR = para.prim_right(2);
    double& wR = para.prim_right(3);
    double& pR = para.prim_right(4);
    para.c_right = sqrt(para.gamma_right * pR / rhoR);
    para.v2_right = uR * uR + vR * vR + wR * wR;
    para.vn_right = uR * para.norm(0) + vR * para.norm(1) + wR * para.norm(2) + para.nt;
    para.energy_right = 0.5 * para.v2_right + pR / (rhoR * (para.gamma_right - 1));
    para.enthalpy_right = para.energy_right + pR / rhoR;
    //FR
    para.flux[0] += rhoR * para.vn_right;
    para.flux[1] += rhoR * para.vn_right * uR + pR * para.norm(0);
    para.flux[2] += rhoR * para.vn_right * vR + pR * para.norm(1);
    para.flux[3] += rhoR * para.vn_right * wR + pR * para.norm(2);
    para.flux[4] += para.vn_right * (rhoR * para.energy_right + pR) - pR * para.nt;

    delta_rho = rhoR - rhoL;
    delta_u = uR - uL;
    delta_v = vR - vL;
    delta_w = wR - wL;
    delta_p = pR - pL;

    //roe平均
    rho_roe = sqrt(rhoL * rhoR);
    u_roe = (uL + uR * sqrt(rhoR / rhoL)) / (1 + sqrt(rhoR / rhoL));
    v_roe = (vL + vR * sqrt(rhoR / rhoL)) / (1 + sqrt(rhoR / rhoL));
    w_roe = (wL + wR * sqrt(rhoR / rhoL)) / (1 + sqrt(rhoR / rhoL));
    enthalpy_roe = (para.enthalpy_left + para.enthalpy_right * sqrt(rhoR / rhoL)) / (1 + sqrt(rhoR / rhoL));
    c_roe = (para.gamma_left - 1) * (enthalpy_roe - 0.5 * (u_roe * u_roe + v_roe * v_roe + w_roe * w_roe));
    c_roe = sqrt(c_roe);
    vn_roe = u_roe * para.norm(0) + v_roe * para.norm(1) + w_roe * para.norm(2) + para.nt;
    delta_vn = para.vn_right - para.vn_left;
    a1 = abs(vn_roe) * (delta_rho - delta_p/ (c_roe * c_roe));
    a2 = abs(vn_roe + c_roe) * (delta_p + rho_roe * c_roe * delta_vn ) / (2 * c_roe * c_roe);
    a3 = abs(vn_roe - c_roe) * (delta_p - rho_roe * c_roe * delta_vn ) / (2 * c_roe * c_roe);
    a4 = a1 + a2 + a3;
    a5 = c_roe * (a2 - a3);
    a6 = abs(vn_roe) * (delta_rho * delta_u - para.norm(0) * delta_p * delta_vn );
    a7 = abs(vn_roe) * (delta_rho * delta_v - para.norm(1) * delta_p * delta_vn );
    a8 = abs(vn_roe) * (delta_rho * delta_w - para.norm(2) * delta_p * delta_vn );

    para.flux[0] -=  a4;
    para.flux[1] -=  a4 * u_roe + a5 * para.norm(0)  + a6;
    para.flux[2] -=  a4 * v_roe + a5 * para.norm(1)  + a7;
    para.flux[3] -=  a4 * w_roe + a5 * para.norm(2)  + a8;
    para.flux[4] -= a4 * enthalpy_roe + a5 * (vn_roe - para.nt) + a6 * u_roe + a7 * v_roe + a8 * w_roe - c_roe * c_roe * a1 / (para.gamma_left - 1);

    para.flux[0] *= 0.5*deltaNorm;
    para.flux[1] *= 0.5*deltaNorm;
    para.flux[2] *= 0.5*deltaNorm;
    para.flux[3] *= 0.5*deltaNorm;
    para.flux[4] *= 0.5*deltaNorm;

}
