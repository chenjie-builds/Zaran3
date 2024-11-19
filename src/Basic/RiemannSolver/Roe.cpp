#include "Roe.h"
#include "MathBasic.h"
#include "Log.h"
using namespace zaran;
void Roe::Solver(RiemannSolverPara &para)
{
    double &rhoL = para.prim_left(0);
    double &uL = para.prim_left(1);
    double &vL = para.prim_left(2);
    double &wL = para.prim_left(3);
    double &pL = para.prim_left(4);
    double deltaNorm = para.norm.norm();
    para.norm /= deltaNorm;
    para.nt /= deltaNorm;
    para.c_left = sqrt(para.gamma_left * pL / rhoL);
    para.v2_left = uL * uL + vL * vL + wL * wL;
    para.vn_left = uL * para.norm(0) + vL * para.norm(1) + wL * para.norm(2) + para.nt;
    para.energy_left = 0.5 * para.v2_left + pL / (rhoL * (para.gamma_left - 1));
    para.enthalpy_left = para.energy_left + pL / rhoL;
    // FL
    para.flux[0] = rhoL * para.vn_left;
    para.flux[1] = rhoL * para.vn_left * uL + pL * para.norm(0);
    para.flux[2] = rhoL * para.vn_left * vL + pL * para.norm(1);
    para.flux[3] = rhoL * para.vn_left * wL + pL * para.norm(2);
    para.flux[4] = para.vn_left * (rhoL * para.energy_left + pL) - pL * para.nt;

    double &rhoR = para.prim_right(0);
    double &uR = para.prim_right(1);
    double &vR = para.prim_right(2);
    double &wR = para.prim_right(3);
    double &pR = para.prim_right(4);
    para.c_right = sqrt(para.gamma_right * pR / rhoR);
    para.v2_right = uR * uR + vR * vR + wR * wR;
    para.vn_right = uR * para.norm(0) + vR * para.norm(1) + wR * para.norm(2) + para.nt;
    para.energy_right = 0.5 * para.v2_right + pR / (rhoR * (para.gamma_right - 1));
    para.enthalpy_right = para.energy_right + pR / rhoR;
    // FR
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

    // roe平均
    rho_roe = sqrt(rhoL * rhoR);
    u_roe = (uL + uR * sqrt(rhoR / rhoL)) / (1 + sqrt(rhoR / rhoL));
    v_roe = (vL + vR * sqrt(rhoR / rhoL)) / (1 + sqrt(rhoR / rhoL));
    w_roe = (wL + wR * sqrt(rhoR / rhoL)) / (1 + sqrt(rhoR / rhoL));
    enthalpy_roe = (para.enthalpy_left + para.enthalpy_right * sqrt(rhoR / rhoL)) / (1 + sqrt(rhoR / rhoL));
    c_roe = (para.gamma_left - 1) * (enthalpy_roe - 0.5 * (u_roe * u_roe + v_roe * v_roe + w_roe * w_roe));
    c_roe = sqrt(c_roe);
    vn_roe = u_roe * para.norm(0) + v_roe * para.norm(1) + w_roe * para.norm(2) + para.nt;
    double v2_roe = u_roe * u_roe + v_roe * v_roe + w_roe * w_roe;
    delta_vn = para.vn_right - para.vn_left;
    double lamda1 = abs(vn_roe - c_roe);
    double lamda3 = abs(vn_roe);
    double lamda5 = abs(vn_roe + c_roe);
    // Harten(lamda1, 2.0 * Max(para.vn_right - para.vn_left, 0.0));
    // Harten(lamda2, 2.0 * Max(para.vn_right + para.c_right - para.vn_left - para.c_left, 0.0));
    // Harten(lamda3, 2.0 * Max(para.vn_right - para.c_right - para.vn_left + para.c_left, 0.0));
    Harten(lamda1, 0.01);
    Harten(lamda3, 0.01);
    Harten(lamda5, 0.01);
    lamda1 *= (delta_p - rho_roe * c_roe * delta_vn) / (2 * c_roe * c_roe);
    lamda5 *= (delta_p + rho_roe * c_roe * delta_vn) / (2 * c_roe * c_roe);
    para.flux[0] -= lamda1 + lamda5;
    para.flux[1] -= lamda1 * (u_roe - para.norm(0) * c_roe) + lamda5 * (u_roe + para.norm(0) * c_roe);
    para.flux[2] -= lamda1 * (v_roe - para.norm(1) * c_roe) + lamda5 * (v_roe + para.norm(1) * c_roe);
    para.flux[3] -= lamda1 * (w_roe - para.norm(2) * c_roe) + lamda5 * (w_roe + para.norm(2) * c_roe);
    para.flux[4] -= lamda1 * (enthalpy_roe - c_roe * vn_roe) + lamda5 * (enthalpy_roe + c_roe * vn_roe);

    double dlp = lamda3 * (delta_rho - delta_p / (c_roe * c_roe));
    double dlr = lamda3 * rho_roe;
    para.flux[0] -= dlp;
    para.flux[1] -= dlp * u_roe + dlr * (delta_u - delta_vn * para.norm(0));
    para.flux[2] -= dlp * v_roe + dlr * (delta_v - delta_vn * para.norm(1));
    para.flux[3] -= dlp * w_roe + dlr * (delta_w - delta_vn * para.norm(2));
    para.flux[4] -= dlp * (0.5 * v2_roe) + dlr * (u_roe * delta_u + v_roe * delta_v + w_roe * delta_w - delta_vn * vn_roe);

    para.flux[0] *= 0.5 * deltaNorm;
    para.flux[1] *= 0.5 * deltaNorm;
    para.flux[2] *= 0.5 * deltaNorm;
    para.flux[3] *= 0.5 * deltaNorm;
    para.flux[4] *= 0.5 * deltaNorm;
    if (abs(rhoL - 1.0) < 1e-8 && abs(rhoR - 0.125) < 1e-8)
    {
        Log::info("Roe: prim_left: {},{},{},{},{},{} prim_right: {},{},{},{},{},{} flux: {},{},{},{},{}", para.prim_left(0), para.prim_left(1), para.prim_left(2), para.prim_left(3), para.prim_left(4), para.prim_left(5), para.prim_right(0), para.prim_right(1), para.prim_right(2), para.prim_right(3), para.prim_right(4), para.prim_right(5), para.flux(0), para.flux(1), para.flux(2), para.flux(3), para.flux(4));
        exit(0);
    }

}

void zaran::Roe::Harten(double &lamda, double delta)
{
    if (lamda <= delta)
    {
        lamda = (lamda * lamda + delta * delta) / (2 * delta);
    }
}
