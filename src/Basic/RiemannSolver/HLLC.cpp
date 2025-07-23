#include "HLLC.h"
#include "Log.h"
#include "MathBasic.h"
using namespace zaran;
// HLLC
void HLLC::Solver(RiemannSolverPara &para)
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

    double &rhoR = para.prim_right(0);
    double &uR = para.prim_right(1);
    double &vR = para.prim_right(2);
    double &wR = para.prim_right(3);
    double &pR = para.prim_right(4);
    para.c_right = sqrt(para.gamma_right * pR / rhoR);
    para.v2_right = uR * uR + vR * vR + wR * wR;
    para.vn_right = uR * para.norm(0) + vR * para.norm(1) + wR * para.norm(2) + para.nt;
    para.energy_right = 0.5 * para.v2_right + pR / (rhoR * (para.gamma_right - 1));

    double dl = sqrt(rhoL);
    double dr = sqrt(rhoR);
    double ita = 0.5 * dl * dr / pow(dl + dr, 2); // 计算ita
    double u_roe = (dl * uL + dr * uR) / (dl + dr);
    double v_roe = (dl * vL + dr * vR) / (dl + dr);
    double w_roe = (dl * wL + dr * wR) / (dl + dr);
    double vn_roe = u_roe * para.norm(0) + v_roe * para.norm(1) + w_roe * para.norm(2) + para.nt; // 速度的roe平均值
    double c_roe = (dl * para.c_left * para.c_left + dr * para.c_right * para.c_right) / (dl + dr) +
                   ita * (std::pow(uL - uR, 2) + std::pow(vL - vR, 2) + std::pow(wL - wR, 2)); // 声速的roe平均值
    c_roe = sqrt(c_roe); // 计算声速的roe平均值
    double s_left = Min(para.vn_left - para.c_left, vn_roe - c_roe);
    double s_right = Max(para.vn_right + para.c_right, vn_roe + c_roe);
    // double s_left = para.vn_left - para.c_left;
    // double s_right = para.vn_right + para.c_right;
    double sl = rhoL * (s_left - para.vn_left);
    double sr = rhoR * (s_right - para.vn_right);
    double s_middle = (sr * para.vn_right - sl * para.vn_left + pL - pR) / (sr - sl);

    double p_star = -sl * (para.vn_left - s_middle) + pL;
    double p_star_R = -sr * (para.vn_right - s_middle) + pR;
    if (abs(p_star - p_star_R) > 1e-10)
    {
        Log::error("HLLC: p_star and p_star_R are not equal, p_star = {}, p_star_R = {}", p_star, p_star_R);
        Log::error("HLLC: s_left = {}, s_right = {}, s_middle = {}", s_left, s_right, s_middle);
        Log::error("HLLC: vn_left = {}, vn_right = {}, c_left = {}, c_right = {}", para.vn_left, para.vn_right,
                   para.c_left, para.c_right);
        Log::error("HLLC: rhoL = {}, rhoR = {}, uL = {}, uR = {}, vL = {}, vR = {}, wL = {}, wR = {}, pL = {}, pR = {}",
                   rhoL, rhoR, uL, uR, vL, vR, wL, wR, pL, pR);
        Log::error("HLLC: norm = {}, {}, {}, nt = {}", para.norm(0), para.norm(1), para.norm(2), para.nt);
        Log::error("HLLC: gamma_left = {}, gamma_right = {}", para.gamma_left, para.gamma_right);
        Log::error("HLLC: energy_left = {}, energy_right = {}", para.energy_left, para.energy_right);
        Log::error("HLLC: vn_left = {}, vn_right = {}, v2_left = {}, v2_right = {}", para.vn_left, para.vn_right,
                   para.v2_left, para.v2_right);
        Log::error("HLLC: s_left = {}, s_right = {}, s_middle = {}", s_left, s_right, s_middle);
        Log::error("HLLC: rhoL = {}, rhoR = {}, uL = {}, uR = {}, vL = {}, vR = {}, wL = {}, wR = {}, pL = {}, pR = {}",
                   rhoL, rhoR, uL, uR, vL, vR, wL, wR, pL, pR);
        exit(0);
    }
    double cons_star_L[5];
    double par1 = (s_left - para.vn_left);
    double par2 = 1.0 / (s_left - s_middle);
    cons_star_L[0] = par2 * par1 * rhoL;
    cons_star_L[1] = par2 * par1 * rhoL * uL + par2 * para.norm(0) * (p_star - pL);
    cons_star_L[2] = par2 * par1 * rhoL * vL + par2 * para.norm(1) * (p_star - pL);
    cons_star_L[3] = par2 * par1 * rhoL * wL + par2 * para.norm(2) * (p_star - pL);
    cons_star_L[4] = par2 * par1 * rhoL * para.energy_left + par2 * (s_middle * p_star - para.vn_left * pL);
    double cons_star_R[5];
    par1 = (s_right - para.vn_right);
    par2 = 1.0 / (s_right - s_middle);
    cons_star_R[0] = par2 * par1 * rhoR;
    cons_star_R[1] = par2 * par1 * rhoR * uR + par2 * para.norm(0) * (p_star - pR);
    cons_star_R[2] = par2 * par1 * rhoR * vR + par2 * para.norm(1) * (p_star - pR);
    cons_star_R[3] = par2 * par1 * rhoR * wR + par2 * para.norm(2) * (p_star - pR);
    cons_star_R[4] = par2 * par1 * rhoR * para.energy_right + par2 * (s_middle * p_star - para.vn_right * pR);
    if (s_left > 0)
    {
        para.flux[0] = rhoL * para.vn_left;
        para.flux[1] = rhoL * para.vn_left * uL + pL * para.norm(0);
        para.flux[2] = rhoL * para.vn_left * vL + pL * para.norm(1);
        para.flux[3] = rhoL * para.vn_left * wL + pL * para.norm(2);
        para.flux[4] = para.vn_left * (rhoL * para.energy_left + pL) - pL * para.nt;
    }
    else if (s_left <= 0 && 0 < s_middle)
    {
        para.flux[0] = s_middle * cons_star_L[0];
        para.flux[1] = s_middle * cons_star_L[1] + p_star * para.norm(0);
        para.flux[2] = s_middle * cons_star_L[2] + p_star * para.norm(1);
        para.flux[3] = s_middle * cons_star_L[3] + p_star * para.norm(2);
        para.flux[4] = s_middle * cons_star_L[4] + s_middle * p_star - p_star * para.nt;
    }
    else if (s_middle <= 0 && 0 <= s_right)
    {
        para.flux[0] = s_middle * cons_star_R[0];
        para.flux[1] = s_middle * cons_star_R[1] + p_star * para.norm(0);
        para.flux[2] = s_middle * cons_star_R[2] + p_star * para.norm(1);
        para.flux[3] = s_middle * cons_star_R[3] + p_star * para.norm(2);
        para.flux[4] = s_middle * cons_star_R[4] + s_middle * p_star - p_star * para.nt;
    }
    else if (s_right <= 0)
    {
        para.flux[0] = rhoR * para.vn_right;
        para.flux[1] = rhoR * para.vn_right * uR + pR * para.norm(0);
        para.flux[2] = rhoR * para.vn_right * vR + pR * para.norm(1);
        para.flux[3] = rhoR * para.vn_right * wR + pR * para.norm(2);
        para.flux[4] = para.vn_right * (rhoR * para.energy_right + pR) - pR * para.nt;
    }
    for (int i = 0; i < 5; ++i)
    {
        para.flux[i] = para.flux[i] * deltaNorm;
    }
    // if (abs(rhoL - 1.0) < 1e-8 && abs(rhoR - 0.125) < 1e-8)
    // {
    //     Log::info("HLLC: prim_left: {},{},{},{},{},{} prim_right: {},{},{},{},{},{} flux: {},{},{},{},{}",
    //     para.prim_left(0), para.prim_left(1), para.prim_left(2), para.prim_left(3), para.prim_left(4),
    //     para.prim_left(5), para.prim_right(0), para.prim_right(1), para.prim_right(2), para.prim_right(3),
    //     para.prim_right(4), para.prim_right(5), para.flux(0), para.flux(1), para.flux(2), para.flux(3),
    //     para.flux(4)); exit(0);
    // }
    // if (isnan(para.flux(0)) || isnan(para.flux(1)) || isnan(para.flux(2)) || isnan(para.flux(3)) ||
    // isnan(para.flux(4)))
    // {
    //     Log::error("HLLC: flux is nan");
    //     Log::error("HLLC: prim_left: {},{},{},{},{} prim_right: {},{},{},{},{} flux: {},{},{},{},{}",
    //     para.prim_left(0), para.prim_left(1), para.prim_left(2), para.prim_left(3), para.prim_left(4),
    //     para.prim_right(0), para.prim_right(1), para.prim_right(2), para.prim_right(3), para.prim_right(4),
    //     para.flux(0), para.flux(1), para.flux(2), para.flux(3), para.flux(4)); Log::error("s_left: {}, s_star: {},
    //     s_right: {}", s_left, s_star, s_right); Log::error("rhoL: {}, rhoR: {}", rhoL, rhoR); Log::error("uL: {}, uR:
    //     {}", uL, uR); Log::error("vL: {}, vR: {}", vL, vR); Log::error("wL: {}, wR: {}", wL, wR); Log::error("pL: {},
    //     pR: {}", pL, pR); Log::error("gamma_left: {}, gamma_right: {}", para.gamma_left, para.gamma_right);
    //     Log::error("c_left: {}, c_right: {}, c_roe: {}", para.c_left, para.c_right, c_roe);
    //     Log::error("vn_left: {}, vn_right: {}, vn_roe: {}", para.vn_left, para.vn_right, vn_roe);
    //     Log::error("cons_star_L: {},{},{},{},{}", cons_star_L[0], cons_star_L[1], cons_star_L[2], cons_star_L[3],
    //     cons_star_L[4]); Log::error("cons_star_R: {},{},{},{},{}", cons_star_R[0], cons_star_R[1], cons_star_R[2],
    //     cons_star_R[3], cons_star_R[4]); Log::error("s1: {}, s2: {}, s_star: {}", s1, s2, s_star);
    //     Log::error("p_star: {}", p_star);
    //     exit(0);
    // }
}
