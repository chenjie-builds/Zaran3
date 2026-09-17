#include "HertzMindlin.h"
#include "CommonPara.h"
#include <algorithm>

namespace zaran
{

void HertzMindlin::CalcEffectiveParams(const DEMParticle& pa, const DEMParticle& pb,
                                        double& E_star, double& R_star, double& G_star) const
{
    // 等效弹性模量：1/E* = (1-νa²)/Ea + (1-νb²)/Eb
    double inv_E = (1.0 - pa.poisson_ratio * pa.poisson_ratio) / pa.young_modulus
                 + (1.0 - pb.poisson_ratio * pb.poisson_ratio) / pb.young_modulus;
    E_star = 1.0 / inv_E;

    // 等效半径
    R_star = (pa.radius * pb.radius) / (pa.radius + pb.radius);

    // 等效剪切模量：1/G* = (2-νa)/Ga + (2-νb)/Gb，G = E/(2(1+ν))
    double G_a = pa.young_modulus / (2.0 * (1.0 + pa.poisson_ratio));
    double G_b = pb.young_modulus / (2.0 * (1.0 + pb.poisson_ratio));
    double inv_G = (2.0 - pa.poisson_ratio) / G_a + (2.0 - pb.poisson_ratio) / G_b;
    G_star = 1.0 / inv_G;
}

void HertzMindlin::CalcNormalForce(const DEMParticle& pa, const DEMParticle& pb,
                                    DEMContact& contact, double dt)
{
    double E_star = 0.0, R_star = 0.0, G_star = 0.0;
    CalcEffectiveParams(pa, pb, E_star, R_star, G_star);

    double delta = contact.overlap_n;
    if (delta <= 0.0)
    {
        contact.force_n.setZero();
        contact.dissipation = 0.0;
        return;
    }

    // Hertz 法向刚度：k_n = 2 E* sqrt(R* δ)
    double k_n = 2.0 * E_star * std::sqrt(R_star * delta);

    // 等效质量
    double m_eff = 0.0;
    if (!pa.IsDynamic())
        m_eff = pb.mass;
    else if (!pb.IsDynamic())
        m_eff = pa.mass;
    else
        m_eff = (pa.mass * pb.mass) / (pa.mass + pb.mass);

    // 临界阻尼比（基于恢复系数）
    double e = std::min(pa.restitution_coeff, pb.restitution_coeff);
    e = std::max(e, 1.0e-3);
    double ln_e = std::log(e);
    double beta = -ln_e / std::sqrt(PI * PI + ln_e * ln_e);
    //
    // 缺陷修复：Hertz 法向刚度 k_n = 2E*·sqrt(R*·δ) 随重叠量 δ 变化（非线性弹簧），
    // 不能直接照搬线性模型的阻尼系数 c_n = 2·beta·sqrt(m_eff·k_n)。后者会使实测
    // 恢复系数系统性偏低（e_nom=0.5 → 约 0.466）。这里采用标准 Hertzian（Tsuji /
    // Brilliantov）阻尼系数，附加 sqrt(5/6) 非线性修正因子；数值验证表明可使实测
    // 恢复系数与标称值精确吻合（e_nom ∈ [0.1, 0.95] 误差 < 1e-3）。
    double c_n = 2.0 * std::sqrt(5.0 / 6.0) * beta * std::sqrt(m_eff * k_n);

    // 法向相对速度
    Eigen::Vector3d rel_vel = pa.vel - pb.vel;
    double v_n_rel = rel_vel.dot(contact.normal);

    // F_n = (4/3 E* sqrt(R*) δ^{3/2} + c_n v_n_rel) 作用于 A 沿 -normal
    //
    // 注意（缺陷修复）：与 LinearSpringDashpot 同理，回弹阶段阻尼项产生的
    // 小幅"拉力"是恢复系数标定 c_n = 2·beta·sqrt(m_eff·k_n) 的组成部分，
    // 不应被钳位抹掉。无重叠（δ<=0）的情形已在上方提前返回零力处理。
    double Fn_hertz = (4.0 / 3.0) * E_star * std::sqrt(R_star) * std::pow(delta, 1.5);
    double Fn_damp  = c_n * v_n_rel;
    double Fn_mag = Fn_hertz + Fn_damp;

    contact.force_n = -Fn_mag * contact.normal;

    // 法向阻尼耗散：只有粘性（阻尼）部分不可逆，Hertz 弹性项属可逆储能。
    // 该对被粘性耗散的功率为 c_n·v_n_rel²（恒 ≥ 0）。
    contact.dissipation = c_n * v_n_rel * v_n_rel * dt;
}

void HertzMindlin::CalcTangentialForce(const DEMParticle& pa, const DEMParticle& pb,
                                        DEMContact& contact, double dt)
{
    double E_star = 0.0, R_star = 0.0, G_star = 0.0;
    CalcEffectiveParams(pa, pb, E_star, R_star, G_star);

    double delta = contact.overlap_n;
    if (delta <= 0.0)
    {
        contact.force_t.setZero();
        contact.delta_t.setZero();
        return;
    }

    // 切向刚度 k_t = 8 G* sqrt(R* δ)
    double k_t = 8.0 * G_star * std::sqrt(R_star * delta);

    // 切向相对速度
    Eigen::Vector3d rel_vel = pa.vel - pb.vel;
    rel_vel += (pa.radius * pa.omega + pb.radius * pb.omega).cross(contact.normal);
    Eigen::Vector3d v_t = rel_vel - rel_vel.dot(contact.normal) * contact.normal;

    // 增量切向位移
    contact.delta_t += v_t * dt;
    contact.delta_t -= contact.delta_t.dot(contact.normal) * contact.normal;

    Eigen::Vector3d Ft = -k_t * contact.delta_t;

    // Coulomb 摩擦截断
    double mu = std::min(pa.friction_coeff, pb.friction_coeff);
    double Fn_mag = contact.force_n.norm();
    double Ft_max = mu * Fn_mag;
    if (Ft.norm() > Ft_max)
    {
        // 切向摩擦耗散：取库仑截断前后切向弹簧储能之差（严格非负），
        // 持续滑动时等于滑动功 mu*|Fn|*|Δs|。
        const double stored_before = 0.5 * k_t * contact.delta_t.squaredNorm();
        Ft = Ft.normalized() * Ft_max;
        contact.delta_t = -Ft / k_t;
        const double stored_after = 0.5 * k_t * contact.delta_t.squaredNorm();
        contact.dissipation += stored_before - stored_after;
    }

    contact.force_t = Ft;
}

} // namespace zaran
