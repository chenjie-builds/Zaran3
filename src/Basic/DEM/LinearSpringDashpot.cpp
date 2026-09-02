#include "LinearSpringDashpot.h"
#include "CommonPara.h"
#include <algorithm>

namespace zaran
{

void LinearSpringDashpot::CalcNormalStiffness(const DEMParticle& pa, const DEMParticle& pb,
                                               double& k_n, double& c_n) const
{
    // 等效杨氏模量：1/E* = (1-νa²)/Ea + (1-νb²)/Eb
    double E_a = pa.young_modulus * (1.0 - pa.poisson_ratio * pa.poisson_ratio);
    double E_b = pb.young_modulus * (1.0 - pb.poisson_ratio * pb.poisson_ratio);
    double E_star = 1.0 / (1.0 / E_a + 1.0 / E_b);

    // 等效半径
    double R_star = (pa.radius * pb.radius) / (pa.radius + pb.radius);

    // 法向刚度（简单线弹簧使用弦接触宽度的量级）
    k_n = 2.0 * E_star * R_star;

    // 等效质量
    double m_eff = (pa.mass * pb.mass) / (pa.mass + pb.mass);

    // 阻尼系数 c = 2 * beta * sqrt(m_eff * k_n)，使用恢复系数换算 beta
    // beta = -ln(e) / sqrt(pi² + ln²(e))
    double e = std::min(pa.restitution_coeff, pb.restitution_coeff);
    e = std::max(e, 1.0e-3); // 防止 e=0
    double ln_e = std::log(e);
    double beta = -ln_e / std::sqrt(PI * PI + ln_e * ln_e);
    c_n = 2.0 * beta * std::sqrt(m_eff * k_n);
}

void LinearSpringDashpot::CalcNormalForce(const DEMParticle& pa, const DEMParticle& pb,
                                           DEMContact& contact, double /*dt*/)
{
    double k_n = 0.0, c_n = 0.0;
    CalcNormalStiffness(pa, pb, k_n, c_n);

    // 法向相对速度 v_n_rel = (vel_a - vel_b) · n
    Eigen::Vector3d rel_vel = pa.vel - pb.vel;
    double v_n_rel = rel_vel.dot(contact.normal);

    // F_n = -(k_n * δ + c_n * v_n_rel) * n  (排斥力，方向沿 -normal 作用于 A)
    double Fn_mag = k_n * contact.overlap_n - c_n * v_n_rel;
    if (Fn_mag < 0.0) Fn_mag = 0.0; // 不允许拉力

    contact.force_n = -Fn_mag * contact.normal;
}

void LinearSpringDashpot::CalcTangentialForce(const DEMParticle& pa, const DEMParticle& pb,
                                               DEMContact& contact, double dt)
{
    double k_n = 0.0, c_n = 0.0;
    CalcNormalStiffness(pa, pb, k_n, c_n);
    double k_t = 0.5 * k_n; // 切向刚度取法向的 1/2（常用经验值）

    // 切向相对速度（去除法向分量）
    Eigen::Vector3d rel_vel = pa.vel - pb.vel;
    // 旋转贡献
    rel_vel += (pa.radius * pa.omega + pb.radius * pb.omega).cross(contact.normal);
    Eigen::Vector3d v_t = rel_vel - rel_vel.dot(contact.normal) * contact.normal;

    // 增量切向位移
    contact.delta_t += v_t * dt;
    // 将 delta_t 投影回切向平面（防止坐标旋转导致的漂移）
    contact.delta_t -= contact.delta_t.dot(contact.normal) * contact.normal;

    Eigen::Vector3d Ft = -k_t * contact.delta_t;

    // Coulomb 摩擦截断
    double mu = std::min(pa.friction_coeff, pb.friction_coeff);
    double Fn_mag = contact.force_n.norm();
    double Ft_max = mu * Fn_mag;
    if (Ft.norm() > Ft_max)
    {
        Ft = Ft.normalized() * Ft_max;
        // 滑动时重置弹簧位移
        contact.delta_t = -Ft / k_t;
    }

    contact.force_t = Ft;
}

} // namespace zaran
