#include "LinearSpringDashpot.h"
#include "CommonPara.h"
#include <algorithm>

namespace zaran
{

void LinearSpringDashpot::CalcNormalStiffness(const DEMParticle& pa, const DEMParticle& pb,
                                               double& k_n, double& c_n) const
{
    // 等效杨氏模量：1/E* = (1-νa²)/Ea + (1-νb²)/Eb
    double inv_E = (1.0 - pa.poisson_ratio * pa.poisson_ratio) / pa.young_modulus
                 + (1.0 - pb.poisson_ratio * pb.poisson_ratio) / pb.young_modulus;
    double E_star = 1.0 / inv_E;

    // 等效半径
    double R_star = (pa.radius * pb.radius) / (pa.radius + pb.radius);

    // 法向刚度（简单线弹簧使用弦接触宽度的量级）
    k_n = 2.0 * E_star * R_star;

    // 等效质量
    double m_eff = 0.0;
    if (!pa.IsDynamic())
        m_eff = pb.mass;
    else if (!pb.IsDynamic())
        m_eff = pa.mass;
    else
        m_eff = (pa.mass * pb.mass) / (pa.mass + pb.mass);

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
    //
    // 注意（缺陷修复）：
    //   标准线性弹簧-阻尼模型（Cundall-Strack）在回弹阶段允许阻尼项产生
    //   小幅"拉力"（c_n·v_n_rel 使 F_n 为负），这正是
    //       c_n = 2·beta·sqrt(m_eff·k_n),  beta = -ln(e)/sqrt(pi² + ln²e)
    //   这一恢复系数标定公式成立的前提。旧实现中的
    //       if (Fn_mag < 0.0) Fn_mag = 0.0;   // 不允许拉力
    //   会在回弹阶段抹掉这部分耗能，使实测恢复系数系统性高于标称值
    //   （例如 e=0.5 实测约 0.55）。接触的建立与解除由几何重叠判据
    //   （dist < sum_r）负责，无需用力钳位来防止颗粒粘连。
    double Fn_mag = k_n * contact.overlap_n + c_n * v_n_rel;

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
