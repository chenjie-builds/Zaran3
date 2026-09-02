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
                                    DEMContact& contact, double /*dt*/)
{
    double E_star = 0.0, R_star = 0.0, G_star = 0.0;
    CalcEffectiveParams(pa, pb, E_star, R_star, G_star);

    double delta = contact.overlap_n;
    if (delta <= 0.0)
    {
        contact.force_n.setZero();
        return;
    }

    // Hertz 法向刚度：k_n = 2 E* sqrt(R* δ)
    double k_n = 2.0 * E_star * std::sqrt(R_star * delta);

    // 等效质量
    double m_eff = (pa.mass * pb.mass) / (pa.mass + pb.mass);

    // 临界阻尼比（基于恢复系数）
    double e = std::min(pa.restitution_coeff, pb.restitution_coeff);
    e = std::max(e, 1.0e-3);
    double ln_e = std::log(e);
    double beta = -ln_e / std::sqrt(PI * PI + ln_e * ln_e);
    double c_n = 2.0 * beta * std::sqrt(m_eff * k_n);

    // 法向相对速度
    Eigen::Vector3d rel_vel = pa.vel - pb.vel;
    double v_n_rel = rel_vel.dot(contact.normal);

    // F_n = (4/3 E* sqrt(R*) δ^{3/2} + c_n v_n_rel) 作用于 A 沿 -normal
    double Fn_hertz = (4.0 / 3.0) * E_star * std::sqrt(R_star) * std::pow(delta, 1.5);
    double Fn_damp  = c_n * v_n_rel;
    double Fn_mag = Fn_hertz + Fn_damp;
    if (Fn_mag < 0.0) Fn_mag = 0.0;

    contact.force_n = -Fn_mag * contact.normal;
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
        Ft = Ft.normalized() * Ft_max;
        contact.delta_t = -Ft / k_t;
    }

    contact.force_t = Ft;
}

} // namespace zaran
