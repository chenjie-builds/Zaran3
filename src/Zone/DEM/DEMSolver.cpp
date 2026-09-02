#include "DEMSolver.h"
#include "LinearSpringDashpot.h"
#include "HertzMindlin.h"
#include "KDTree.h"
#include "Log.h"
#include "ZaranError.h"
#include "CommonPara.h"
#include <fstream>

namespace zaran
{

DEMSolver::DEMSolver(index_type index,
                     const std::string& name,
                     shared_ptr<DEMSolverParam> para,
                     shared_ptr<DEMFieldData>   dem_data)
    : FieldSolver(index, name, para, nullptr, nullptr), // DEM 无背景网格，无传统 DataManager
      m_dem_data(std::move(dem_data))
{
}

DEMSolverParam* DEMSolver::GetDEMParam() const
{
    return static_cast<DEMSolverParam*>(GetPara());
}

void DEMSolver::Init()
{
    InitSolver();
    InitField();
}

void DEMSolver::InitSolver()
{
    auto* para = GetDEMParam();
    para->Init();

    const std::string& model_name = para->GetContactModel();
    if (model_name == "LinearSpringDashpot" || model_name == "Linear")
    {
        m_contact_model = make_unique<LinearSpringDashpot>();
        Log::info("DEMSolver: contact model = LinearSpringDashpot");
    }
    else if (model_name == "HertzMindlin" || model_name == "Hertz")
    {
        m_contact_model = make_unique<HertzMindlin>();
        Log::info("DEMSolver: contact model = HertzMindlin");
    }
    else
    {
        Log::warn("DEMSolver: unknown contact model '{}', fallback to HertzMindlin", model_name);
        m_contact_model = make_unique<HertzMindlin>();
    }
}

void DEMSolver::InitField()
{
    auto* para = GetDEMParam();
    // 将默认材料参数应用到每个粒子（若粒子未显式指定则用全局值）
    for (auto& p : m_dem_data->GetParticles())
    {
        if (p.mass <= 0.0)
        {
            // 默认球体：m = 4/3 π r³ ρ，这里 ρ 由 young_modulus 量级代替（仅占位）
            // 实际密度应由粒子文件给出；这里只是防御性补零
            double rho = 2500.0; // kg/m³
            p.mass    = (4.0 / 3.0) * PI * std::pow(p.radius, 3) * rho;
            p.inertia = 0.4; // I = 2/5 m r²
        }
        p.young_modulus     = para->GetYoungModulus();
        p.poisson_ratio     = para->GetPoissonRatio();
        p.friction_coeff    = para->GetFrictionCoeff();
        p.restitution_coeff = para->GetRestitutionCoeff();
    }
    Log::info("DEMSolver InitField: {} particles", m_dem_data->GetParticleNum());
}

void DEMSolver::Preprocess()
{
    // 当前版本无前处理工作（扩展时可加：MPI 通信等）
}

void DEMSolver::Postprocess()
{
    // 当前版本无后处理工作
}

void DEMSolver::Solve()
{
    double dt = GetDEMParam()->GetTimeStep();

    ZeroForce();
    ContactDetection();
    CalcContactForce();
    CalcWallForce();
    CalcGravity();
    Integrate();
}

void DEMSolver::ZeroForce()
{
    for (auto& p : m_dem_data->GetParticles())
    {
        p.force.setZero();
        p.torque.setZero();
    }
}

void DEMSolver::ContactDetection()
{
    auto& particles = m_dem_data->GetParticles();
    const index_type N = particles.size();

    m_dem_data->ClearContacts();

    if (N == 0) return;

    // 构建 KDTree（仅含 active 粒子）
    point_vec pts;
    index_vec active_indices;
    pts.reserve(N);
    active_indices.reserve(N);
    for (index_type i = 0; i < N; ++i)
    {
        if (!particles[i].active) continue;
        pts.push_back({ particles[i].pos.x(), particles[i].pos.y(), particles[i].pos.z() });
        active_indices.push_back(i);
    }

    if (active_indices.size() < 2) return;

    KDTree kd(pts);

    // 对每个粒子 i，查询半径 = r_i + r_max_neighbor 内的邻居
    double r_max = 0.0;
    for (auto i : active_indices)
        r_max = std::max(r_max, particles[i].radius);

    for (index_type ki = 0; ki < active_indices.size(); ++ki)
    {
        index_type i = active_indices[ki];
        const DEMParticle& pa = particles[i];
        coord_vec pt_i = { pa.pos.x(), pa.pos.y(), pa.pos.z() };

        double search_r = pa.radius + r_max;
        auto neighbors = kd.NeighborhoodIndices(pt_i, search_r);

        for (index_type kj : neighbors)
        {
            index_type j = active_indices[kj];
            if (j <= i) continue; // 每对只算一次

            const DEMParticle& pb = particles[j];
            Eigen::Vector3d d = pb.pos - pa.pos;
            double dist = d.norm();
            double sum_r = pa.radius + pb.radius;
            if (dist >= sum_r || dist < 1.0e-15) continue;

            DEMContact c;
            c.type      = ContactType::ParticleParticle;
            c.idx_a     = i;
            c.idx_b     = j;
            c.overlap_n = sum_r - dist;
            c.normal    = d / dist; // n 从 A 指向 B
            c.contact_point = pa.pos + (pa.radius - 0.5 * c.overlap_n) * c.normal;
            m_dem_data->AddContact(c);
        }
    }
}

void DEMSolver::CalcContactForce()
{
    auto& particles = m_dem_data->GetParticles();
    double dt = GetDEMParam()->GetTimeStep();

    for (auto& c : m_dem_data->GetContacts())
    {
        if (c.type != ContactType::ParticleParticle) continue;
        DEMParticle& pa = particles[c.idx_a];
        DEMParticle& pb = particles[c.idx_b];

        m_contact_model->CalcNormalForce(pa, pb, c, dt);
        m_contact_model->CalcTangentialForce(pa, pb, c, dt);

        Eigen::Vector3d F = c.force_n + c.force_t;

        // 作用-反作用
        if (pa.active) { pa.force += F; pa.torque += c.contact_point.cross(c.force_t) - pa.pos.cross(c.force_t); }
        if (pb.active) { pb.force -= F; pb.torque -= c.contact_point.cross(c.force_t) - pb.pos.cross(c.force_t); }
    }
}

void DEMSolver::CalcWallForce()
{
    auto& particles = m_dem_data->GetParticles();
    auto& walls     = m_dem_data->GetWalls();
    double dt = GetDEMParam()->GetTimeStep();

    for (auto& wall : walls)
    {
        for (auto& pa : particles)
        {
            if (!pa.active) continue;
            double d = wall.SignedDist(pa.pos);
            double overlap = pa.radius - d;
            if (overlap <= 0.0) continue;

            // 构造虚拟墙粒子用于接触力计算
            DEMParticle pb_wall;
            pb_wall.radius           = 1.0e10; // 极大半径 → 平面
            pb_wall.mass             = 1.0e30;
            pb_wall.vel.setZero();
            pb_wall.omega.setZero();
            pb_wall.young_modulus    = wall.young_modulus;
            pb_wall.poisson_ratio    = wall.poisson_ratio;
            pb_wall.friction_coeff   = wall.friction_coeff;
            pb_wall.restitution_coeff= wall.restitution_coeff;
            // 墙粒子质心沿法向偏移极大半径
            pb_wall.pos = pa.pos - (d - pa.radius + 1.0e10) * wall.normal;

            DEMContact c;
            c.type      = ContactType::ParticleWall;
            c.idx_a     = 0;
            c.idx_b     = wall.id;
            c.overlap_n = overlap;
            c.normal    = wall.normal; // 法向从 A 指向墙
            c.contact_point = pa.pos - d * wall.normal;

            m_contact_model->CalcNormalForce(pa, pb_wall, c, dt);
            m_contact_model->CalcTangentialForce(pa, pb_wall, c, dt);

            pa.force  += c.force_n + c.force_t;
            pa.torque += c.contact_point.cross(c.force_t) - pa.pos.cross(c.force_t);
        }
    }
}

void DEMSolver::CalcGravity()
{
    const Eigen::Vector3d& g = GetDEMParam()->GetGravity();
    for (auto& p : m_dem_data->GetParticles())
    {
        if (p.active)
            p.force += p.mass * g;
    }
}

void DEMSolver::Integrate()
{
    double dt = GetDEMParam()->GetTimeStep();
    for (auto& p : m_dem_data->GetParticles())
    {
        if (!p.active) continue;
        // Velocity-Verlet（半步速度已在上一步末尾完成；此处简化为前向 Euler，
        //  完整 Verlet 需分两半步—可在后续迭代中完善）
        Eigen::Vector3d acc   = p.force  / p.mass;
        Eigen::Vector3d alpha = p.torque / (p.inertia * p.mass * p.radius * p.radius);

        p.vel   += acc   * dt;
        p.omega += alpha * dt;
        p.pos   += p.vel * dt;
    }
}

void DEMSolver::BackupField(std::string& back_folder)
{
    BackupField(static_cast<const std::string&>(back_folder));
}

void DEMSolver::BackupField(const std::string& back_folder) const
{
    std::string fname = back_folder + "/particles.csv";
    std::ofstream fout(fname);
    if (!fout.is_open())
    {
        Log::warn("DEMSolver::BackupField: cannot open {}", fname);
        return;
    }
    fout << "id,group,radius,mass,px,py,pz,vx,vy,vz,ox,oy,oz\n";
    for (const auto& p : m_dem_data->GetParticles())
    {
        fout << p.id    << ","
             << p.group << ","
             << p.radius << ","
             << p.mass  << ","
             << p.pos.x()   << "," << p.pos.y()   << "," << p.pos.z()   << ","
             << p.vel.x()   << "," << p.vel.y()   << "," << p.vel.z()   << ","
             << p.omega.x() << "," << p.omega.y() << "," << p.omega.z() << "\n";
    }
    Log::info("DEMSolver::BackupField: {} particles written to {}", m_dem_data->GetParticleNum(), fname);
}

} // namespace zaran
