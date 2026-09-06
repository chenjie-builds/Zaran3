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
    return static_cast<DEMSolverParam*>(const_cast<DEMSolver*>(this)->GetPara());
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
        throw ZaranError("Unknown DEM contact model: " + model_name);
    }
}

void DEMSolver::InitField()
{
    auto* para = GetDEMParam();
    // 将默认材料参数应用到每个粒子（若粒子未显式指定则用全局值）
    for (auto& p : m_dem_data->GetParticles())
    {
        if (p.radius <= 0.0 || p.inertia <= 0.0)
            throw ZaranError("DEM particle radius and inertia factor must be positive");
        if (p.mass <= 0.0)
        {
            // 默认球体：m = 4/3 π r³ ρ
            double rho = para->GetDensity();
            p.mass    = (4.0 / 3.0) * PI * std::pow(p.radius, 3) * rho;
            p.inertia = 0.4; // I = 2/5 m r²
        }
        if (!p.pos.allFinite() || !p.vel.allFinite() || !p.omega.allFinite())
            throw ZaranError("DEM particle state contains a non-finite value");
        // 仅在粒子未从文件显式设置材料参数时，才使用全局参数覆盖
        if (!p.material_from_file)
        {
            p.young_modulus     = para->GetYoungModulus();
            p.poisson_ratio     = para->GetPoissonRatio();
            p.friction_coeff    = para->GetFrictionCoeff();
            p.restitution_coeff = para->GetRestitutionCoeff();
        }
    }
    for (const auto& bond : m_dem_data->GetBonds())
    {
        if (bond.idx_a >= m_dem_data->GetParticleNum()
            || bond.idx_b >= m_dem_data->GetParticleNum()
            || bond.idx_a == bond.idx_b || bond.rest_length <= 0.0
            || bond.normal_stiffness < 0.0 || bond.tangential_stiffness < 0.0
            || bond.normal_damping < 0.0 || bond.tangential_damping < 0.0)
            throw ZaranError("DEM bond data are invalid");
    }
    Log::info("DEMSolver InitField: {} particles, {} bonds",
              m_dem_data->GetParticleNum(), m_dem_data->GetBondNum());
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
    m_contacts_this_step.clear();
    ZeroForce();
    CalcBondForce();
    ContactDetection();
    CalcContactForce();
    CalcWallForce();
    PruneContactHistory();
    CalcGravity();
    Integrate();
}

void DEMSolver::CalcBondForce()
{
    auto& particles = m_dem_data->GetParticles();
    const double dt = GetDEMParam()->GetTimeStep();
    for (auto& bond : m_dem_data->GetBonds())
    {
        if (!bond.active) continue;
        DEMParticle& pa = particles[bond.idx_a];
        DEMParticle& pb = particles[bond.idx_b];
        const Eigen::Vector3d separation = pb.pos - pa.pos;
        const double length = separation.norm();
        if (length <= 1.0e-15)
            throw ZaranError("A DEM bond has coincident endpoints");

        const Eigen::Vector3d normal = separation / length;
        const Eigen::Vector3d arm_a = 0.5 * length * normal;
        const Eigen::Vector3d arm_b = -arm_a;
        Eigen::Vector3d relative_velocity = pa.vel + pa.omega.cross(arm_a)
                                          - pb.vel - pb.omega.cross(arm_b);
        const double velocity_n = relative_velocity.dot(normal);
        const Eigen::Vector3d velocity_t = relative_velocity - velocity_n * normal;

        bond.delta_t += velocity_t * dt;
        bond.delta_t -= bond.delta_t.dot(normal) * normal;
        bond.extension = length - bond.rest_length;

        const Eigen::Vector3d force_n = (bond.normal_stiffness * bond.extension
                                          - bond.normal_damping * velocity_n) * normal;
        const Eigen::Vector3d force_t = -bond.tangential_stiffness * bond.delta_t
                                        - bond.tangential_damping * velocity_t;
        bond.force_a = force_n + force_t;

        pa.force += bond.force_a;
        pa.torque += arm_a.cross(force_t);
        pb.force -= bond.force_a;
        pb.torque += arm_b.cross(-force_t);
    }
}

DEMSolver::ContactKey DEMSolver::MakeContactKey(const DEMContact& contact) const
{
    return std::make_tuple(static_cast<int>(contact.type), contact.idx_a, contact.idx_b);
}

void DEMSolver::RestoreContactHistory(DEMContact& contact) const
{
    const auto it = m_tangential_history.find(MakeContactKey(contact));
    if (it != m_tangential_history.end())
        contact.delta_t = it->second;
}

void DEMSolver::SaveContactHistory(const DEMContact& contact)
{
    const ContactKey key = MakeContactKey(contact);
    m_tangential_history[key] = contact.delta_t;
    m_contacts_this_step.insert(key);
}

void DEMSolver::PruneContactHistory()
{
    for (auto it = m_tangential_history.begin(); it != m_tangential_history.end();)
    {
        if (m_contacts_this_step.find(it->first) == m_contacts_this_step.end())
            it = m_tangential_history.erase(it);
        else
            ++it;
    }
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

    // 固定粒子也进入邻域搜索；它们不积分，但仍可向活动粒子施加接触力。
    point_vec pts;
    pts.reserve(N);
    for (index_type i = 0; i < N; ++i)
    {
        pts.push_back({ particles[i].pos.x(), particles[i].pos.y(), particles[i].pos.z() });
    }

    if (N < 2) return;

    KDTree kd(pts);

    // 对每个粒子 i，查询半径 = r_i + r_max_neighbor 内的邻居
    double r_max = 0.0;
    for (const auto& particle : particles)
        r_max = std::max(r_max, particle.radius);

    for (index_type i = 0; i < N; ++i)
    {
        const DEMParticle& pa = particles[i];
        coord_vec pt_i = { pa.pos.x(), pa.pos.y(), pa.pos.z() };

        double search_r = pa.radius + r_max;
        auto neighbors = kd.NeighborhoodIndices(pt_i, search_r);

        for (index_type kj : neighbors)
        {
            index_type j = kj;
            if (j <= i) continue; // 每对只算一次

            const DEMParticle& pb = particles[j];
            if (!pa.IsDynamic() && !pb.IsDynamic()) continue;
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
            RestoreContactHistory(c);
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
        SaveContactHistory(c);

        Eigen::Vector3d F = c.force_n + c.force_t;

        // 作用-反作用
        pa.force += F;
        pa.torque += c.contact_point.cross(c.force_t) - pa.pos.cross(c.force_t);
        pb.force -= F;
        pb.torque -= c.contact_point.cross(c.force_t) - pb.pos.cross(c.force_t);
    }
}

void DEMSolver::CalcWallForce()
{
    auto& particles = m_dem_data->GetParticles();
    auto& walls     = m_dem_data->GetWalls();
    double dt = GetDEMParam()->GetTimeStep();

    for (auto& wall : walls)
    {
        for (index_type pi = 0; pi < particles.size(); ++pi)
        {
            DEMParticle& pa = particles[pi];
            if (!pa.IsDynamic()) continue;
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
            c.idx_a     = pi;
            c.idx_b     = wall.id;
            c.overlap_n = overlap;
            c.normal    = -wall.normal; // 法向从 A 指向墙
            c.contact_point = pa.pos - pa.radius * wall.normal;
            RestoreContactHistory(c);

            m_contact_model->CalcNormalForce(pa, pb_wall, c, dt);
            m_contact_model->CalcTangentialForce(pa, pb_wall, c, dt);
            SaveContactHistory(c);

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
        if (p.IsDynamic())
            p.force += p.mass * g;
    }
}

void DEMSolver::Integrate()
{
    double dt = GetDEMParam()->GetTimeStep();
    for (auto& p : m_dem_data->GetParticles())
    {
        if (!p.active) continue;
        if (p.kinematic)
        {
            p.rotation += p.omega * dt;
            p.pos += p.vel * dt;
            continue;
        }
        // 半隐式 Euler：先更新速度，再用新速度更新位置。
        Eigen::Vector3d acc   = p.force  / p.mass;
        Eigen::Vector3d alpha = p.torque / (p.inertia * p.mass * p.radius * p.radius);

        p.vel   += acc   * dt;
        p.omega += alpha * dt;
        p.rotation += p.omega * dt;
        p.pos   += p.vel * dt;
        
    }
}

void DEMSolver::BackupField(std::string& back_folder)
{
    BackupField(static_cast<const std::string&>(back_folder));
}

void DEMSolver::BackupField(const std::string& back_folder) const
{
    std::string fname = back_folder + "/particles.dat";
    std::ofstream fout(fname);
    if (!fout.is_open())
    {
        Log::warn("DEMSolver::BackupField: cannot open {}", fname);
        return;
    }
    fout << "id,group,radius,mass,px,py,pz,vx,vy,vz,rx,ry,rz,ox,oy,oz,fx,fy,fz,motion_type\n";
    for (const auto& p : m_dem_data->GetParticles())
    {
        fout << p.id    << ","
             << p.group << ","
             << p.radius << ","
             << p.mass  << ","
             << p.pos.x()   << "," << p.pos.y()   << "," << p.pos.z()   << ","
             << p.vel.x()   << "," << p.vel.y()   << "," << p.vel.z()   << ","
             << p.rotation.x() << "," << p.rotation.y() << "," << p.rotation.z() << ","
             << p.omega.x() << "," << p.omega.y() << "," << p.omega.z() << ","
             << p.force.x() << "," << p.force.y() << "," << p.force.z() << ","
             << (p.kinematic ? 2 : (p.active ? 1 : 0)) << "\n";
    }
    fout.close();

    const std::string bond_name = back_folder + "/bonds.dat";
    std::ofstream bond_out(bond_name);
    if (!bond_out.is_open())
    {
        Log::warn("DEMSolver::BackupField: cannot open {}", bond_name);
        return;
    }
    bond_out << "id,particle_a_id,particle_b_id,rest_length,extension,fx,fy,fz,active\n";
    const auto& particles = m_dem_data->GetParticles();
    for (const auto& bond : m_dem_data->GetBonds())
    {
        bond_out << bond.id << "," << particles[bond.idx_a].id << ","
                 << particles[bond.idx_b].id << "," << bond.rest_length << ","
                 << bond.extension << "," << bond.force_a.x() << ","
                 << bond.force_a.y() << "," << bond.force_a.z() << ","
                 << (bond.active ? 1 : 0) << "\n";
    }
    Log::info("DEMSolver::BackupField: {} particles and {} bonds written to {}",
              m_dem_data->GetParticleNum(), m_dem_data->GetBondNum(), back_folder);
}

} // namespace zaran
