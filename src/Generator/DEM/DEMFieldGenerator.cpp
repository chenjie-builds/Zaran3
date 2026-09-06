#include "DEMFieldGenerator.h"
#include "ReadDEMParticle.h"
#include "ReadDEMBond.h"
#include "GlobalData.h"
#include "Log.h"
#include "ZaranError.h"
#include "File.h"
#include "Rand.h"

namespace zaran
{

std::string DEMFieldGenerator::GetParticleFilePath() const
{
    std::string work_dir  = GlobalData::GetString("work_dir");
    std::string para_file = GlobalData::IsExist("dem.particle_file")
                              ? GlobalData::GetString("dem.particle_file")
                              : "particles.csv";
    return work_dir + "/" + para_file;
}

std::string DEMFieldGenerator::GetBondFilePath() const
{
    if (!GlobalData::IsExist("dem.bond_file")) return {};
    return GlobalData::GetString("work_dir") + "/" + GlobalData::GetString("dem.bond_file");
}

shared_ptr<FieldManager> DEMFieldGenerator::Create()
{
    shared_ptr<FieldManager> field_manager = make_shared<FieldManager>();
    shared_ptr<DEMField> dem_field = make_shared<DEMField>();
    dem_field->SetIdx(0);
    dem_field->Allocate();

    std::string full_path = GetParticleFilePath();

    if (IsFileExist(full_path))
    {
        LoadParticlesFromFile(dem_field);
    }
    else
    {
        Log::info("DEMFieldGenerator: particle file '{}' not found, generating from box.", full_path);
        GenerateParticlesFromBox(dem_field);
    }

    if (GlobalData::IsExist("dem.bond_file"))
        LoadBondsFromFile(dem_field);

    AddBoxWalls(dem_field);

    field_manager->AddField(dem_field, nullptr);
    Log::info("DEMFieldGenerator: created DEMField with {} particles, {} bonds, {} walls",
              dem_field->GetDEMData()->GetParticleNum(),
              dem_field->GetDEMData()->GetBondNum(),
              dem_field->GetDEMData()->GetWallNum());
    return field_manager;
}

void DEMFieldGenerator::LoadBondsFromFile(const shared_ptr<DEMField>& field) const
{
    const std::string full_path = GetBondFilePath();
    if (!IsFileExist(full_path))
        throw ZaranError("DEM bond file not found: " + full_path);

    ReadDEMBond reader;
    std::vector<DEMBond> bonds;
    reader.ReadCSV(full_path, field->GetDEMData()->GetParticles(), bonds);
    for (const auto& bond : bonds)
        field->GetDEMData()->AddBond(bond);
}

void DEMFieldGenerator::LoadParticlesFromFile(const shared_ptr<DEMField>& field) const
{
    std::string full_path = GetParticleFilePath();

    ReadDEMParticle reader;
    std::vector<DEMParticle> particles;
    reader.ReadCSV(full_path, particles);

    auto dem_data = field->GetDEMData();
    for (auto& p : particles)
        dem_data->AddParticle(p);
}

void DEMFieldGenerator::GenerateParticlesFromBox(const shared_ptr<DEMField>& field) const
{
    // 读取包围盒参数（复用 zaran.box.* 或 dem.box.*）
    auto getD = [](const std::string& key, double def) -> double {
        return GlobalData::IsExist(key) ? GlobalData::GetDouble(key) : def;
    };
    double x_min = getD("dem.box.x_min", -1.0);
    double x_max = getD("dem.box.x_max",  1.0);
    double y_min = getD("dem.box.y_min", -1.0);
    double y_max = getD("dem.box.y_max",  1.0);
    double z_min = getD("dem.box.z_min", -1.0);
    double z_max = getD("dem.box.z_max",  1.0);

    double r    = getD("dem.particle_radius", 0.05);
    int    num  = GlobalData::IsExist("dem.particle_num") ? GlobalData::GetInt("dem.particle_num") : 100;

    if (r <= 0.0)
        throw ZaranError("dem.particle_radius must be positive");
    if (num < 0)
        throw ZaranError("dem.particle_num must be non-negative");
    if (x_max - x_min < 2.0 * r || y_max - y_min < 2.0 * r || z_max - z_min < 2.0 * r)
        throw ZaranError("DEM box must be at least one particle diameter in every direction");

    auto dem_data = field->GetDEMData();
    for (int i = 0; i < num; ++i)
    {
        DEMParticle p;
        p.id     = static_cast<index_type>(i);
        p.group  = 0;
        p.radius = r;
        // 质量在 DEMSolver::InitField 中由 dem.density 计算。
        p.mass   = 0.0;
        p.inertia = 0.4;
        bool placed = false;
        for (int attempt = 0; attempt < 10000 && !placed; ++attempt)
        {
            p.pos.x() = x_min + r + (x_max - x_min - 2.0 * r) * RandDouble(0.0, 1.0);
            p.pos.y() = y_min + r + (y_max - y_min - 2.0 * r) * RandDouble(0.0, 1.0);
            p.pos.z() = z_min + r + (z_max - z_min - 2.0 * r) * RandDouble(0.0, 1.0);
            placed = true;
            for (const auto& existing : dem_data->GetParticles())
            {
                if ((p.pos - existing.pos).squaredNorm()
                    < (p.radius + existing.radius) * (p.radius + existing.radius))
                {
                    placed = false;
                    break;
                }
            }
        }
        if (!placed)
            throw ZaranError("Unable to place non-overlapping DEM particles in the configured box");
        dem_data->AddParticle(p);
    }
    Log::info("DEMFieldGenerator: generated {} particles in box [{},{},{} → {},{},{}]",
              num, x_min, y_min, z_min, x_max, y_max, z_max);
}

void DEMFieldGenerator::AddBoxWalls(const shared_ptr<DEMField>& field) const
{
    auto getD = [](const std::string& key, double def) -> double {
        return GlobalData::IsExist(key) ? GlobalData::GetDouble(key) : def;
    };
    double x_min = getD("dem.box.x_min", -1.0);
    double x_max = getD("dem.box.x_max",  1.0);
    double y_min = getD("dem.box.y_min", -1.0);
    double y_max = getD("dem.box.y_max",  1.0);
    double z_min = getD("dem.box.z_min", -1.0);
    double z_max = getD("dem.box.z_max",  1.0);

    auto dem_data = field->GetDEMData();
    index_type wid = 0;

    auto make_wall = [&](Eigen::Vector3d n, Eigen::Vector3d pt) {
        DEMWall w;
        w.id     = wid++;
        w.normal = n.normalized();
        w.point  = pt;
        dem_data->AddWall(w);
    };

    // 六面墙：法向指向域内
    make_wall({ 1.0, 0.0, 0.0}, {x_min, 0.0, 0.0}); // x-min
    make_wall({-1.0, 0.0, 0.0}, {x_max, 0.0, 0.0}); // x-max
    make_wall({ 0.0, 1.0, 0.0}, {0.0, y_min, 0.0}); // y-min
    make_wall({ 0.0,-1.0, 0.0}, {0.0, y_max, 0.0}); // y-max
    make_wall({ 0.0, 0.0, 1.0}, {0.0, 0.0, z_min}); // z-min
    make_wall({ 0.0, 0.0,-1.0}, {0.0, 0.0, z_max}); // z-max
}

} // namespace zaran
