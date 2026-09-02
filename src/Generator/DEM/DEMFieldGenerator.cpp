#include "DEMFieldGenerator.h"
#include "ReadDEMParticle.h"
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

    AddBoxWalls(dem_field);

    field_manager->AddField(dem_field, nullptr);
    Log::info("DEMFieldGenerator: created DEMField with {} particles, {} walls",
              dem_field->GetDEMData()->GetParticleNum(),
              dem_field->GetDEMData()->GetWallNum());
    return field_manager;
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

    auto dem_data = field->GetDEMData();
    const double rho = 2500.0; // kg/m³

    for (int i = 0; i < num; ++i)
    {
        DEMParticle p;
        p.id     = static_cast<index_type>(i);
        p.group  = 0;
        p.radius = r;
        p.mass   = (4.0 / 3.0) * 3.141592653589793 * r * r * r * rho;
        p.inertia = 0.4;
        p.pos.x() = x_min + (x_max - x_min) * RandDouble(0.0, 1.0);
        p.pos.y() = y_min + (y_max - y_min) * RandDouble(0.0, 1.0);
        p.pos.z() = z_min + (z_max - z_min) * RandDouble(0.0, 1.0);
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
