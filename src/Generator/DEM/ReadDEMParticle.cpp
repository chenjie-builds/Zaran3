#include "ReadDEMParticle.h"
#include "Log.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace zaran
{

void ReadDEMParticle::ReadCSV(const std::string& filename,
                               std::vector<DEMParticle>& particles) const
{
    std::ifstream fin(filename);
    if (!fin.is_open())
    {
        Log::error("ReadDEMParticle: cannot open file '{}'", filename);
        throw std::runtime_error("ReadDEMParticle: cannot open " + filename);
    }

    std::string line;
    // 跳过首行（表头）
    std::getline(fin, line);

    index_type auto_id = 0;
    while (std::getline(fin, line))
    {
        if (line.empty() || line[0] == '#') continue;

        // 替换逗号为空格后解析
        for (char& c : line) if (c == ',') c = ' ';
        std::istringstream ss(line);

        DEMParticle p;
        ss >> p.id >> p.group >> p.radius >> p.mass;
        if (ss.fail())
        {
            // 最少需要 id group radius mass
            Log::warn("ReadDEMParticle: skipping malformed line: {}", line);
            continue;
        }
        double px = 0, py = 0, pz = 0;
        double vx = 0, vy = 0, vz = 0;
        double ox = 0, oy = 0, oz = 0;
        ss >> px >> py >> pz;
        ss >> vx >> vy >> vz;
        ss >> ox >> oy >> oz;

        p.pos   = Eigen::Vector3d(px, py, pz);
        p.vel   = Eigen::Vector3d(vx, vy, vz);
        p.omega = Eigen::Vector3d(ox, oy, oz);

        if (p.mass <= 0.0)
        {
            const double rho = 2500.0; // 默认密度 (kg/m³)
            p.mass = (4.0 / 3.0) * 3.141592653589793 * p.radius * p.radius * p.radius * rho;
        }
        p.inertia = 0.4; // 2/5 (球体)

        particles.push_back(p);
        ++auto_id;
    }
    Log::info("ReadDEMParticle: loaded {} particles from '{}'", particles.size(), filename);
}

} // namespace zaran
