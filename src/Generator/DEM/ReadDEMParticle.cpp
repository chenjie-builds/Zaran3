#include "ReadDEMParticle.h"
#include "Log.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <set>

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

    std::set<index_type> ids;
    size_t line_number = 1;
    while (std::getline(fin, line))
    {
        ++line_number;
        if (line.empty() || line[0] == '#') continue;

        // 替换逗号为空格后解析
        for (char& c : line) if (c == ',') c = ' ';
        std::istringstream ss(line);

        DEMParticle p;
        ss >> p.id >> p.group >> p.radius >> p.mass;
        if (ss.fail())
        {
            // 最少需要 id group radius mass
            throw std::runtime_error("ReadDEMParticle: malformed required fields at line "
                                     + std::to_string(line_number));
        }
        double px = 0, py = 0, pz = 0;
        double vx = 0, vy = 0, vz = 0;
        double ox = 0, oy = 0, oz = 0;
        ss >> px >> py >> pz;
        ss >> vx >> vy >> vz;
        ss >> ox >> oy >> oz;

        if (ss.fail())
            throw std::runtime_error("ReadDEMParticle: expected 13 fields at line "
                                     + std::to_string(line_number));

        p.pos   = Eigen::Vector3d(px, py, pz);
        p.vel   = Eigen::Vector3d(vx, vy, vz);
        p.omega = Eigen::Vector3d(ox, oy, oz);

        std::vector<double> optional;
        double optional_value = 0.0;
        while (ss >> optional_value)
            optional.push_back(optional_value);
        if (!ss.eof() || (optional.size() != 0 && optional.size() != 1 && optional.size() != 5))
            throw std::runtime_error("ReadDEMParticle: expected 13, 14, or 18 fields at line "
                                     + std::to_string(line_number));
        if (!optional.empty())
        {
            if (optional[0] != 0.0 && optional[0] != 1.0 && optional[0] != 2.0)
                throw std::runtime_error("ReadDEMParticle: motion type must be 0, 1, or 2 at line "
                                         + std::to_string(line_number));
            p.active = optional[0] != 0.0;
            p.kinematic = optional[0] == 2.0;
        }
        if (optional.size() == 5)
        {
            p.young_modulus = optional[1];
            p.poisson_ratio = optional[2];
            p.friction_coeff = optional[3];
            p.restitution_coeff = optional[4];
            p.material_from_file = true;
        }

        if (p.radius <= 0.0 || p.mass < 0.0 || !p.pos.allFinite()
            || !p.vel.allFinite() || !p.omega.allFinite())
            throw std::runtime_error("ReadDEMParticle: invalid particle values at line "
                                     + std::to_string(line_number));
        if (!ids.insert(p.id).second)
            throw std::runtime_error("ReadDEMParticle: duplicate particle id at line "
                                     + std::to_string(line_number));
        if (p.material_from_file &&
            (p.young_modulus <= 0.0 || p.poisson_ratio <= -1.0 || p.poisson_ratio >= 0.5
             || p.friction_coeff < 0.0 || p.restitution_coeff <= 0.0
             || p.restitution_coeff > 1.0))
            throw std::runtime_error("ReadDEMParticle: invalid material values at line "
                                     + std::to_string(line_number));
        p.inertia = 0.4; // 2/5 (球体)

        particles.push_back(p);
    }
    Log::info("ReadDEMParticle: loaded {} particles from '{}'", particles.size(), filename);
}

} // namespace zaran
