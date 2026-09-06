#include "ReadDEMBond.h"
#include "Log.h"
#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace zaran
{
void ReadDEMBond::ReadCSV(const std::string& filename,
                          const std::vector<DEMParticle>& particles,
                          std::vector<DEMBond>& bonds) const
{
    std::ifstream fin(filename);
    if (!fin.is_open())
        throw std::runtime_error("ReadDEMBond: cannot open " + filename);

    std::map<index_type, index_type> particle_indices;
    for (index_type i = 0; i < particles.size(); ++i)
        particle_indices.emplace(particles[i].id, i);

    std::set<index_type> bond_ids;
    std::set<std::pair<index_type, index_type>> endpoint_pairs;
    std::string line;
    std::getline(fin, line); // header
    size_t line_number = 1;
    while (std::getline(fin, line))
    {
        ++line_number;
        if (line.empty() || line[0] == '#') continue;
        for (char& c : line) if (c == ',') c = ' ';

        index_type particle_a_id = 0, particle_b_id = 0;
        int active = 1;
        DEMBond bond;
        std::istringstream ss(line);
        ss >> bond.id >> particle_a_id >> particle_b_id >> bond.rest_length
           >> bond.normal_stiffness >> bond.tangential_stiffness
           >> bond.normal_damping >> bond.tangential_damping >> active;
        if (ss.fail())
            throw std::runtime_error("ReadDEMBond: expected 9 fields at line "
                                     + std::to_string(line_number));
        double extra = 0.0;
        if (ss >> extra)
            throw std::runtime_error("ReadDEMBond: too many fields at line "
                                     + std::to_string(line_number));
        if (active != 0 && active != 1)
            throw std::runtime_error("ReadDEMBond: active must be 0 or 1 at line "
                                     + std::to_string(line_number));
        if (particle_a_id == particle_b_id)
            throw std::runtime_error("ReadDEMBond: a bond cannot connect a particle to itself");
        const auto a = particle_indices.find(particle_a_id);
        const auto b = particle_indices.find(particle_b_id);
        if (a == particle_indices.end() || b == particle_indices.end())
            throw std::runtime_error("ReadDEMBond: unknown particle id at line "
                                     + std::to_string(line_number));
        if (bond.rest_length < 0.0 || bond.normal_stiffness < 0.0
            || bond.tangential_stiffness < 0.0 || bond.normal_damping < 0.0
            || bond.tangential_damping < 0.0)
            throw std::runtime_error("ReadDEMBond: bond parameters must be non-negative");
        if (!bond_ids.insert(bond.id).second)
            throw std::runtime_error("ReadDEMBond: duplicate bond id");

        bond.idx_a = a->second;
        bond.idx_b = b->second;
        const auto pair = std::minmax(bond.idx_a, bond.idx_b);
        if (!endpoint_pairs.emplace(pair.first, pair.second).second)
            throw std::runtime_error("ReadDEMBond: duplicate particle pair");
        if (bond.rest_length == 0.0)
            bond.rest_length = (particles[bond.idx_b].pos - particles[bond.idx_a].pos).norm();
        if (bond.rest_length <= 0.0)
            throw std::runtime_error("ReadDEMBond: rest length must be positive");
        bond.active = active != 0;
        bonds.push_back(bond);
    }
    Log::info("ReadDEMBond: loaded {} bonds from '{}'", bonds.size(), filename);
}
} // namespace zaran
