/**
 * Zaran - A Totally Automatic CFD Software
 * \file ReadDEMBond.h
 * \brief Read inert elastic bonds from CSV.
 */
#pragma once
#include "DEMBond.h"
#include "DEMParticle.h"
#include <string>
#include <vector>

namespace zaran
{
    class ReadDEMBond
    {
    public:
        /// CSV: id,particle_a_id,particle_b_id,rest_length,kn,kt,cn,ct,active
        void ReadCSV(const std::string& filename,
                     const std::vector<DEMParticle>& particles,
                     std::vector<DEMBond>& bonds) const;
    };
} // namespace zaran
