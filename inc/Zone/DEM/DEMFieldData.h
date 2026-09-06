/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DEMFieldData.h
 * \brief DEMFieldData class, storing the particle array for a DEM field.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "BasicType.h"
#include "DEMParticle.h"
#include "DEMWall.h"
#include "DEMContact.h"
#include "DEMBond.h"
namespace zaran
{
    /// @brief DEM 场数据：粒子、墙面、当前接触对和永久弹性键合
    class DEMFieldData
    {
    public:
        DEMFieldData() = default;
        ~DEMFieldData() = default;

        // --- 粒子管理 ---
        void AddParticle(const DEMParticle& p) { m_particles.push_back(p); }
        void RemoveParticle(index_type idx)
        {
            if (idx < m_particles.size())
                m_particles.erase(m_particles.begin() + static_cast<long>(idx));
        }
        DEMParticle& GetParticle(index_type idx) { return m_particles[idx]; }
        const DEMParticle& GetParticle(index_type idx) const { return m_particles[idx]; }
        index_type GetParticleNum() const { return m_particles.size(); }
        dynamic_array<DEMParticle>& GetParticles() { return m_particles; }
        const dynamic_array<DEMParticle>& GetParticles() const { return m_particles; }

        // --- 墙面管理 ---
        void AddWall(const DEMWall& w) { m_walls.push_back(w); }
        DEMWall& GetWall(index_type idx) { return m_walls[idx]; }
        const DEMWall& GetWall(index_type idx) const { return m_walls[idx]; }
        index_type GetWallNum() const { return m_walls.size(); }
        dynamic_array<DEMWall>& GetWalls() { return m_walls; }

        // --- 接触对管理 ---
        void ClearContacts() { m_contacts.clear(); }
        void AddContact(const DEMContact& c) { m_contacts.push_back(c); }
        dynamic_array<DEMContact>& GetContacts() { return m_contacts; }
        const dynamic_array<DEMContact>& GetContacts() const { return m_contacts; }

        // --- 永久键合管理 ---
        void AddBond(const DEMBond& bond) { m_bonds.push_back(bond); }
        index_type GetBondNum() const { return m_bonds.size(); }
        dynamic_array<DEMBond>& GetBonds() { return m_bonds; }
        const dynamic_array<DEMBond>& GetBonds() const { return m_bonds; }

    private:
        dynamic_array<DEMParticle> m_particles;
        dynamic_array<DEMWall>     m_walls;
        dynamic_array<DEMContact>  m_contacts;
        dynamic_array<DEMBond>     m_bonds;
    };
} // namespace zaran
