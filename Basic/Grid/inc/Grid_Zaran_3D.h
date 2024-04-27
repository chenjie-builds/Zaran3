//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Grid_Zaran_3D.h 													||
//*	@brief	三维扎染网格类                                  						||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "Grid.h"
#include "ZaranBoundPatch.h"
#include "CellTopoInfoZaran.h"
namespace zaran
{
    class Grid_Zaran_3D :public Grid
    {
    public:
        Grid_Zaran_3D()
            : m_ni(0), m_nj(0), m_nk(0)
        {
            GridBase();
            SetDimension(Dimension::three);
            m_node_topo = new NodeTopo();
            m_face_topo = new FaceTopo();
            m_cell_topo = new CellTopoZaran();
            m_bound_map = new BoundaryMap();
            m_bound_patch = new ZaranBoundPatch();
        }
    public:
        void SetNi(int ni);
        void SetNj(int nj);
        void SetNk(int nk);
        int GetNi();
        int GetNj();
        int GetNk();
        void SetNodeNum(int ni, int nj, int nk);
        void GetNodeNum(int& ni, int& nj, int& nk);
        int GetNodeIndex(int i, int j, int k);
        void GetNodeIndex(int index, int& i, int& j, int& k);
        int GetCellIndex(int i, int j, int k);
        void GetCellIndex(int index, int& i, int& j, int& k);
        void SetBox(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
        void GetBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);
        CellTopoZaran* GetCellTopo()override;
        ZaranBoundPatch* GetBoundPatch();
        void GetRange(int& iStart, int& iEnd, int& jStart, int& jEnd, int& kStart, int& kEnd);
    private:
        int m_ni, m_nj, m_nk;//三个方向的节点个数
        ZaranBoundPatch* m_bound_patch;
        double m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, m_zmax;
        double m_dx, m_dy, m_dz;
    };
} // namespace zaran

