#include "Grid_Zaran_3D.h"
namespace zaran
{
    void Grid_Zaran_3D::SetBoundPatch(const ZaranBoundPatch& bound_patch)
    {
        m_bound_patch = std::make_shared<ZaranBoundPatch>(bound_patch);
    }

    void Grid_Zaran_3D::SetBox(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
    {
        m_xmin = xmin;
        m_xmax = xmax;
        m_ymin = ymin;
        m_ymax = ymax;
        m_zmin = zmin;
        m_zmax = zmax;
    }

    void Grid_Zaran_3D::GetBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax)
    {
        xmin = m_xmin;
        xmax = m_xmax;
        ymin = m_ymin;
        ymax = m_ymax;
        zmin = m_zmin;
        zmax = m_zmax;
    }

    ZaranBoundPatch& Grid_Zaran_3D::GetBoundPatch()
    {
        return *m_bound_patch;
    }

void Grid_Zaran_3D::GetRange(int& iStart, int& iEnd, int& jStart, int& jEnd, int& kStart, int& kEnd)
    {
        iStart = 0;
        iEnd = m_ni - 1;
        jStart = 0;
        jEnd = m_nj - 1;
        kStart = 0;
        kEnd = m_nk - 1;
    }


    
void Grid_Zaran_3D::SetNi(int ni)
{
    m_ni = ni;
}

void Grid_Zaran_3D::SetNj(int nj)
{
    m_nj = nj;
}

void Grid_Zaran_3D::SetNk(int nk)
{
    m_nk = nk;
}

int Grid_Zaran_3D::GetNi()
{
    return m_ni;
}

int Grid_Zaran_3D::GetNj()
{
    return m_nj;
}

int Grid_Zaran_3D::GetNk()
{
    return m_nk;
}

void Grid_Zaran_3D::SetNodeNum(int ni, int nj, int nk)
{
    m_ni = ni;
    m_nj = nj;
    m_nk = nk;
}

void Grid_Zaran_3D::GetNodeNum(int& ni, int& nj, int& nk)
{
    ni = m_ni;
    nj = m_nj;
    nk = m_nk;
}

int Grid_Zaran_3D::GetNodeIndex(int i, int j, int k)
{
    return i + j * m_ni + k * m_ni * m_nj;
}

void Grid_Zaran_3D::GetNodeIndex(int index, int& i, int& j, int& k)
{
    i = index % m_ni;
    j = (index / m_ni) % m_nj;
    k = index / (m_ni * m_nj);
}

int Grid_Zaran_3D::GetCellIndex(int i, int j, int k)
{
     return i + j * (m_ni - 1) + k * (m_ni - 1) * (m_nj - 1);
}

void Grid_Zaran_3D::GetCellIndex(int index, int& i, int& j, int& k)
{
    i = index % (m_ni - 1);
    j = (index / (m_ni - 1)) % (m_nj - 1);
    k = index / ((m_ni - 1) * (m_nj - 1));
}

} // namespace zaran
