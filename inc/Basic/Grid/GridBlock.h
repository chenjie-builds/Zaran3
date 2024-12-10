#pragma once
#include "GridStruct.h"
#include "CommonPara.h"
namespace zaran
{
    enum class IBlank
    {
        Fluid = 0,
        Solid = 1,
        Trans = 2,
        Unset = -1,
    };
    class GridBlock : public GridStruct
    {
    public:
        GridBlock(const string &name, index_type index, index_type dim);
        ~GridBlock();
        void Allocate(index_type ni, index_type nj, index_type nk, index_type ghost_level) override;
        void SetDx(double dx) { m_dx = dx; }
        void SetDy(double dy) { m_dy = dy; }
        void SetDz(double dz) { m_dz = dz; }
        void SetBoundBox(const Box &box) { m_bound_box = box; }
        void SetIBlank(index_type i, index_type j, index_type k, IBlank iblank);
        const double &GetDx() const { return m_dx; }
        const double &GetDy() const { return m_dy; }
        const double &GetDz() const { return m_dz; }
        const Box &GetBoundBox() const { return m_bound_box; }
        const IBlank &GetIBlank(index_type i, index_type j, index_type k) const;

    public:
    private:
        double m_dx, m_dy, m_dz;
        Box m_bound_box;
		dynamic_array<IBlank> m_iblank;
    };
}