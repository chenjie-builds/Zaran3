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
        GridBlock(const string &name, Id index, Id dim);
        ~GridBlock();
        void Allocate(Id ni, Id nj, Id nk, Id ghost_level) override;
        void SetDx(double dx) { m_dx = dx; }
        void SetDy(double dy) { m_dy = dy; }
        void SetDz(double dz) { m_dz = dz; }
        void SetBoundBox(const Box &box) { m_bound_box = box; }
        void SetIBlank(Id i, Id j, Id k, IBlank iblank);
        const double &GetDx() const { return m_dx; }
        const double &GetDy() const { return m_dy; }
        const double &GetDz() const { return m_dz; }
        const Box &GetBoundBox() const { return m_bound_box; }
        const IBlank &GetIBlank(Id i, Id j, Id k) const;

    public:
    private:
        double m_dx, m_dy, m_dz;
        Box m_bound_box;
		Array<IBlank> m_iblank;
    };
}