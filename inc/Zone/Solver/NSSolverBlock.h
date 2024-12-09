#pragma once
#include "NSSolverStructDEER.h"
#include "GridBlock.h"
namespace zaran
{
    class NSSolverBlock : public NSSolverStructDEER
    {
    public:
        NSSolverBlock(Id index, string name, std::shared_ptr<FlowSolverParam> para, std::shared_ptr < GridBlock> grid, std::shared_ptr < DataManagerNSStruct>data_manager);
        ~NSSolverBlock();

    protected:
        std::shared_ptr<GridBlock>GetGrid();

        void FluxDifference2nd() override;
        /// @brief mid node value interpolate
        /// @param idx_left left node index
        /// @param idx_right right node index
        /// @param lef_coord coordinate of left point
        /// @param mid_coord coordinate of mid point
        /// @param right_coord coordinate of right point
        /// @param value_left interpolated value of left side at mid point
        /// @param value_right interpolated value of right side at mid point
        virtual void InterMidNodePrim_Grad(int idx_left, int idx_right, const double *lef_coord, const double *mid_coord, const double *right_coord, double *value_left, double *value_right)override;
    };
} // namespace zaran