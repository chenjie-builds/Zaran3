#pragma once
#include "BasicType.h"
namespace zaran
{
    /// @brief bounary information for structred grid
    class BoundStruct
    {
    public:
        /// @brief constructor
        /// @param i_bnd the index of the boundary in i direction
        /// @param j_bnd the index of the boundary in j direction
        /// @param k_bnd the index of the boundary in k direction
        /// @param direction a 3-element array, each element is -1, 0 or 1
        /// @param norm a 3-element array, representing the norm vector of the boundary in 3D space
        /// @param i_bnd_tgt the index of the boundary in i direction of the target block
        /// @param j_bnd_tgt the index of the boundary in j direction of the target block
        /// @param k_bnd_tgt the index of the boundary in k direction of the target block
        /// @param direction_tgt a 3-element array, each element is -1, 0 or 1
        BoundStruct(int i_bnd, int j_bnd, int k_bnd, const int *direction, const double *norm, int i_bnd_tgt = -1, int j_bnd_tgt = -1, int k_bnd_tgt = -1, const int *direction_tgt = nullptr, int target_block = -1);
        BoundStruct();
        ~BoundStruct();
        /// @brief set the index of the boundary
        /// @param i_bnd the index of the boundary in i direction
        /// @param j_bnd the index of the boundary in j direction
        /// @param k_bnd the index of the boundary in k direction
        void SetIdx(const int &i_bnd, const int &j_bnd, const int &k_bnd);
        /// @brief set direction vector of the boundary
        /// @param direction a 3-element array, each element is -1, 0 or 1
        /// @details -1: left boundary, 0: not boundary, 1: right boundary
        void SetDirection(const int *direction);
        /// @brief set norm vector of the boundary
        /// @param norm a 3-element array, representing the norm vector of the boundary in 3D space
        void SetNorm(const double *norm);
        /// @brief get the index of the boundary
        /// @param i_bnd the index of the boundary in i direction
        /// @param j_bnd the index of the boundary in j direction
        /// @param k_bnd the index of the boundary in k direction
        void GetIdx(int &i_bnd, int &j_bnd, int &k_bnd);
        /// @brief  get direction vector of the boundary
        /// @return a 3-element array, each element is -1, 0 or 1
        const int *GetDirection();
        /// @brief get norm vector of the boundary
        /// @return a 3-element array, representing the norm vector of the boundary in 3D space
        const double *GetNorm();
        /// @brief compare two boundaries, for adding boundary to the boundary manager
        bool operator==(const BoundStruct &bound);
        /// @brief set the index of the boundary in the target block
        /// @param i_bnd_tgt the index of the boundary in i direction of the target block
        /// @param j_bnd_tgt the index of the boundary in j direction of the target block
        /// @param k_bnd_tgt the index of the boundary in k direction of the target block
        void SetIdxTgt(const int &i_bnd_tgt, const int &j_bnd_tgt, const int &k_bnd_tgt);
        /// @brief set direction vector of the boundary in the target block
        /// @param direction_tgt a 3-element array, each element is -1, 0 or 1
        /// @details -1: left boundary, 0: not boundary, 1: right boundary
        void SetDirectionTgt(const int *direction_tgt);
        /// @brief get the index of the boundary in the target block
        /// @param i_bnd_tgt the index of the boundary in i direction of the target block
        /// @param j_bnd_tgt the index of the boundary in j direction of the target block
        /// @param k_bnd_tgt the index of the boundary in k direction of the target block
        void GetIdxTgt(int &i_bnd_tgt, int &j_bnd_tgt, int &k_bnd_tgt);
        /// @brief  get direction vector of the boundary in the target block
        /// @return a 3-element array, each element is -1, 0 or 1
        const int *GetDirectionTgt();
        void SetTargetBlock(const int &target_block);
        int GetTargetBlock();

    private:
        // 边界的索引，用于标记边界的位置, src(source)表示当前块，tgt(target)表示相邻块
        int m_i_bnd_src, m_j_bnd_src, m_k_bnd_src;
        // norm vector of the boundary
        double m_norm[3];
        // 边界的方向，用于标记边界的方向，只有一个方向为+1或-1，其余为0
        //  -1：左侧边界
        //  +1：右侧边界
        // example: {1,0,0}表示i方向上的右侧边界
        // src(source)表示当前块，tgt(target)表示相邻块
        int m_dir_src[3];

        // target block index
        int m_target_block;
        // i,j,k index of the boundary in the target block
        int m_i_bnd_tgt, m_j_bnd_tgt, m_k_bnd_tgt;
        // 边界的方向，用于标记边界的方向，只有一个方向为+1或-1，其余为0
        //  -1：左侧边界
        //  +1：右侧边界
        // example: {1,0,0}表示i方向上的右侧边界
        // src(source)表示当前块，tgt(target)表示相邻块
        int m_dir_tgt[3];
    };
}