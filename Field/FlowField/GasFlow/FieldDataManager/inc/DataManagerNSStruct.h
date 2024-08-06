#pragma once 
#include"DataManagerNS.h"

namespace zaran
{
    class DataManagerNSStruct :public DataManagerNS
    {
    public:
        DataManagerNSStruct(FieldData* fieldData, int ni, int nj, int nk);
        ~DataManagerNSStruct();
    public:
        void CreateData()override;
        void RegisterData()override;
    public:
        //! @brief 设置半点左侧基本变量
        /// @param iEqu 变量索引
        /// @param iDim 方向索引
        /// @param iNode 节点索引
        /// @param value 输入的值
        void SetMidNodePrimLeft(int iEqu,int iDim, int iNode,  double value);
        //! @brief 设置半点右侧基本变量
        /// @param iEqu 变量索引
        /// @param iDim 方向索引
        /// @param iNode 节点索引
        /// @param value 输入的值
        void SetMidNodePrimRight(int iEqu, int iDim,int iNode,  double value);
        void SetMidNodePrim(int iEqu, int iDim, int iNode, double value_left, double value_right);
        double GetMidNodePrimLeft(int iEqu, int iDim, int iNode);
        double GetMidNodePrimRight(int iEqu, int iDim, int iNode);
    private:
    /// @brief 半点左侧的基本变量，1D变量索引，2D方向索引(1:i+1,2:j+1,3:k+1)，3D节点索引
    double***m_midnode_prim_left;
    /// @brief 半点右侧的基本变量，1D变量索引，2D方向索引(1:i+1,2:j+1,3:k+1)，3D节点索引
    double***m_midnode_prim_right;
    };
}