//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Solver_NS_3D_Zaran.h															||
//*	@brief	三维NS 方程求解器														||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "NSSolver.h"
#include "Grid_Zaran_3D.h"
namespace zaran
{
    class Solver_NS_3D_Zaran :public NSSolver
    {
    public:
        Solver_NS_3D_Zaran() {}
        ~Solver_NS_3D_Zaran() {}
        void InitField()override;
        void ComputeCoordTrans()override;
        void CreateFieldData()override;
        void RegisterFieldData()override;
    protected:
        void ComputeGradientWLS()override;
        void ComputeTimeStepLocal() override;
        // 计算流动通量
        void InviscidFlux()override;
        //计算粘性通量
        void ViscousFlux() override;
        //计算源项
        void SourceFlux() override;
        void ComputeLimiterCoefVK()override;
        void ComputeLimiterCoefBJ()override;
        void ComputeLimiterCoefNoLimiter()override;
        void ComputeLimiterCoefOneOrder()override;
        void ComputeBoundaryLimiterCoef()override;

        void RungeKutta()override;
    protected:
        void ComputeCoordTransStruct();//计算结构网格部分的坐标变换系数
        void ComputeCoordTransMid();//计算过渡网格部分的坐标变换系数
        void ComputeGradientWLSStruct();//计算结构网格部分的梯度
        void ComputeGradientWLSMid();//计算过渡网格部分的梯度
        void ComputeTimeStepLocalStruct();//计算结构网格部分的局部时间步长
        void ComputeTimeStepLocalMid();//计算过渡网格部分的局部时间步长
        void InviscidFluxStruct();//计算结构网格部分的流动通量
        void InviscidFluxMid();//计算过渡网格部分的流动通量
        void ViscousFluxStruct();//计算结构网格部分的粘性通量
        void ViscousFluxMid();//计算过渡网格部分的粘性通量
        void SourceFluxStruct();//计算结构网格部分的源项
        void SourceFluxMid();//计算过渡网格部分的源项

    protected:
        void BoundaryCondition()override;

    private:
        Ptr<Grid_Zaran_3D> GetGrid() { return std::static_pointer_cast<Grid_Zaran_3D>(gridList_->GetGrid(gridIndex_)); }
    private:
        Array<DArray*> m_prim_bound;
        Array<DArray*> m_prim_bound_gradX;
        Array<DArray*> m_prim_bound_gradY;
        Array<DArray*> m_prim_bound_gradZ;
        Array<DArray*> m_limiter_bound;//边界节点限制器系数

    };
}