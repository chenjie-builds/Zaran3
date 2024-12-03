#pragma once
#include "NSField.h"
#include "GridBlock.h"
#include "GridFNFDM.h"
#include "NSFieldFN.h"
#include "FlowSolverStructPara.h"
#include "DataManagerNSStruct.h"
#include "StructIdxProxy.h"
#include "ModelManager.h"
#include "GlobalField.h"
#include <set>
namespace zaran
{
    class NSFieldZaran : public FieldNS
    {
    public:
        NSFieldZaran();
        ~NSFieldZaran();
        void Allocate() override;

    protected:
        void AllocateSolver() override;
        void AllocateDataManager() override;
        void AllocateIdxProxy();
        void AllocateSolverPara() override;

    public:
        void SetModelManager(ModelManager *model_manager);
        GridBlock *GetGrid() override;
        FlowSolverParamStruct *GetSolverPara() override;
        DataManagerNSStruct *GetDataManager() override;
        ModelManager *GetModelManager() { return m_model_manager; }
        NSFieldFNFDM *GetSlaveField() { return m_slave_field; }
        void CalcResidual() override;
        //-----------------------------Master--------------------------------------
    public:
        void DeleteSlaveField(FieldManager *field_manager);
        void CreateSlaveField(FieldManager *field_manager);

    protected:
        /// @brief slave field is managed by master field
        /// @details using voxel information of master field to generate slave field
        /// @details also can be deleted by master struct field
        NSFieldFNFDM *m_slave_field;
        /// @brief
        ModelManager *m_model_manager;
        IdxStruct *m_idx_proxy;
    };
}