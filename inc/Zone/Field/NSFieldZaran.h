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
		void SetModelManager(std::shared_ptr <ModelManager> model_manager);
		std::shared_ptr<GridBlock>GetGrid();
		std::shared_ptr<FlowSolverParamStruct>GetSolverPara();
		std::shared_ptr<DataManagerNSStruct>GetDataManager();
		std::shared_ptr <ModelManager> GetModelManager() { return m_model_manager; }
		std::shared_ptr <NSFieldFNFDM> GetSlaveField() { return m_slave_field; }
		void CalcResidual() override;
		//-----------------------------Master--------------------------------------
	public:
		void DeleteSlaveField(std::shared_ptr<FieldManager> field_manager);
		void CreateSlaveField(std::shared_ptr<FieldManager> field_manager);

	protected:
		/// @brief slave field is managed by master field
		/// @details using voxel information of master field to generate slave field
		/// @details also can be deleted by master struct field
		std::shared_ptr <NSFieldFNFDM> m_slave_field;
		/// @brief
		std::shared_ptr <ModelManager> m_model_manager;
		std::shared_ptr <IdProxyStruct> m_idx_proxy;
	};
}