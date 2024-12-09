#include "NSField.h"
namespace zaran
{
	FieldNS::FieldNS(std::shared_ptr < GridBase> grid, FieldType fieldType) :Field(grid, fieldType)
	{
		m_res_info = nullptr;
	}
	FieldNS::~FieldNS()
	{
	}
	void FieldNS::AllocateSolverPara()
	{
		m_solver_para = std::make_shared<FlowSolverParam>();
		GetSolverPara()->Init();
	}
	void FieldNS::Allocate()
	{
		Field::Allocate();
		AllocateResInfo();
	}

	std::shared_ptr<zaran::FlowSolverParam> FieldNS::GetSolverPara()
	{
		return std::static_pointer_cast<FlowSolverParam>(Field::GetSolverPara());
	}

	std::shared_ptr<zaran::DataManagerNS> FieldNS::GetDataManager()
	{
		return std::static_pointer_cast<DataManagerNS>(Field::GetDataManager());
	}

	void FieldNS::AllocateResInfo()
	{
		if (m_res_info != nullptr)
		{
			delete m_res_info;
		}
		m_res_info = new ResInfo(GetSolverPara()->GetEqNum());
	}

}