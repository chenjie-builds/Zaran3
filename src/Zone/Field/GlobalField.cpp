#include "GlobalField.h"
namespace zaran
{
	void FieldManager::RemoveField(std::shared_ptr<Field> field)
	{
		m_field.erase(m_field.begin() + field->GetIdx());
		m_field_data_comm_info.erase(m_field_data_comm_info.begin() + field->GetIdx());
	}
	FieldManager::FieldManager()
	{

	}
	FieldManager::~FieldManager()
	{

	}
	std::shared_ptr<Field> FieldManager::GetField(Id idx_field)
	{
		return m_field[idx_field];
	}
	std::shared_ptr<FieldDataCommInfo> FieldManager::GetFieldDataCommInfo(Id idx_field)
	{
		return m_field_data_comm_info[idx_field];
	}
	void FieldManager::SetFieldDataCommInfo(Id idx_field, std::shared_ptr<FieldDataCommInfo>field_data_comm_info)
	{
		m_field_data_comm_info[idx_field] = field_data_comm_info;
	}
	Id FieldManager::GetFieldNum()
	{
		return m_field.size();
	}
	void FieldManager::AddField(std::shared_ptr<Field> field, std::shared_ptr<FieldDataCommInfo> field_data_comm_info)
	{
		field->SetIdx(m_field.size());
		m_field.push_back(field);
		m_field_data_comm_info.push_back(field_data_comm_info);
	}
} // namespace zaran