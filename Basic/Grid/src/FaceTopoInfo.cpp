#include "FaceTopoInfo.h"
namespace zaran
{
	FaceTopoInfo::FaceTopoInfo()
	{

	}

	FaceTopoInfo::~FaceTopoInfo()
	{

	}
	void FaceTopoInfo::SetNodeIndex(Array<IArray>& node_index)
	{
		m_node_index = node_index;
	}
	void FaceTopoInfo::SetLeftCell(IArray& left_cell_index)
	{
		m_left_cell_index = left_cell_index;
	}
	void FaceTopoInfo::SetRightCell(IArray& right_cell_index)
	{
		m_right_cell_index = right_cell_index;
	}
}