#include "FaceTopoInfo.h"
namespace zaran
{


	FaceTopoInfo::FaceTopoInfo()
	{

	}

	FaceTopoInfo::~FaceTopoInfo()
	{

	}

	void FaceTopoInfo::SetFace2Node(std::vector<int>& face2node)
	{
		face2node_ = face2node;
	}

	void FaceTopoInfo::SetLeftCell(int& leftCell)
	{
		leftCell_ = leftCell;
	}

	void FaceTopoInfo::SetRightCell(int& rightCell)
	{
		rightCell_ = rightCell;
	}

	std::vector<int>& FaceTopoInfo::GetFace2Node()
	{
		return face2node_;
	}

	int& FaceTopoInfo::GetLeftCell()
	{
		return leftCell_;
	}

	int& FaceTopoInfo::GetRightCell()
	{
		return rightCell_;
	}

}