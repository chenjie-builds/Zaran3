#include "FaceTopoInfo.h"
using namespace zaran;
FaceTopoInfo::FaceTopoInfo()
{

}

FaceTopoInfo::~FaceTopoInfo()
{

}

void FaceTopoInfo::SetFace2Node(IArray& face2node)
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

IArray& FaceTopoInfo::GetFace2Node()
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

