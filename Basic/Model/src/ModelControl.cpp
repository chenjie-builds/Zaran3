#include "ModelControl.h"
#include "Log.h"
#include<fstream>
#include<iostream>
using namespace zaran;
ModelControl::ModelControl(string filename)
{
	SetFileName(filename);
}


void ModelControl::SetFileName(string filename)
{
	m_filename = filename;
}

void ModelControl::Input()
{
	double x, y;
	vector<DVector3D>vertex;

	std::ifstream fin(m_filename);
	//用于扎染演示
	//string tmpstr;
	//while (!fin.eof())
	//{
	//	fin >> x;
	//	fin >> y;
	//	vertex.push_back(Point(x, y));
	//}
	//if (vertex.size() != 0)
	//{
	//	modelNum_ = 1;
	//	modelVec_.resize(modelNum_);
	//	modelVec_[0] = new ModPointCloud(vertex);
	//}
	//else
	//{
	//	modelNum_ = 0;
	//	modelVec_.resize(modelNum_);
	//}
	//delta = 0.001;
	//return;
	//原始设定的不同外形读取
	fin >> modelNum_;
	modelVec_.resize(modelNum_);
	ZaranLog::info("Model Num= {}", modelNum_);
	string modeType;
	int pointNum;
	for (int iModel = 0; iModel < modelNum_; ++iModel)
	{


		fin >> modeType;
		if (modeType == "Point")
		{
			fin >> pointNum;
			vertex.resize(pointNum);
			for (int iPoint = 0; iPoint < pointNum; ++iPoint)
			{
				fin >> x;
				fin >> y;
				vertex[iPoint] = { x,y,0 };
			}
			modelVec_[iModel] = new PointCloudModel(vertex);
		}
	}
	delta = 0.01;
}


double ModelControl::Delta()
{
	return delta;
}
