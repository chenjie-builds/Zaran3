//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	ModelContol.h														||
//*	@brief	模型控制类： 输入控制模型的参数，并对模型进行修改						||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include<vector>
#include<string>
#include"Model.h"
#include"segment.h"
#include"PolygonModel.h"
#include"PointCloudModel.h"
using std::vector;
using std::string;
class ModelControl
{
private:
	int modelNum_;//模型个数
	vector<Model*>modelVec_;//模型数组
	double delta;//模型生成散点的间距
	string m_filename;//模型控制文件的文件名
public:
	ModelControl() {};
	ModelControl(string filename);
	void SetFileName(string filename);
	void Input();
	vector<Model*> Model() { return modelVec_; };
	double Delta();
};