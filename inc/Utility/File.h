//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	File.h																||
//*	@brief	文件相关操作															||
//*	@author	Chen Jie.															||
//==============================================================================||
#include<string>
#include<Eigen\Dense>
//判断文件存在
bool IsFileExist(const std::string& fileName);
//判断文件夹是否存在
bool IsFloderExist(const std::string& floderName);
