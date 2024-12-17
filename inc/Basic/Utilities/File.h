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
#include"BasicType.h"
namespace zaran
{
//判断文件存在
bool IsFileExist(const string& fileName);
//判断文件夹是否存在
bool IsFolderExist(const string& folderName);
//创建文件夹
bool CreateFolder(const string& folderName);
//复制文件
bool CopySingleFile(const string& src_file, const string& dest_file);
//删除文件
bool RemoveFile(const string& fileName);
string GetCurPath();

}