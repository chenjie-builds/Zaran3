/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file File.h
 * \brief File related operations.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */

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