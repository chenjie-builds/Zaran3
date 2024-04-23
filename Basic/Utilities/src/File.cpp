#include"File.h"
#include"MathBasic.h"
#include<io.h>
#include <filesystem>
bool IsFileExist(const std::string& fileName)
{
	if (std::filesystem::exists(fileName))
		return true;
	return false;
}
bool IsFolderExist(const std::string& floderName)
{
	return IsFileExist(floderName);
}

bool CreateFolder(const std::string& folderName)
{
	if (IsFolderExist(folderName))
		return true;
	if (std::filesystem::create_directory(folderName))
		return true;
	return false;
}
