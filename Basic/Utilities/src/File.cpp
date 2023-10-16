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
bool IsFloderExist(const std::string& floderName)
{
	return IsFileExist(floderName);
}
