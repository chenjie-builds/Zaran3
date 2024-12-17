#include"File.h"
#include"MathBasic.h"
#include <filesystem>
namespace zaran
{

	bool IsFileExist(const string& fileName)
	{
		if (std::filesystem::exists(fileName))
			return true;
		return false;
	}
	bool IsFolderExist(const string& floderName)
	{
		return IsFileExist(floderName);
	}

	bool CreateFolder(const string& folderName)
	{
		if (IsFolderExist(folderName))
			return true;
		if (std::filesystem::create_directory(folderName))
			return true;
		return false;
	}

	bool CopySingleFile(const string& src_file, const string& dest_file)
	{
		std::filesystem::copy_file(src_file, dest_file);
		return true;
	}
	bool RemoveFile(const string& fileName)
	{
		if (std::filesystem::remove(fileName))
			return true;
		return false;
	}
	string GetCurPath()
	{
		return std::filesystem::current_path().string();
	}
}