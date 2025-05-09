#include "Application.h"
#include "Log.h"
#include "BasicType.h"
#include <iostream>
#include"File.h"
int main(int argc, char** argv)
{
	std::string work_dir;
	std::string control_file = "zaran.ini";
	if (argc == 1)
	{
		work_dir = zaran::GetCurPath();
		control_file = work_dir + "\\" + control_file;
		if (zaran::IsFileExist(control_file))
		{
			Log::info("Use default control file: {}", control_file);
		}
		else
		{
			Log::error("Default control file are not exist!");
			Log::warn("Use \"Zaran.exe path_to_work_dir\"");
			exit(0);
		}
	}
	else
	{
		work_dir = argv[1];
		control_file = work_dir + "\\" + control_file;
		Log::info("Control file:{}", control_file);
	}
	using namespace zaran;
	Logger::Start(work_dir);
	unique_ptr<Application> app = make_unique<Application>(work_dir);
	app->Run();
	for (int i = 10; i > 0; i--)
	{
		std::cout << "\r" << "Program will exit in " << i << " seconds.";
		std::this_thread::sleep_for(std::chrono::seconds(1));
	}
	return 0;

}