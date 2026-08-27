#include "Application.h"
#include "Log.h"
#include "BasicType.h"
#include "ZaranError.h"
#include <iostream>
#include <thread>
#include <chrono>
#include <exception>
#include "File.h"

int main(int argc, char** argv)
{
    try
    {
        std::string work_dir;
        std::string control_file = "zaran.toml";
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
                throw zaran::ZaranError("Default control file zaran.toml not found");
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
    catch (const zaran::ZaranError& e)
    {
        Log::error("Fatal error: {}", e.what());
        std::cerr << "Fatal error: " << e.what() << std::endl;
        return 1;
    }
    catch (const std::exception& e)
    {
        Log::error("Unexpected exception: {}", e.what());
        std::cerr << "Unexpected exception: " << e.what() << std::endl;
        return 1;
    }
}
