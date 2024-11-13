#include "Application.h"
#include "Log.h"
#include "BasicType.h"
#include <iostream>
int main()
{
	using namespace zaran;
	Logger::Start();
	Application* app = new Application();
	app->Run();
	delete app;
	Log::info("Simulation Finished! Program will exit in 10 seconds.");
	for (int i = 10; i > 0; i--)
	{
		spdlog::info("\r{} ", i);
		std::this_thread::sleep_for(std::chrono::seconds(1));
	}
	return 0;

}