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
	Log::info("Simulation Finished!");
	std::cin.get();
}