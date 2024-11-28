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
	for (int i = 10; i > 0; i--)
	{
		std::cout << "\r" << "Program will exit in " << i << " seconds.";
		std::this_thread::sleep_for(std::chrono::seconds(1));
	}
	return 0;

}