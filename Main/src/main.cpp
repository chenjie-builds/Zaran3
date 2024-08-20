#include "Application.h"
#include "Log.h"
#include "BasicType.h"
#include <iostream>
int main()
{
	using namespace zaran;
	Logger::Start();
	Application *app = new Application();
	app->Run();
	delete app;
	Log::info("计算结束！按任意键退出……");
	std::cin.get();
}