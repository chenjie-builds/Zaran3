#include"Simulation.h"
#include "Log.h"
#include"BasicType.h"
#include"Test.h"
#include<iostream>
int main()
{
	using namespace zaran;
	Logger::Start();
	Simulation* simulation = new Simulation();
	simulation->Run();
	delete[] simulation;
	Log::info("计算结束！按任意键退出……");
	std::cin.get();
}