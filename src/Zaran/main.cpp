#include"Simulation.h"
#include<memory>
#include "log.h"
#include"ReadSTL.h"
int main()
{
	using namespace zaran;
	Log::Satrt();
	STLReader reader;
	reader.ReadSTLFile("D:\\hellfire.stl");
	std::shared_ptr< Simulation > simulation = std::make_shared<Simulation>();
	simulation->Start();
	ZaranLog::info("计算结束！按任意键退出……");
	system("pause");
}