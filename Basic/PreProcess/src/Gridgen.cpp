#include "Gridgen.h"
#include <fstream>
using namespace zaran;
Gridgen::Gridgen()
{
}
Gridgen::~Gridgen()
{
}
GridBase* zaran::Gridgen::CreateGrid()
{
	return nullptr;
}
void Gridgen::ReadGridgenCoord()
{

}

void Gridgen::SetGridFileName(std::string& fileName)
{
	gridFileName_ = fileName;
}


