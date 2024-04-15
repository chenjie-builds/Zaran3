#include<cmath>
#include"Gas.h"
Gas::Gas(const double& Mw, const double& gamma, const Dimensionless& refValue)
{
	mw_ = Mw;
	Rm_ = GAS_CONSTANT / mw_ / refValue.GetRefRm();
}
