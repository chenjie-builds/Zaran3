#include "RefValue.h"
#include"CommonPara.h"
#include"MathBasic.h"
#include "Gas.h"
#include<cmath>
RefValue::RefValue()
{
	gama_ = 1.4;
	mw_ = 0.02896;
	length_ = 1.0;
	density_ = 1.0;
	temp_ = 2.488917144919069e+003;
	Rm_ = gama_ * GAS_CONSTANT / mw_;
	velocity_ = sqrt(gama_ * GAS_CONSTANT / mw_ * temp_);//当地声速
	pressure_ = density_ * velocity_ * velocity_;
	time_ = length_ / velocity_;
	mass_ = density_ * length_ * length_ * length_;
	Re_ = density_ * velocity_ * length_ / Southerland(temp_);
	double mu0 = 1.711e-5 / density_ / velocity_ / length_;
	double T0 = 273.0 / temp_;
	double Ts = 110.4 / temp_;
	Re_ = 1.0 / Southerland(1.0, mu0, T0, Ts) * length_;

}
RefValue::RefValue(const double& refDensity, const double& refGamma, const double& refMw, const double& refLength, const double& refTemp)
{
	gama_ = refGamma;
	mw_ = refMw;
	length_ = refLength;
	density_ = refDensity;
	temp_ = refTemp;
	velocity_ = sqrt(gama_ * GAS_CONSTANT / mw_ * temp_);//当地声速
	Rm_ = velocity_ * velocity_ / temp_;
	pressure_ = density_ * velocity_ * velocity_;
	time_ = length_ / velocity_;
	mass_ = density_ * length_ * length_ * length_;
	Re_ = density_ * velocity_ * length_ / Southerland(temp_);
	double mu0 = 1.711e-5 / density_ / velocity_ / length_;
	double T0 = 273.0 / temp_;
	double Ts = 110.4 / temp_;
	Re_ = 1.0 / Southerland(1.0, mu0, T0, Ts) * length_;
}
