#include"TimeCompute.h"
TimeCompute::TimeCompute()
{
	start_ = 0;
	end_ = 0;
}
void TimeCompute::Update(const double& delta_t)
{
	delta_ = delta_t;
	current_ += delta_;
}