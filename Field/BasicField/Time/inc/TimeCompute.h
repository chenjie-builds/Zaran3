//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	TimeCompute.h														||
//*	@brief	计算时间类, 用于统计计算物理时间(后面可能需要修改)						||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
class TimeCompute
{
public:
	TimeCompute();
	TimeCompute(const double& start, const double& end) :start_(start), end_(end),delta_(0),current_(0) {};
	void Update(const double& delta_t);
	double GetCurrentPhysicsTime()const { return current_; }
	const double& GetDeltaT()const { return delta_; }
	double GetEndTime()const { return end_; }
private:
	double start_;
	double end_;
	double delta_;
	double current_;
};