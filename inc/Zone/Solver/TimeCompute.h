/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file TimeCompute.h
 * \brief TimeCompute class, used to compute time.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
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