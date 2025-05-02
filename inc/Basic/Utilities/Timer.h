/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file Timer.h
 * \brief Timer class, used to calculate the time of program.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */

#pragma once
#include<chrono>
class Timer
{
private:
	std::chrono::high_resolution_clock::time_point start;
	std::chrono::high_resolution_clock::time_point record;
	std::chrono::high_resolution_clock::time_point end;
public:
	Timer() {}
	void StartRecord();
	float DurationFromLastRecord();
	float DurationFromStart();
	void Record();
};