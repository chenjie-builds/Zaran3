//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Timer.h																||
//*	@brief	计时器类, 用于计算程序计算时间											||
//*	@author	Chen Jie.															||
//==============================================================================||
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