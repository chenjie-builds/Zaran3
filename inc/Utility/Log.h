//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Log.h																||
//*	@brief	使用spdlog进行日志记录												||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
//对spdlog进行初始化
//TODO 对spdlog库更熟悉后： 设置输出位置，log等级自定义
namespace ZaranLog = spdlog;
class Log
{
public:
	static Log& Satrt();
private:
	Log();
	~Log();
};




