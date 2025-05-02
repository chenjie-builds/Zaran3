/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file Log.h
 * \brief Using spdlog for logging.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */

#pragma once
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
//对spdlog进行初始化
//TODO 对spdlog库更熟悉后： 设置输出位置，log等级自定义
namespace Log = spdlog;
class Logger
{
public:
	static Logger& Start(const std::string& work_dir);
private:
	Logger(const std::string& work_dir);
	~Logger();
};




