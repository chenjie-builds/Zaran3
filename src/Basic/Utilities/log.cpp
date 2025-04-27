#include "Log.h"
#include <cstdio>
#include <chrono>
Logger& Logger::Start(const std::string& work_dir)
{
	static Logger log(work_dir);
	return log;
}
Logger::Logger(const std::string& work_dir)
{
	auto console_sink = std::make_shared<Log::sinks::stdout_color_sink_mt>();
	console_sink->set_level(Log::level::debug);
	console_sink->set_pattern("[%D %H:%M:%S] [%^%l%$] %v");
	auto file_sink = std::make_shared<Log::sinks::basic_file_sink_mt>(work_dir+"\\log.txt", false);
	file_sink->set_level(Log::level::trace);
	file_sink->set_pattern("[%D %H:%M:%S] [%^%l%$] %v");
	Log::logger logger("multi_sink", { console_sink, file_sink });
	logger.set_level(Log::level::trace);
	logger.info("Start Logger");
	logger.info("Log file path: {}", file_sink->filename());
	Log::set_default_logger(std::make_shared<Log::logger>(logger));
	Log::flush_every(std::chrono::milliseconds(100));//每隔0.1s刷新一次log文件
}

Logger::~Logger()
{
	Log::drop_all();
	Log::shutdown();
}