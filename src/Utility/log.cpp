#include "log.h"
#include <cstdio>
#include <chrono>
Log& Log::Satrt()
{
	static Log log;
	return log;
}
Log::Log()
{
	auto console_sink = std::make_shared<ZaranLog::sinks::stdout_color_sink_mt>();
	console_sink->set_level(ZaranLog::level::debug);
	console_sink->set_pattern("[%D %H:%M:%S] [%^%l%$] %v");
	auto file_sink = std::make_shared<ZaranLog::sinks::basic_file_sink_mt>("log/log.txt", true);
	file_sink->set_level(ZaranLog::level::trace);
	file_sink->set_pattern("[%D %H:%M:%S] [%^%l%$] %v");
	ZaranLog::logger logger("multi_sink", { console_sink, file_sink });
	logger.set_level(ZaranLog::level::trace);
	logger.info("Start Logger, Log file path: ./{}", file_sink->filename());
	ZaranLog::set_default_logger(std::make_shared<ZaranLog::logger>(logger));
	ZaranLog::flush_every(std::chrono::seconds(1));//每隔1s刷新一次log文件
}

Log::~Log()
{
	ZaranLog::drop_all();
	ZaranLog::shutdown();
}