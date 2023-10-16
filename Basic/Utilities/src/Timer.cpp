#include"Timer.h"
void Timer::StartRecord()
{
	start = std::chrono::high_resolution_clock::now();
	record = end = start;
}
float Timer::DurationFromLastRecord()
{
	std::chrono::duration<float>duration;
	duration = end - record;
	float s = duration.count();
	return s;
}
float Timer::DurationFromStart()
{
	std::chrono::duration<float>duration;
	duration = end - start;
	float s = duration.count();
	return s;
}
void Timer::Record()
{
	record = end;
	end = std::chrono::high_resolution_clock::now();
}