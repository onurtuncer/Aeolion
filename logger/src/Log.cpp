// Logger/Log.cpp
//
// Builds the two named loggers declared in Log.h on a shared pair of
// sinks: a coloured console sink and a truncating file sink (Aeolion.log
// in the working directory). Debug builds log everything from trace up;
// release builds start at info.

#include "Aeolion/Logger/Log.h"

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>

#include <vector>

namespace Aeolion::Logger {

std::shared_ptr<spdlog::logger> Log::s_CoreLogger;
std::shared_ptr<spdlog::logger> Log::s_ClientLogger;

void Log::Init()
{
    if (s_CoreLogger)
        return;

    std::vector<spdlog::sink_ptr> logSinks;
    logSinks.emplace_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
    logSinks.emplace_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>("Aeolion.log", true));

    logSinks[0]->set_pattern("%^[%T] %n: %v%$");
    logSinks[1]->set_pattern("[%T] [%l] %n: %v");

#ifdef NDEBUG
    constexpr auto runtimeLevel = spdlog::level::info;
#else
    constexpr auto runtimeLevel = spdlog::level::trace;
#endif

    s_CoreLogger = std::make_shared<spdlog::logger>("AEOLION", begin(logSinks), end(logSinks));
    spdlog::register_logger(s_CoreLogger);
    s_CoreLogger->set_level(runtimeLevel);
    s_CoreLogger->flush_on(runtimeLevel);

    s_ClientLogger = std::make_shared<spdlog::logger>("APP", begin(logSinks), end(logSinks));
    spdlog::register_logger(s_ClientLogger);
    s_ClientLogger->set_level(runtimeLevel);
    s_ClientLogger->flush_on(runtimeLevel);
}

} // namespace Aeolion::Logger
