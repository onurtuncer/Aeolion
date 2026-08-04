// Logger/Log.h
//
// Thin wrapper around spdlog giving Aeolion two named loggers: a core
// logger for library-internal messages (solver, panelbuilder) and a
// client logger for application-side messages (viewer, app). Call
// Log::Init() once at program startup before using any logging macro;
// messages then go to a coloured console sink and to Aeolion.log.

#pragma once

#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>

#include <memory>

namespace Aeolion::Logger {

class Log
{
public:
    // Creates and registers both loggers. Safe to call more than once;
    // every call after the first is a no-op.
    static void Init();

    static std::shared_ptr<spdlog::logger>& GetCoreLogger() { return s_CoreLogger; }
    static std::shared_ptr<spdlog::logger>& GetClientLogger() { return s_ClientLogger; }

private:
    static std::shared_ptr<spdlog::logger> s_CoreLogger;
    static std::shared_ptr<spdlog::logger> s_ClientLogger;
};

} // namespace Aeolion::Logger

// Macros rather than inline functions so they stay callable unqualified
// from any namespace without a using-declaration, and so a future
// distribution build can compile them away entirely.

// Core logging macros (library-internal).
#define AE_CORE_TRACE(...)    ::Aeolion::Logger::Log::GetCoreLogger()->trace(__VA_ARGS__)
#define AE_CORE_INFO(...)     ::Aeolion::Logger::Log::GetCoreLogger()->info(__VA_ARGS__)
#define AE_CORE_WARN(...)     ::Aeolion::Logger::Log::GetCoreLogger()->warn(__VA_ARGS__)
#define AE_CORE_ERROR(...)    ::Aeolion::Logger::Log::GetCoreLogger()->error(__VA_ARGS__)
#define AE_CORE_CRITICAL(...) ::Aeolion::Logger::Log::GetCoreLogger()->critical(__VA_ARGS__)

// Client logging macros (application-side).
#define AE_TRACE(...)         ::Aeolion::Logger::Log::GetClientLogger()->trace(__VA_ARGS__)
#define AE_INFO(...)          ::Aeolion::Logger::Log::GetClientLogger()->info(__VA_ARGS__)
#define AE_WARN(...)          ::Aeolion::Logger::Log::GetClientLogger()->warn(__VA_ARGS__)
#define AE_ERROR(...)         ::Aeolion::Logger::Log::GetClientLogger()->error(__VA_ARGS__)
#define AE_CRITICAL(...)      ::Aeolion::Logger::Log::GetClientLogger()->critical(__VA_ARGS__)
