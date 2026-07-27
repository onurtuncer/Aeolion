// TestLog.cpp -- standing regression check on the logger's contract:
// Init() creates and registers both named loggers, calling it again is a
// harmless no-op (spdlog would throw on a duplicate registration), the
// AE_* / AE_CORE_* macros compile and run from client code, and messages
// actually reach the Aeolion.log file sink.
#include "Aeolion/Logger/Log.h"

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

using Aeolion::Logger::Log;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

int main() {
    Log::Init();
    CHECK(Log::GetCoreLogger() != nullptr, "core logger created");
    CHECK(Log::GetClientLogger() != nullptr, "client logger created");
    CHECK(spdlog::get("AEOLION") == Log::GetCoreLogger(), "core logger registered as AEOLION");
    CHECK(spdlog::get("APP") == Log::GetClientLogger(), "client logger registered as APP");

    auto coreBefore = Log::GetCoreLogger();
    Log::Init();
    CHECK(Log::GetCoreLogger() == coreBefore, "second Init() is a no-op");

    AE_CORE_INFO("core message {}", 42);
    AE_INFO("client message {}", 3.5);
    Log::GetCoreLogger()->flush();
    Log::GetClientLogger()->flush();

    std::ifstream logFile("Aeolion.log");
    CHECK(logFile.good(), "Aeolion.log created in the working directory");
    std::string contents((std::istreambuf_iterator<char>(logFile)), std::istreambuf_iterator<char>());
    CHECK(contents.find("core message 42") != std::string::npos, "core message reached the file sink");
    CHECK(contents.find("client message 3.5") != std::string::npos, "client message reached the file sink");

    if (failures == 0)
        std::cout << "TestLog: all checks passed\n";
    return failures;
}
