# Aeolion logger: two named spdlog loggers (core = library-internal,
# client = application) behind the AE_CORE_* / AE_* macros in Logger/Log.h.
#
# STATIC, same rationale as panelbuilder/CMakeLists.txt: no ABI-stability
# or plugin boundary to justify a DLL's __declspec plumbing.

find_package(spdlog CONFIG REQUIRED)

add_library(aeolion_logger STATIC
    src/Log.cpp
)

target_include_directories(aeolion_logger PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# spdlog is PUBLIC because Log.h exposes spdlog types and the macros expand
# to spdlog calls inside the consumer's own translation unit. aeolion is
# PUBLIC for the same reason as in panelbuilder/: it carries the warnings/standard
# flags on its INTERFACE.
target_link_libraries(aeolion_logger PUBLIC aeolion spdlog::spdlog)

if(AEOLION_BUILD_TESTS)
    add_subdirectory(tests)
endif()
