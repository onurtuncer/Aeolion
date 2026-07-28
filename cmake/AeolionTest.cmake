# AeolionTest.cmake
#
# aeolion_add_test(<name>) -- registers <name>.cpp (in the CALLING
# directory) as a ctest executable linked against the header-only `aeolion`
# target. Shared by every module's tests/CMakeLists.txt (solver/,
# panelbuilder/, and the top-level tests/ for the core header-only
# library) so the boilerplate -- and the fixture-directory path -- is
# defined once rather than once per module.
#
# Fixture JSON lives in ONE place (tests/Data/, at the project root) no
# matter which module's test consumes it, so AEOLION_TEST_DATA_DIR is
# anchored at CMAKE_SOURCE_DIR (the top project source dir) rather than
# CMAKE_CURRENT_SOURCE_DIR (which would resolve to a different, wrong
# directory depending on which module's tests/CMakeLists.txt called this).
function(aeolion_add_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE aeolion)
    target_compile_definitions(${name} PRIVATE
        AEOLION_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/Data")
    add_test(NAME ${name} COMMAND ${name})
endfunction()
