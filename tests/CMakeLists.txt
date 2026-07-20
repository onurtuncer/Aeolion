function(aeolion_add_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE aeolion)
    # Fixture files (e.g. the JSON handoff) are read from the source tree, so
    # tests do not depend on where the binary ends up.
    target_compile_definitions(${name} PRIVATE
        AEOLION_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/Data")
    add_test(NAME ${name} COMMAND ${name})
endfunction()

aeolion_add_test(TestSolverCore)
aeolion_add_test(TestBEMT)
aeolion_add_test(TestPropVane)
aeolion_add_test(TestDenseSolve)
aeolion_add_test(TestHandoffContract)

aeolion_add_test(TestSourcePanel)

aeolion_add_test(TestPanelBuilder)
# The only test that exercises compiled code rather than headers alone.
target_link_libraries(TestPanelBuilder PRIVATE aeolion_panelbuilder)
