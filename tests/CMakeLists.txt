function(aeolion_add_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE aeolion)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

aeolion_add_test(TestVlmCore)
aeolion_add_test(TestPanelMethod)
aeolion_add_test(TestBoundaryLayer)
aeolion_add_test(TestBemt)
aeolion_add_test(TestPropVane)
aeolion_add_test(TestSolverBackend)

# TestMeshSlice pulls in the shared fixtures/SyntheticWing.h
add_executable(TestMeshSlice TestMeshSlice.cpp)
target_link_libraries(TestMeshSlice PRIVATE aeolion)
target_include_directories(TestMeshSlice PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME TestMeshSlice COMMAND TestMeshSlice)
