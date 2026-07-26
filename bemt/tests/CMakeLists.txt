aeolion_add_test(TestBEMT)
aeolion_add_test(TestPropContract)
aeolion_add_test(TestPropVane)
# These exercise compiled code rather than headers alone.
target_link_libraries(TestBEMT PRIVATE aeolion_bemt)
target_link_libraries(TestPropContract PRIVATE aeolion_bemt)
target_link_libraries(TestPropVane PRIVATE aeolion_bemt)
