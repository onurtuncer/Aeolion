aeolion_add_test(TestLog)
# Exercises compiled code rather than headers alone.
target_link_libraries(TestLog PRIVATE aeolion_logger)
