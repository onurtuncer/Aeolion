# The DAVE-ML flight model's verifier.
#
# The .dml is a GENERATED artifact, and a generated file that nothing
# re-reads is the one thing in this repo that can rot silently. So the
# verifier runs as an ordinary ctest: it parses the model from scratch,
# checks its structural invariants, evaluates every gridded table, and
# runs each checkData staticShot against its stated tolerance.
#
# It shares no code with the assembler deliberately -- an evaluator built
# on the writer's own helpers would agree with it by construction.
#
# The test is registered only when the model has actually been assembled;
# a checkout that has never run build-daveml.py should not fail its test
# suite over an artifact it was never asked to produce.

find_package(Python3 COMPONENTS Interpreter QUIET)

if(NOT Python3_Interpreter_FOUND)
  message(STATUS "DAVE-ML: no Python 3 interpreter; skipping TestDaveMLModel")
  return()
endif()

if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/AetherionFlightModel.dml")
  message(STATUS "DAVE-ML: model not assembled yet; skipping TestDaveMLModel "
                 "(run models/build-daveml.py)")
  return()
endif()

add_test(NAME TestDaveMLModel
         COMMAND ${Python3_EXECUTABLE}
                 "${CMAKE_CURRENT_SOURCE_DIR}/verify-daveml.py"
                 "${CMAKE_CURRENT_SOURCE_DIR}/AetherionFlightModel.dml")
