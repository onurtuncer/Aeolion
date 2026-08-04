# Aeolion Solver: the VLM core. Header-only -- solver/include is added to
# the `aeolion` INTERFACE target's include path directly from the root
# CMakeLists.txt, the same way include/ itself is, since there is no
# compiled part here to build a separate target for (contrast
# panelbuilder/, which does have one). This file's only job is the solver's
# own test folder.

if(AEOLION_BUILD_TESTS)
    add_subdirectory(tests)
endif()
