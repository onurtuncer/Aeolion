# CompilerWarnings.cmake
#
# aeolion_enable_warnings(<interface-target>) -- turns on the project's
# warning set on a header-only INTERFACE library, honouring the
# AEOLION_WARNINGS option. Callers of that library inherit the flags.

function(aeolion_enable_warnings target)
    if(NOT AEOLION_WARNINGS)
        return()
    endif()
    if(MSVC)
        target_compile_options(${target} INTERFACE /W4)
    else()
        target_compile_options(${target} INTERFACE -Wall -Wextra)
    endif()
endfunction()
