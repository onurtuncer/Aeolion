# CompilerWarnings.cmake
#
# aeolion_enable_warnings(<interface-target>) -- turns on the project's
# warning set on the header-only INTERFACE library, honouring the
# AEOLION_WARNINGS option. Everything that links it inherits the flags,
# including Aeolion's own compiled libraries, so their sources are held to
# the same standard as consumer code.

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
