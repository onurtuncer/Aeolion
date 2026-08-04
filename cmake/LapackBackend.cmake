# LapackBackend.cmake
#
# aeolion_link_lapack(<interface-target>) -- finds a LAPACK provider and
# links it into the given INTERFACE library.
#
# The dense linear solve calls LAPACK's dgetrf/dgetrs through their Fortran
# ABI directly -- no LAPACKE C header needed, so a LAPACK *library* is the
# only hard dependency. CMake's FindLAPACK locates OpenBLAS, reference
# LAPACK, MKL, Accelerate, etc. and pulls in the matching BLAS.
function(aeolion_link_lapack target)
  find_package(LAPACK REQUIRED)
  target_link_libraries(${target} INTERFACE LAPACK::LAPACK)
  message(STATUS "aeolion: using LAPACK (${LAPACK_LIBRARIES})")
endfunction()
