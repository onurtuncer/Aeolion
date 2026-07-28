# The core header-only library's (Math/Lattice/Geometry/DragEstimate) own
# tests. Every other module's tests live beside the module they test --
# solver/tests/, panelbuilder/tests/ -- registered through the shared
# aeolion_add_test() in cmake/AeolionTest.cmake.
aeolion_add_test(TestHandoffContract)
aeolion_add_test(TestPropeller)
