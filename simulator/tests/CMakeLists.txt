add_executable(test_sinc_quadrature test_sinc_quadrature.cpp)

target_link_libraries(test_sinc_quadrature PRIVATE sophus::sophus)

target_compile_options(test_sinc_quadrature PRIVATE -UNDEBUG)

add_test(NAME SincQuadratureTest COMMAND test_sinc_quadrature)

add_executable(test_fermi_dirac test_fermi_dirac.cpp)

target_link_libraries(test_fermi_dirac PRIVATE sophus::sophus)

target_compile_options(test_fermi_dirac PRIVATE -UNDEBUG)

add_test(NAME FermiDiracTest COMMAND test_fermi_dirac)
