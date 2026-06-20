add_executable(test_sinc_quadrature test_sinc_quadrature.cpp)

target_link_libraries(test_sinc_quadrature PRIVATE sophus::sophus)

# Keep physical guardrail asserts active under Release builds.
target_compile_options(test_sinc_quadrature PRIVATE -UNDEBUG)

add_test(NAME SincQuadratureTest COMMAND test_sinc_quadrature)
