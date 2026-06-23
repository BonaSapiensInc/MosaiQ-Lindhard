add_executable(test_sinc_quadrature test_sinc_quadrature.cpp)

target_link_libraries(test_sinc_quadrature PRIVATE mosaiq::mosaiq)

target_compile_options(test_sinc_quadrature PRIVATE -UNDEBUG)

add_test(NAME SincQuadratureTest COMMAND test_sinc_quadrature)

add_executable(test_fermi_dirac test_fermi_dirac.cpp)

target_link_libraries(test_fermi_dirac PRIVATE mosaiq::mosaiq)

target_compile_options(test_fermi_dirac PRIVATE -UNDEBUG)

add_test(NAME FermiDiracTest COMMAND test_fermi_dirac)

add_executable(test_lindhard test_lindhard.cpp)

target_link_libraries(test_lindhard PRIVATE mosaiq::mosaiq)

target_compile_options(test_lindhard PRIVATE -UNDEBUG)

add_test(NAME LindhardTest COMMAND test_lindhard)

add_executable(test_rpa test_rpa.cpp)

target_link_libraries(test_rpa PRIVATE mosaiq::mosaiq)

target_compile_options(test_rpa PRIVATE -UNDEBUG)

add_test(NAME RPATest COMMAND test_rpa)

add_executable(test_plasmon_pole_extractor test_plasmon_pole_extractor.cpp)

target_link_libraries(test_plasmon_pole_extractor PRIVATE mosaiq::mosaiq)

target_compile_options(test_plasmon_pole_extractor PRIVATE -UNDEBUG)

add_test(NAME PlasmonPoleExtractorTest COMMAND test_plasmon_pole_extractor)

add_executable(test_sum_rules test_sum_rules.cpp)

target_link_libraries(test_sum_rules PRIVATE mosaiq::mosaiq)

target_compile_options(test_sum_rules PRIVATE -UNDEBUG)

add_test(NAME SumRulesTest COMMAND test_sum_rules)

add_executable(test_t0_limit test_t0_limit.cpp)

target_link_libraries(test_t0_limit PRIVATE mosaiq::mosaiq)

target_compile_options(test_t0_limit PRIVATE -UNDEBUG)

add_test(
    NAME T0LimitTest
    COMMAND test_t0_limit
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/..
)
