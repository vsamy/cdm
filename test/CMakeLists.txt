enable_testing()

find_package(RBDyn REQUIRED)

macro(addTest testName)
    add_executable(${testName} ${testName}.cpp)
    target_compile_definitions(${testName} PUBLIC DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS)
    # if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    #     target_compile_options(${testName} PUBLIC /W4)
    # else()
    #     target_compile_options(${testName} PUBLIC -Wall)
    # endif()
    target_link_libraries(${testName} PUBLIC ${PROJECT_NAME} Eigen3::Eigen RBDyn::RBDyn)
    # Adding a project configuration file (for MSVC only)
    generate_msvc_dot_user_file(${testName})

    add_test(${testName}Unit ${testName})
endmacro(addTest)

addTest(model_test)
addTest(forward_kinematics_test)
addTest(inverse_dynamics_test)
addTest(forward_dynamics_test)
addTest(jacobians_test)
