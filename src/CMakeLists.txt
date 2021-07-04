
set(SRC 
    Body.cpp
    Joint.cpp
    Model.cpp
    ModelConstructor.cpp
)

set(HEADERS 
    ../include/${PROJECT_NAME}/API.hpp
    ../include/${PROJECT_NAME}/BasicJacobian.hpp
    ../include/${PROJECT_NAME}/Body.hpp
    ../include/${PROJECT_NAME}/FD.hpp
    ../include/${PROJECT_NAME}/FK.hpp
    ../include/${PROJECT_NAME}/ID.hpp
    ../include/${PROJECT_NAME}/Joint.hpp
    ../include/${PROJECT_NAME}/JointForceJacobian.hpp
    ../include/${PROJECT_NAME}/JointMomentumJacobian.hpp
    ../include/${PROJECT_NAME}/LinkForceJacobian.hpp
    ../include/${PROJECT_NAME}/LinkMomentumJacobian.hpp
    # ${PROJECT_NAME}/LowerBlockTriangularMatrix.hpp
    ../include/${PROJECT_NAME}/math_utility.hpp
    ../include/${PROJECT_NAME}/Model.hpp
    ../include/${PROJECT_NAME}/ModelConfig.hpp
    ../include/${PROJECT_NAME}/ModelConstructor.hpp
    ../include/${PROJECT_NAME}/typedefs.hpp
)

add_library(${PROJECT_NAME} SHARED ${SRC})
target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../include> $<INSTALL_INTERFACE:include>)
target_link_libraries(${PROJECT_NAME} PUBLIC coma::coma)
set_target_properties(${PROJECT_NAME} PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR} VERSION ${PROJECT_VERSION})
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DCDM_EXPORTS")

install(
    TARGETS ${PROJECT_NAME}
    EXPORT "${TARGETS_EXPORT_NAME}"
    LIBRARY DESTINATION "lib"
    ARCHIVE DESTINATION "lib"
    RUNTIME DESTINATION "bin"
)
install(FILES ${HEADERS} DESTINATION ${INCLUDE_INSTALL_DESTINATION})