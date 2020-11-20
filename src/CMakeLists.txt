
set(SRC 
    Body.cpp
    Joint.cpp
    Model.cpp
    ModelConstructor.cpp
)

set(HEADERS 
    ${PROJECT_NAME}/Core
    ${PROJECT_NAME}/API.hpp
    ${PROJECT_NAME}/Body.hpp
    # ${PROJECT_NAME}/FD.hpp
    ${PROJECT_NAME}/FK.hpp
    ${PROJECT_NAME}/ID.hpp
    ${PROJECT_NAME}/Joint.hpp
    ${PROJECT_NAME}/Model.hpp
    ${PROJECT_NAME}/ModelConfig.hpp
    ${PROJECT_NAME}/ModelConstructor.hpp
    ${PROJECT_NAME}/typedefs.hpp
)

add_library(${PROJECT_NAME} SHARED ${SRC})
target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../include> $<INSTALL_INTERFACE:include>)
target_link_libraries(${PROJECT_NAME} PUBLIC coma::coma)
set_target_properties(${PROJECT_NAME} PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR} VERSION ${PROJECT_VERSION})
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DCDM_EXPORTS")

install(
    TARGETS cdm
    EXPORT "${TARGETS_EXPORT_NAME}"
    LIBRARY DESTINATION "lib"
    ARCHIVE DESTINATION "lib"
    RUNTIME DESTINATION "bin"
)
install(FILES ${HEADERS} DESTINATION ${INCLUDE_INSTALL_DESTINATION})