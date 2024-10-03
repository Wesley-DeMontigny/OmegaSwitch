file(GLOB ncl_INC "*.h")
file(GLOB ncl_SRC "*.cpp")

# Adjust to include the current source directory and its parent
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)

add_library(ncl_shared SHARED ${ncl_SRC})
add_library(ncl_static STATIC ${ncl_SRC})

target_include_directories(ncl_shared PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_include_directories(ncl_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

set_target_properties(ncl_shared PROPERTIES OUTPUT_NAME ncl)
set_target_properties(ncl_static PROPERTIES OUTPUT_NAME ncl)

install(TARGETS ncl_shared DESTINATION "lib/ncl")
install(TARGETS ncl_static DESTINATION "lib/ncl")
install(FILES ${ncl_INC} DESTINATION "include/ncl")