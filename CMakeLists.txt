cmake_minimum_required(VERSION 3.10)

project(xpbd)

file(GLOB ELE_FILES "${CMAKE_SOURCE_DIR}/data/*.ele")

foreach(ELE ${ELE_FILES})
	file(COPY ${ELE} DESTINATION ${CMAKE_BINARY_DIR})
endforeach()

file(GLOB FACE_FILES "${CMAKE_SOURCE_DIR}/data/*.face")

foreach(FACE ${FACE_FILES})
	file(COPY ${FACE} DESTINATION ${CMAKE_BINARY_DIR})
endforeach()

file(GLOB NODE_FILES "${CMAKE_SOURCE_DIR}/data/*.node")

foreach(NODE ${NODE_FILES})
	file(COPY ${NODE} DESTINATION ${CMAKE_BINARY_DIR})
endforeach()

file(GLOB SMESH_FILES "${CMAKE_SOURCE_DIR}/data/*.smesh")

foreach(SMESH ${SMESH_FILES})
	file(COPY ${NODE} DESTINATION ${CMAKE_BINARY_DIR})
endforeach()



file(GLOB Shader "${CMAKE_CURRENT_SOURCE_DIR}/shader")

file(COPY ${Shader} DESTINATION ${CMAKE_BINARY_DIR})


set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 11)

include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/include)

link_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/lib)

add_executable(xpbd_cloth "xpbd_cloth.cpp" )

add_executable(xpbd_ball "xpbd_ball.cpp"   )

add_executable(xpbd_spring "xpbd_spring.cpp"  )

add_executable(xpbd_soft "xpbd_soft.cpp" )

add_executable(xpbd_col "xpbd_collision.cpp" )

add_definitions(-DFREEGLUT_STATIC)
target_link_libraries(xpbd_cloth freeglut_static.lib)
target_link_libraries(xpbd_ball glfw3.lib)
target_link_libraries(xpbd_spring glfw3.lib)
target_link_libraries(xpbd_soft glfw3.lib)
target_link_libraries(xpbd_col glfw3.lib)