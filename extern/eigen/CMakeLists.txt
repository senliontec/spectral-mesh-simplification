add_library(eigen INTERFACE)
target_include_directories(eigen INTERFACE ${CMAKE_CURRENT_LIST_DIR}/eigen)
if(MSVC)
	target_sources(eigen INTERFACE ${CMAKE_CURRENT_LIST_DIR}/eigen/debug/msvc/eigen.natvis)
endif()
