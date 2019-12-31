cmake_minimum_required(VERSION 3.11)
project(SpectralCollapsing)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "cmake")

add_subdirectory(extern/eigen)
add_subdirectory(extern/spectra)
add_subdirectory(extern/libigl)

add_executable(spectral-collapsing)
	file(GLOB_RECURSE SOURCES source/*.h source/*.cpp)
	source_group(TREE ${CMAKE_CURRENT_LIST_DIR} FILES ${SOURCES})
	target_sources(spectral-collapsing PUBLIC ${SOURCES})
	target_compile_features(spectral-collapsing PUBLIC cxx_std_17)
	target_link_libraries(spectral-collapsing
		PUBLIC eigen
		PUBLIC spectra
		PUBLIC libigl)
	set_target_properties(spectral-collapsing PROPERTIES
		VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin
		RUNTIME_OUTPUT_DIRECTORY      $<1:${CMAKE_CURRENT_LIST_DIR}/bin>)

	if(MSVC)
		target_compile_options(spectral-collapsing PRIVATE /W4)
	else()
		target_compile_options(spectral-collapsing PRIVATE -Wall -Wextra -pedantic)
	endif()
	set_property(GLOBAL PROPERTY VS_STARTUP_PROJECT spectral-collapsing)
