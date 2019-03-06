cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(OfflineMapping_basic)

set (CMAKE_CXX_STANDARD 11)
# ===================== dependencies ==========================
# PCL
find_package(PCL 1.2 REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

# local
INCLUDE_DIRECTORIES(./)
# file(GLOB_RECURSE HEADERS  ${CMAKE_CURRENT_LIST_DIR}  "*.h" "*.hpp")


# ===================== executables ==========================
AUX_SOURCE_DIRECTORY(. DIR_SRCS)
add_executable (${PROJECT_NAME} ${DIR_SRCS})
target_link_libraries (${PROJECT_NAME} ${PCL_LIBRARIES})
