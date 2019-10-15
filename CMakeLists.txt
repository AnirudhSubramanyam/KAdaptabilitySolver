# CHECK CMAKE VERSION
cmake_minimum_required(VERSION 2.8)

# SPECIFY PROJECT NAME
project(K_ADAPTABILITY)

# SPECIFY SOURCE FILES (.CPP, .C)
set(SOURCE_FILES
    ./src/uncertainty.cpp
    ./src/problemInfo.cpp
    ./src/problemInfo_knp.cpp
    ./src/problemInfo_spp.cpp
    ./src/problemInfo_psp.cpp
    ./src/indexingTools.cpp
    ./src/robustSolver.cpp
    ./src/test.cpp
)

# SPECIFY CUSTOM HEADER FILE DIRECTORIES (OR SET "DEFAULT" IF NONE)
set(HDRDIR ./inc)

# SPECIFY NAME OF EXECUTABLE
set(EXEC kadaptability)

# SET COMPILER DIRECTORY (OR SET "DEFAULT")
set(COMPILER_DIR /usr/local/bin/g++)

# SET CPLEX DIRECTORY
set(CPLEX_DIR /opt/ibm/ILOG/CPLEX_Studio127)

# CHECK OS
if(APPLE)
  set(mySystem "osx")
elseif(UNIX)
  set(mySystem "linux")
else()
  message(FATAL_ERROR "ONLY UNIX IS SUPPORTED AT THIS TIME")
endif()

# SET COMPILER DIRECTORY
if(NOT "${COMPILER_DIR}" STREQUAL "DEFAULT")
  set(CMAKE_CXX_COMPILER ${COMPILER_DIR})
endif()

# SET COMPILER FLAGS
set(COMP_FLAGS_DEBUG " -g -std=c++14 -m64 -O0 -DIL_STD -pthread -D_GLIBCXX_DEBUG -Wall -Wextra -fPIC -fno-strict-aliasing -fexceptions")
set(COMP_FLAGS_RELEASE " -std=c++14 -m64 -O2 -DIL_STD -pthread -DNDEBUG -Wall -Wextra -fPIC -fno-strict-aliasing -fexceptions")
set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS} ${COMP_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${COMP_FLAGS_RELEASE}")

# SET INCLUDE DIRECTORIES (.h .hpp)
include_directories(${EXEC} ${CPLEX_DIR}/cplex/include ${CPLEX_DIR}/concert/include)

# INCLUDE THE CUSTOM INC DIRECTORIES IF DEFINED
if(NOT "${HDRDIR}" STREQUAL "DEFAULT")
  include_directories(${EXEX} ${HDRDIR})
endif()

# SET LINK DIRECTORIES
link_directories(${CPLEX_DIR}/cplex/lib/x86-64_${mySystem}/static_pic ${CPLEX_DIR}/concert/lib/x86-64_${mySystem}/static_pic)

# CREATE THE TARGET (executable or library)
add_executable(${EXEC} ${SOURCE_FILES})

# LINK THE DESIRED LIBRARIES
target_link_libraries(${EXEC} concert ilocplex cplex m pthread)
