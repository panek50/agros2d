project(plugins)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Add QT.
SET(QT_USE_QTOPENGL TRUE)
SET(QT_USE_QTUITOOLS TRUE)
SET(QT_USE_QTNETWORK TRUE)
SET(QT_USE_QTOPENGL TRUE)
SET(QT_USE_QTSQL TRUE)
SET(QT_USE_QTXML TRUE)
SET(QT_USE_QTSVG TRUE)
SET(QT_USE_QTTEST TRUE)
SET(QT_USE_QTDBUS TRUE)
SET(QT_USE_QTSCRIPT TRUE)
SET(QT_USE_QTWEBKIT TRUE)
SET(QT_USE_QTXMLPATTERNS TRUE)
SET(QT_USE_PHONON TRUE)
find_package(Qt4 REQUIRED)

INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})
ADD_DEFINITIONS(-DQT_PLUGIN)
ADD_DEFINITIONS(-DQT_SHARED)
ADD_DEFINITIONS(-DQT_DLL)

# TODO is it correct place?
SET(PARALUTION_LIBRARY agros2d_3dparty_paralution)
set(WITH_PARALUTION YES)
SET(PARALUTION_LIBRARY ${PARALUTION_LIBRARY})
set(PARALUTION_INCLUDE_DIR ${CMAKE_HOME_DIRECTORY}/../3rdparty/paralution/src)

# Set global compiler parameters.
find_package(OpenMP REQUIRED)
IF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  INCLUDE_DIRECTORIES(/usr/include/google)
  IF(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  INCLUDE_DIRECTORIES(omp)
  ENDIF(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

  INCLUDE_DIRECTORIES(/usr/include)
  INCLUDE_DIRECTORIES(/usr/include/python2.7)
  
  EXECUTE_PROCESS(COMMAND "python" "-c" "import distutils.sysconfig; print distutils.sysconfig.get_config_var('LOCALMODLIBS')" OUTPUT_VARIABLE PYTHON_MODLIBS)
  EXECUTE_PROCESS(COMMAND "python" "-c" "from distutils import sysconfig; print '-lpython'+sysconfig.get_config_var('VERSION')" OUTPUT_VARIABLE PYTHON_LIB)
  STRING(STRIP ${PYTHON_MODLIBS} PYTHON_MODLIBS)
  STRING(STRIP ${PYTHON_LIB} PYTHON_LIB)
  SET(CMAKE_SHARED_LINKING_FLAGS "${CMAKE_SHARED_LINKING_FLAGS} ${PYTHON_MODLIBS} ${PYTHON_LIB}")
  SET(CMAKE_MODULE_LINKING_FLAGS "${CMAKE_MODULE_LINKING_FLAGS} ${PYTHON_MODLIBS} ${PYTHON_LIB}")
  SET(CMAKE_EXE_LINKING_FLAGS "${CMAKE_EXE_LINKING_FLAGS} ${PYTHON_MODLIBS} ${PYTHON_LIB}")
ENDIF()

IF(MSVC)
  INCLUDE_DIRECTORIES(c:/hpfem/hermes/dependencies/include)
  INCLUDE_DIRECTORIES(d:/hpfem/hermes/dependencies/include)
  INCLUDE_DIRECTORIES(c:/Python27/include)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /openmp /Zc:wchar_t")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} /NODEFAULTLIB:libcmtd /NODEFAULTLIB:libcmt")
ENDIF(MSVC)

SET(CMAKE_AGROS_DIRECTORY "${CMAKE_HOME_DIRECTORY}/../")
# Include OUR header files location
include(${CMAKE_AGROS_DIRECTORY}/IncludeSubdirs.cmake)

INCLUDE(${CMAKE_AGROS_DIRECTORY}/hermes/CMake.vars)
SET(CMAKE_MODULE_PATH ${CMAKE_AGROS_DIRECTORY}/hermes/cmake)

# Look for UMFPACK, etc.
# Solvers.
IF(WITH_UMFPACK)
  FIND_PACKAGE(UMFPACK REQUIRED)
  INCLUDE_DIRECTORIES(${UMFPACK_INCLUDE_DIRS})
ENDIF(WITH_UMFPACK)

IF(WITH_PARALUTION)
  FIND_PACKAGE(PARALUTION REQUIRED)
  INCLUDE_DIRECTORIES(${PARALUTION_INCLUDE_DIR})
ENDIF()

IF(WITH_MUMPS)
  FIND_PACKAGE(MUMPS REQUIRED)
  INCLUDE_DIRECTORIES(${PARALUTION_INCLUDE_DIR})
ENDIF()

{{#SOURCE}}
ADD_SUBDIRECTORY({{ID}})
{{/SOURCE}}

