SET(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty)
SET(GPRO_DIR ${THIRD_PARTY_DIR}/pargo/gpro)
#SET(BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/block-base)
SET(BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/parallel-base)
SET(COMMON_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/common-base)

# TauDEM Parallel IO source code
SET(TAUDEM_SRCS ${THIRD_PARTY_DIR}/TauDEM_IO/CommonLib.cpp ${THIRD_PARTY_DIR}/TauDEM_IO/tiffIO.cpp ${THIRD_PARTY_DIR}/TauDEM_IO/createpart.cpp)
# Base functions for iPSM Project
file(GLOB_RECURSE BASE_SRCS ${BASE_DIR}/*.c ${BASE_DIR}/*.cpp ${BASE_DIR}/*.h ${COMMON_BASE_DIR}/*.cpp ${COMMON_BASE_DIR}/*.h ${GPRO_DIR}/*.cpp ${GPRO_DIR}/*.h)
#file(GLOB_RECURSE BASE_SRCS ${BASE_DIR}/*.c ${BASE_DIR}/*.cpp ${BASE_DIR}/*.h ${COMMON_BASE_DIR}/*.cpp ${COMMON_BASE_DIR}/*.h))
# TinyXML source code for integrating inference knowledge stored in xml format 2021/3/11 Zhaofh
file(GLOB XML_SRCS ${THIRD_PARTY_DIR}/tixml/*.cpp ${THIRD_PARTY_DIR}/tixml/*.h)
# SET(ALL_BASE_SRCS ${TAUDEM_SRCS} ${BASE_SRCS} ${UTIL_SRCS})
# Integrating Pargo
SET(ALL_BASE_SRCS ${BASE_SRCS})

# Adds include directories in such a way that directories from the TauDEM source tree go first.
geo_include_directories(${BASE_DIR} ${THIRD_PARTY_DIR} ${THIRD_PARTY_DIR}/tixml ${GDAL_INCLUDE_DIR} ${MPI_INCLUDE_PATH} ${COMMON_BASE_DIR} ${THIRD_PARTY_DIR}/pargo/gpro)

# ADD_SUBDIRECTORY(./utils) # No need to build? 2019/08/06 ZHULJ
ADD_SUBDIRECTORY(./apps)
#ADD_SUBDIRECTORY(./3rdparty)

