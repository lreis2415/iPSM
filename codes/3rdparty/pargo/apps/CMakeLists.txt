# GPRO source code
SET(GPRO_SRCS ${GPRO_DIR}/application.cpp)
# Adds include directories in such a way that directories from the TauDEM source tree go first.
geo_include_directories(${GPRO_DIR} ${GDAL_INCLUDE_DIR} ${MPI_INCLUDE_PATH})
STATUS("  input apps")
ADD_SUBDIRECTORY(./hydrology)
ADD_SUBDIRECTORY(./morphology)
#ADD_SUBDIRECTORY(./spatial) # TODO: 'class GPRO::RasterLayer<double>' has no member named 'rowWriteFile'

