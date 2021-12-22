# The pqReader.xml file contains xml defining readers with their file
# extensions and descriptions.  The Server Manager XML file defines the
# properties of the reader.  This example uses the already existing PNG reader
# but gives it another identity in the Server Manager.  Adding your own reader
# to ParaView would mean the SourceProxy class points to your class.

find_library (CFITSIO_LIB NAMES cfitsio PATHS ${CFITSIO_LIB_PATH})
include_directories(${CFITSIO_INC_PATH})
include_directories(${BOOST_INC_PATH})

paraview_add_plugin(FitsReader
  REQUIRED_ON_SERVER
  VERSION "1.0"
  MODULES FitsReader::vtkFitsReader
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/vtkFitsReader/vtk.module"
)

