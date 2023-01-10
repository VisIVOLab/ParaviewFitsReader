paraview_add_plugin(FitsReader
  REQUIRED_ON_SERVER
  VERSION "1.0"
  MODULES FitsReader::vtkFitsReader
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/vtkFitsReader/vtk.module"
)