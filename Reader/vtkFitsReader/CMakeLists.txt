find_package(Boost REQUIRED)
find_package(CFITSIO REQUIRED)

SET(classes
  vtkFitsReader)

vtk_module_add_module(FitsReader::vtkFitsReader
  FORCE_STATIC # Using FORCE_STATIC build the vtk module statically into the plugin library, to avoid confusion when loading                      
  CLASSES ${classes})

paraview_add_server_manager_xmls(
  XMLS  FitsReaderServerManager.xml)

target_include_directories(vtkFitsReader PRIVATE ${Boost_INCLUDE_DIRS})

target_link_libraries(vtkFitsReader PRIVATE CFITSIO::CFITSIO)
