# The pqReader.xml file contains xml defining readers with their file
# extensions and descriptions.  The Server Manager XML file defines the
# properties of the reader.  This example uses the already existing PNG reader
# but gives it another identity in the Server Manager.  Adding your own reader
# to ParaView would mean the SourceProxy class points to your class.



SET(classes
  vtkFitsReader)

    
vtk_module_add_module(FitsReader::vtkFitsReader
  FORCE_STATIC # Using FORCE_STATIC build the vtk module statically into the plugin library, to avoid confusion when loading                      
  CLASSES ${classes})

paraview_add_server_manager_xmls(
  XMLS  FitsReaderServerManager.xml)


target_link_libraries(vtkFitsReader PRIVATE ${CFITSIO_LIB})
