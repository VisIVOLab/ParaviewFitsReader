// .NAME vtkFitsReader - read structured points from FITS file.
// .SECTION Description
// vtkFitsReader is a source object that reads FITS data files
// .SECTION Caveats
// Uses CFITSIO v4.1.0 (http://heasarc.gsfc.nasa.gov/docs/software/fitsio)

#ifndef __vtkFitsReader_h
#define __vtkFitsReader_h

#include <vtkMPIImageReader.h>
#include <vtkNew.h>
#include <vtkTable.h>

#include <map>
#include <string>
#include <vector>

class vtkTable;

class VTK_EXPORT vtkFitsReader : public vtkMPIImageReader
{
public:
    static vtkFitsReader *New();
    vtkTypeMacro(vtkFitsReader, vtkMPIImageReader);
    void PrintSelf(ostream &os, vtkIndent indent) override;

    int CanReadFile(VTK_FILEPATH const char *fname) override;

    /**
     * Get the file extensions for this format.
     * Returns a string with a space separated list of extensions in
     * the format .extension
     */
    const char *GetFileExtensions() override { return ".fits"; }

    /**
     * Return a descriptive name for the file format that might be useful in a GUI.
     */
    const char *GetDescriptiveName() override { return "FITS"; }

    /// Point data field type
    vtkSetMacro(PointDataType, int);
    vtkGetMacro(PointDataType, int);

    /// Set the data type: int, float....
    vtkSetMacro(DataType, int);
    vtkGetMacro(DataType, int);

    /// Number of components
    vtkSetMacro(NumberOfComponents, int);
    vtkGetMacro(NumberOfComponents, int);

protected:
    vtkFitsReader();
    ~vtkFitsReader() override;

    int PointDataType;
    int DataType;
    int NumberOfComponents;

    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *outVec) override;
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *outVec) override;
    int FillOutputPortInformation(int port, vtkInformation *info) override;

private:
    vtkFitsReader(const vtkFitsReader &) = delete;
    vtkFitsReader &operator=(const vtkFitsReader &) = delete;

    /**
     * @brief FITS Header
     *
     */
    vtkNew<vtkTable> table;

    /**
     * @brief   Read the header and store the key-value pairs in g.
     *          This function ignores HISTORY, COMMENT and empty keywords.
     *
     * @return  0 on success, greater than 0 otherwise.
     */
    int ReadFITSHeader();
};
#endif
