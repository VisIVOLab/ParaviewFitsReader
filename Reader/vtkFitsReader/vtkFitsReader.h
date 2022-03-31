// .NAME vtkFitsReader - read structured points from FITS file.
// .SECTION Description
// vtkFitsReader is a source object that reads FITS data files
// .SECTION Caveats
// Uses CFITSIO v2.0 (http://heasarc.gsfc.nasa.gov/docs/software/fitsio)

#ifndef __vtkFitsReader_h
#define __vtkFitsReader_h

#include "vtkIOImageModule.h" // For export macro
#include "vtkMPIImageReader.h"
#include "vtkStructuredPoints.h"
#include "vtkFloatArray.h"
#include "vtkMultiProcessController.h"

extern "C"
{
    #include "fitsio.h"
}

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

    // Description:
    /// Report the status of the reading process.
    /// If this is different than zero, there have been some error
    /// parsing the complete header information.
    vtkGetMacro(ReadStatus, int);

    /// Point data field type
    vtkSetMacro(PointDataType, int);
    vtkGetMacro(PointDataType, int);

    ///
    /// Set the data type: int, float....
    vtkSetMacro(DataType, int);
    vtkGetMacro(DataType, int);

    ///
    /// Number of components
    vtkSetMacro(NumberOfComponents, int);
    vtkGetMacro(NumberOfComponents, int);

protected:
    vtkFitsReader();
    ~vtkFitsReader() override;

    fitsfile *fptr;
    int ReadStatus;
    int PointDataType;
    int DataType;
    int NumberOfComponents;
    long naxes[3];

    int RequestInformation(vtkInformation*, vtkInformationVector**,  vtkInformationVector* outVec) override;
  
    int RequestData(vtkInformation*, vtkInformationVector**,vtkInformationVector* outVec)override;

private:
    vtkFitsReader(const vtkFitsReader &) = delete;
    void operator=(const vtkFitsReader &) = delete;
    void printerror(int status); // from fitsio distribution

    int dataExtent[6] = {0};

    vtkMultiProcessController* controller;
};
#endif
