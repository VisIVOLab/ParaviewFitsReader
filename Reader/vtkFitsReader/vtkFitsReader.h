// .NAME vtkFitsReader - read structured points from FITS file.
// .SECTION Description
// vtkFitsReader is a source object that reads FITS data files
// .SECTION Caveats
// Uses CFITSIO v4.1.0 (http://heasarc.gsfc.nasa.gov/docs/software/fitsio)

#ifndef __vtkFitsReader_h
#define __vtkFitsReader_h

#include <vtkMPIImageReader.h>
#include <vtkNew.h>

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
    const char *GetFileExtensions() override
    {
        return ".fits";
    }

    /**
     * Return a descriptive name for the file format that might be useful in a
     * GUI.
     */
    const char *GetDescriptiveName() override
    {
        return "FITS";
    }

    vtkGetMacro(ReadSubExtent, bool);
    vtkSetMacro(ReadSubExtent, bool);

    vtkGetVector6Macro(SubExtent, int);
    vtkSetVector6Macro(SubExtent, int);

    vtkGetMacro(ScaleFactor, int);
    vtkSetMacro(ScaleFactor, int);

  protected:
    vtkFitsReader();
    ~vtkFitsReader() override;

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

    /**
     * @brief This property specifies if the reader must read a subset of the
     * data.
     *
     */
    bool ReadSubExtent;

    /**
     * @brief This property specifies the sub-extent to read. It is ignored if
     * ReadSubExtent is disabled.
     *
     */
    int SubExtent[6];

    /**
     * @brief This property can be used to read only every inc-th pixel along the
     * dimensions of the image.
     *
     */
    int ScaleFactor;
};
#endif
