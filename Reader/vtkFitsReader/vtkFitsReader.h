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
    enum imageType { EMPTY, FITS2DIMAGE, FITS3DIMAGE };
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
        return ".fits .fit .fts";
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

    vtkGetMacro(AutoScale, bool);
    vtkSetMacro(AutoScale, bool);

    vtkGetMacro(CubeMaxSize, int);
    vtkSetMacro(CubeMaxSize, int);

    vtkGetMacro(ScaleFactor, int);
    vtkSetMacro(ScaleFactor, int);

    vtkGetMacro(ImgType, int);

    vtkGetMacro(ReadAsPVSlice, bool);
    vtkSetMacro(ReadAsPVSlice, bool);

    void SetStartPoint(int, int);
    void SetEndPoint(int, int);
    void SetZBounds(int newZMin, int newZMax);
  protected:
    vtkFitsReader();
    ~vtkFitsReader() override;

    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *outVec) override;
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *outVec) override;
    int FillOutputPortInformation(int port, vtkInformation *info) override;

  private:
    vtkFitsReader(const vtkFitsReader &) = delete;
    vtkFitsReader &operator=(const vtkFitsReader &) = delete;

    int ReadPVSliceData(int ProcID, vtkInformationVector* outVec);

    /**
     * @brief FITS Header
     *
     */
    vtkNew<vtkTable> table;

    /**
     * @brief This property specifies if the FITS file is an image (2D) or a cube (3D).
     *
     */
    imageType ImgType;

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
     * @brief This property specifies whether to use a ScaleFactor by default.
     *
     */
    bool AutoScale;

    /**
     * @brief This property can be used along with AutoScale to use at most MaxCubeSize (MB)
     *  for reading the cube.
     *
     */
    int CubeMaxSize;

    /**
     * @brief This property can be used to read only every inc-th pixel along the
     * dimensions of the image.
     *
     */
    int ScaleFactor;

    /**
     * @brief This property can be used to read a position-velocity slice from the cube
     * instead of reading the entire cube. The start and end points must be set as well.
     */
    bool ReadAsPVSlice = false;

    /**
     * @brief This property specifies the start coordinate for the PV slice.
     */
    std::pair<int, int> PVStart;

    /**
     * @brief This property specifies the end coordinate for the PV slice.
     */
    std::pair<int, int> PVEnd;

    /**
     * @brief This property specifies the bounds in the z-axis for the PV slice.
     */
    std::pair<int, int> PVZSubset;

    /**
     * @brief Convert2DIndexToLinear
     *        Utility function to convert a 2D index to a linear array index, assuming row-major storage.
     * @param row The row component of the 2D index
     * @param col The column component of the 2D index
     * @param rowWidth The width of the rows in the storage
     * @return The index in the linear array corresponding to (row, col).
     */
    int Convert2DIndexToLinear(const int row, const int col, const int rowWidth){
        return rowWidth * row + col;
    }
};

/**
 * Linear interpolation
 */
template<typename T>
static const T interpolate(const T& fromPointA, const T& toPointB, float factor){
    return ((1.0 - factor) * fromPointA) + (factor * toPointB);
}

/**
 * Bilinear interpolation
 */
template<typename T>
static const T interpolate(const T& fromPointA1, const T& fromPointA2, const T& toPointB1, const T& toPointB2, float factor1, float factor2){
    return interpolate(interpolate(fromPointA1, fromPointA2, factor1), interpolate(toPointB1, toPointB2, factor1), factor2);
}
#endif
