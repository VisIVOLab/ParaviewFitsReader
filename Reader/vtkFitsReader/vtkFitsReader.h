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

<<<<<<< Updated upstream
  public:
    enum imageType { EMPTY, FITS2DIMAGE, FITS3DIMAGE };
    static vtkFitsReader *New();
    vtkTypeMacro(vtkFitsReader, vtkMPIImageReader);
    void PrintSelf(ostream &os, vtkIndent indent) override;
=======
    public:
        enum imageType { EMPTY, FITS2DIMAGE, FITS3DIMAGE };
        enum readerType { RAW, MOMENTMAP, PVSLICE };
        static vtkFitsReader *New();
        vtkTypeMacro(vtkFitsReader, vtkMPIImageReader);
        void PrintSelf(ostream &os, vtkIndent indent) override;
>>>>>>> Stashed changes

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

<<<<<<< Updated upstream
    vtkGetMacro(ReadSubExtent, bool);
    vtkSetMacro(ReadSubExtent, bool);
=======
        vtkGetMacro(ReadType, int);
        /**
         * @brief SetReadType
         * Function to set the ReadType property.
         * @param type
         * 0 is reading the file raw.
         * 1 is reading a moment map (see MomentOrder).
         * 2 is reading a position-velocity slice.
         */
        void SetReadType(int type) {this->ReadType = (readerType) type;};
>>>>>>> Stashed changes

    vtkGetVector6Macro(SubExtent, int);
    vtkSetVector6Macro(SubExtent, int);

<<<<<<< Updated upstream
    vtkGetMacro(AutoScale, bool);
    vtkSetMacro(AutoScale, bool);
=======
        void SetStartPoint(int x, int y);
        void SetEndPoint(int x, int y);
        void SetZBounds(int z1, int z2);

        vtkGetMacro(ReadSubExtent, bool);
        vtkSetMacro(ReadSubExtent, bool);
>>>>>>> Stashed changes

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

<<<<<<< Updated upstream
    /**
     * @brief This property specifies if the FITS file is an image (2D) or a cube (3D).
     *
     */
    imageType ImgType;
=======
        /**
         * @brief ReadType This property specifies what algorithm the reader is using to read the FITS file.
         * 0 is reading the file raw.
         * 1 is reading a moment map (see MomentOrder).
         * 2 is reading a position-velocity slice.
         */
        readerType ReadType;
>>>>>>> Stashed changes

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

<<<<<<< Updated upstream
    /**
     * @brief This property specifies whether to use a ScaleFactor by default.
     *
     */
    bool AutoScale;
=======
        /**
         * @brief StartPoint This property specifies the start point for a PV slice.
         */
        std::pair<int, int> StartPoint;

        /**
         * @brief EndPoint This property specifies the end point for a PV slice.
         */
        std::pair<int, int> EndPoint;

        /**
         * @brief ZBounds This property specifies the z-bounds for a PV slice.
         */
        std::pair<int, int> ZBounds;

        /**
         * @brief This property specifies if the reader must read a subset of the
         * data.
         *
         */
        bool ReadSubExtent;
>>>>>>> Stashed changes

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

<<<<<<< Updated upstream
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
=======
        vtkNew<vtkFloatArray>CalculateMoment(int order);
        int MomentMapRequestInfo(int ProcId, vtkInformationVector* outVec);
        int MomentMapRequestData(int ProcId, vtkInformationVector* outVec);
        int PVSliceRequestInfo(int ProcId, vtkInformationVector* outVec);
        int PVSliceRequestData(int ProcId, vtkInformationVector* outVec);

        double CDELT3, initSlice;

        int Convert2DIndexToLinear(const int row, const int col, const int rowWidth){
            return rowWidth * row + col;
        }
>>>>>>> Stashed changes
};

/**
 * Linear interpolation
 */
template<typename T>
static const T interpolate(const T& fromPointA, const T& toPointB, float factor){
<<<<<<< Updated upstream
    return ((1.0 - factor) * fromPointA) + (factor * toPointB);
=======
        return ((1.0 - factor) * fromPointA) + (factor * toPointB);
>>>>>>> Stashed changes
}

/**
 * Bilinear interpolation
 */
template<typename T>
static const T interpolate(const T& fromPointA1, const T& fromPointA2, const T& toPointB1, const T& toPointB2, float factor1, float factor2){
<<<<<<< Updated upstream
    return interpolate(interpolate(fromPointA1, fromPointA2, factor1), interpolate(toPointB1, toPointB2, factor1), factor2);
}
#endif
=======
        return interpolate(interpolate(fromPointA1, fromPointA2, factor1), interpolate(toPointB1, toPointB2, factor1), factor2);
}

#endif //__vtkFitsReader_h
>>>>>>> Stashed changes
