// .NAME vtkFitsReader - read structured points from FITS file.
// .SECTION Description
// vtkFitsReader is a source object that reads FITS data files
// .SECTION Caveats
// Uses CFITSIO v4.1.0 (http://heasarc.gsfc.nasa.gov/docs/software/fitsio)

#ifndef __vtkFitsReader_h
#define __vtkFitsReader_h

#include <vtkFloatArray.h>
#include <vtkMPIImageReader.h>
#include <vtkNew.h>

class vtkTable;

class VTK_EXPORT vtkFitsReader : public vtkMPIImageReader
{

    public:
        enum imageType { EMPTY, FITS2DIMAGE, FITS3DIMAGE };
        enum readerType { RAW, MOMENTMAP };
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
            return "Flexible Image Transport System (FITS)";
        }

        vtkGetMacro(ReadType, int);
        /**
         * @brief SetReadType
         * Function to set the ReadType property.
         * @param type
         * 0 is reading the file raw.
         * 1 is reading a moment map (see MomentOrder).
         */
        void SetReadType(int type) {this->ReadType = (readerType) type;};

        vtkGetMacro(MomentOrder, int);
        vtkSetMacro(MomentOrder, int);

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

        int GetImgType() {return this->ImgType;};

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
         * @brief ReadType This property specifies what algorithm the reader is using to read the FITS file.
         * 0 is reading the file raw.
         * 1 is reading a moment map (see MomentOrder).
         */
        readerType ReadType;

        /**
         * @brief This property specifies if the FITS file is an image (2D) [1] or a cube (3D) [2].
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
         * @brief MomentOrder This property specifies which order of the moment map to read.
         */
        int MomentOrder;

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

        vtkNew<vtkFloatArray>CalculateMoment(int order);
        int MomentMapRequestInfo(int ProcId, vtkInformationVector* outVec);
        int MomentMapRequestData(int ProcId, vtkInformationVector* outVec);

        double CDELT3, initSlice;
};

#endif //__vtkFitsReader_h
