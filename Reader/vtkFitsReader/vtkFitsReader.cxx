#include "vtkFitsReader.h"
#include "vtkCommand.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkFloatArray.h"
#include <cmath>
#include "vtkPointData.h"

#include <stdlib.h> /* atof */
#include <string>
#include <vector>
#include <sstream>
#include <vtkErrorCode.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

// vtkCxxRevisionMacro(vtkFitsReader, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkFitsReader);

//----------------------------------------------------------------------------
vtkFitsReader::vtkFitsReader()
{
    this->fptr = nullptr;
    this->ReadStatus = 0;
    for (int i = 0; i < 3; i++)
    {
        naxes[i] = 10;
    }
}

//----------------------------------------------------------------------------
vtkFitsReader::~vtkFitsReader()
{
    this->SetFileName(0);
}

//------------------------------------------------------------------------------
void vtkFitsReader::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}

int vtkFitsReader::CanReadFile(const char *fname)
{
    return 1;
}

void vtkFitsReader::ExecuteInformation()
{
    cerr << "ExecuteInformation " << FileName << " " << endl;
    if (fits_open_data(&this->fptr, FileName, READONLY, &ReadStatus))
    {
        vtkErrorMacro("vtkFITSReader::ExecuteInformation: ERROR IN CFITSIO! Error reading"
                      " "
                      << FileName << ": \n");
        fits_report_error(stderr, this->ReadStatus);
        return;
    }

    this->SetPointDataType(vtkDataSetAttributes::SCALARS);
    this->SetNumberOfComponents(1);

    this->SetDataType(VTK_FLOAT);
    this->SetDataScalarType(VTK_FLOAT);

    // Set axis information
    int dataExtent[6] = {0};
    double spacings[3] = {0.};
    double origin[3] = {0.};
    char naxis[10];
    int nfound = 0;
    char comment[120];

    if (fits_read_key(this->fptr, TSTRING, "NAXIS", naxis, comment, &ReadStatus))
        printerror(this->ReadStatus);

    cerr << "naxis " << naxis << endl;
    cerr << "comment " << comment << endl;

    long naxes[boost::lexical_cast<int>(naxis)];
    if (fits_read_keys_lng(this->fptr, "NAXIS", 1, 3, naxes, &nfound, &ReadStatus))
        printerror(this->ReadStatus);

    // calculate the dataExtent and setting the Spacings and Origin
    for (unsigned int axii = 0; axii < boost::lexical_cast<int>(naxis); axii++)
    {
        dataExtent[2 * axii] = 0;
        dataExtent[2 * axii + 1] = naxes[axii] - 1; // StringToInt(this->GetHeaderValue(("SlicerAstro.NAXIS" + IntToString(axii + 1)).c_str())) - 1;
        origin[axii] = 0.0;
        spacings[axii] = 1.0;
    }

    this->SetDataExtent(dataExtent);
    this->SetDataSpacing(spacings);
    this->SetDataOrigin(origin);

    this->vtkImageReader2::ExecuteInformation();

    if (fits_close_file(this->fptr, &this->ReadStatus))
    {
        fits_report_error(stderr, this->ReadStatus);
    }
}

//------------------------------------------------------------------------------
// This function reads a data from a file.  The datas extent/axes
// are assumed to be the same as the file extent/order.
void vtkFitsReader::ExecuteDataWithInformation(vtkDataObject *output, vtkInformation *outInfo)
{
    cerr << "ExecuteDataWithInformation" << endl;

    if (this->GetOutputInformation(0))
    {
        this->GetOutputInformation(0)->Set(
            vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),
            this->GetOutputInformation(0)->Get(
                vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()),
            6);
    }
    vtkImageData *data = this->AllocateOutputData(output, outInfo);

    if (data == nullptr)
    {
        vtkErrorMacro(<< "vtkFITSReader::ExecuteDataWithInformation: "
                         "data not allocated.");
        return;
    }

    if (this->GetFileName() == nullptr)
    {
        vtkErrorMacro(<< "vtkFITSReader::ExecuteDataWithInformation: "
                         "Either a FileName or FilePrefix must be specified.");
        return;
    }

    // Open the fits.  Yes, this means that the file is being opened
    // twice: once by ExecuteInformation, and once here
    if (fits_open_data(&this->fptr, this->GetFileName(), READONLY, &this->ReadStatus))
    {
        vtkErrorMacro("vtkFITSReader::ExecuteDataWithInformation: "
                      "ERROR IN CFITSIO! Error reading "
                      << this->GetFileName() << ":\n");
        fits_report_error(stderr, this->ReadStatus);
        return;
    }

    data->GetPointData()->GetScalars()->SetName("FITSImage");
    // Get data pointer
    void *ptr = nullptr;
    ptr = data->GetPointData()->GetScalars()->GetVoidPointer(0);
    this->ComputeDataIncrements();
    unsigned int naxes = data->GetDataDimension();

    cerr << "naxes" << naxes << endl;
    int naxe[naxes];
    data->GetDimensions(naxe);

    int dim = 1;
    for (unsigned int axii = 0; axii < naxes; axii++)
    {
        dim *= naxe[axii];
    }

    cerr << "dim" << dim << endl;

    float nullptrval = NAN;
    int anynullptr;
    if (fits_read_img(this->fptr, TFLOAT, 1, dim, &nullptrval, ptr, &anynullptr, &this->ReadStatus))
    {
        fits_report_error(stderr, this->ReadStatus);
        vtkErrorMacro(<< "vtkFITSReader::ExecuteDataWithInformation: data is nullptr.");
        return;
    }

    if (fits_close_file(this->fptr, &this->ReadStatus))
    {
        vtkErrorMacro("vtkFITSReader::ExecuteDataWithInformation: ERROR IN CFITSIO! Error closing "
                      << this->GetFileName() << ":\n");
        fits_report_error(stderr, this->ReadStatus);
    }
}

// Note: from cookbook.c in fitsio distribution.
void vtkFitsReader::printerror(int status)
{

    cerr << "vtkFitsReader ERROR.";
    if (status)
    {
        fits_report_error(stderr, status); /* print error report */
        exit(status);                      /* terminate the program, returning error status */
    }
    return;
}