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

#include "vtkProcessModule.h"

// vtkCxxRevisionMacro(vtkFitsReader, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkFitsReader);

//----------------------------------------------------------------------------
vtkFitsReader::vtkFitsReader()
{
    this->fptr = nullptr;
    this->ReadStatus = 0;
    this->controller = vtkMultiProcessController::GetGlobalController();

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

int vtkFitsReader::RequestInformation(
  vtkInformation*, vtkInformationVector**,
  vtkInformationVector* outVec)
{
    cerr << "RequestInformation " << FileName << " " << endl;
    vtkInformation* outInfo = outVec->GetInformationObject(0);

    if (fits_open_data(&this->fptr, FileName, READONLY, &ReadStatus))
    {
        vtkErrorMacro("vtkFITSReader::ExecuteInformation: ERROR IN CFITSIO! Error reading"
                      " "
                      << FileName << ": \n");
        fits_report_error(stderr, this->ReadStatus);
        return -1;
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
    outInfo->Set(CAN_PRODUCE_SUB_EXTENT(), 1);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), dataExtent, 6);

    //this->vtkMPIImageReader::ExecuteInformation();

    if (fits_close_file(this->fptr, &this->ReadStatus))
    {
        fits_report_error(stderr, this->ReadStatus);
    }

    return 1;
}

//------------------------------------------------------------------------------
// This function reads a data from a file.  The datas extent/axes
// are assumed to be the same as the file extent/order.
int vtkFitsReader::RequestData(  vtkInformation*, vtkInformationVector**,  vtkInformationVector* outVec)
{
    cerr << "RequestData" << endl;
    vtkInformation* outInfo = outVec->GetInformationObject(0);
    vtkImageData* output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    int extent[6];
    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent);
    output->SetExtent(extent);
    std::stringstream msg;

    vtkProcessModule *ParInfo;
    if (ParInfo->GetPartitionId()==0)
      msg<<"Num of Processors = "<<ParInfo->GetNumberOfLocalPartitions()<<endl;
    msg<< "subext from rank "<< ParInfo->GetPartitionId() << " 0: " <<extent[0]<<" 1: "<<extent[1]<<" 2: "<<extent[2]<<" 3: "<<extent[3]<<" 4: "<<extent[4]<<" 5: "<<extent[5] <<endl;

    vtkImageData *data = this->AllocateOutputData(output, outInfo);

    if (data == nullptr)
    {
        vtkErrorMacro(<< "vtkFITSReader::ExecuteDataWithInformation: "
                         "data not allocated.");
        return -1;
    }

    if (this->GetFileName() == nullptr)
    {
        vtkErrorMacro(<< "vtkFITSReader::ExecuteDataWithInformation: "
                         "Either a FileName or FilePrefix must be specified.");
        return -1;
    }

    // Open the fits.  Yes, this means that the file is being opened
    // twice: once by RequestInformation, and once here
    if (fits_open_data(&this->fptr, this->GetFileName(), READONLY, &this->ReadStatus))
    {
        vtkErrorMacro("vtkFITSReader::ExecuteDataWithInformation: "
                      "ERROR IN CFITSIO! Error reading "
                      << this->GetFileName() << ":\n");
        fits_report_error(stderr, this->ReadStatus);
        return -1;
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
        msg<< naxe[axii]<<endl;
    }
        msg<< ""<<endl;

    //long dim= naxe[0]*naxe[1]*(extent[4]-extent[5]);
    int start_position=(naxe[0]*naxe[1]*extent[4])+1;
    msg<< "i am the PE"<< ParInfo->GetPartitionId() << " and i read " <<dim<<" elements from "<< start_position<<endl;

    cerr<<msg.str()<<endl;

    float nullptrval = NAN;
    int anynullptr;
    if (fits_read_img(this->fptr, TFLOAT, start_position, dim, &nullptrval, ptr, &anynullptr, &this->ReadStatus))
    {
        fits_report_error(stderr, this->ReadStatus);
        vtkErrorMacro(<< "vtkFITSReader::ExecuteDataWithInformation: data is nullptr.");
        return -1;
    }

    if (fits_close_file(this->fptr, &this->ReadStatus))
    {
        vtkErrorMacro("vtkFITSReader::ExecuteDataWithInformation: ERROR IN CFITSIO! Error closing "
                      << this->GetFileName() << ":\n");
        fits_report_error(stderr, this->ReadStatus);
                return -1;

    }

    return 1;
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