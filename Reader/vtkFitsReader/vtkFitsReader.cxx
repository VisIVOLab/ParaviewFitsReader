#include "vtkFitsReader.h"

#include "fitsio.h"

#include "vtkCommand.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkImageData.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkErrorCode.h"
#include "vtkProcessModule.h"

#include "vtkSMPTools.h"



#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>

// vtkCxxRevisionMacro(vtkFitsReader, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkFitsReader);

//----------------------------------------------------------------------------
vtkFitsReader::vtkFitsReader()
{
    this->FileName = NULL;
    vtkSMPTools::Initialize();
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

//------------------------------------------------------------------------------
int vtkFitsReader::CanReadFile(const char *fname)
{
    return 1;
}

//------------------------------------------------------------------------------
int vtkFitsReader::RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *outVec)
{
    vtkProcessModule *ProcInfo;
    cout << "RequestInformation " << FileName << " (#" << ProcInfo->GetPartitionId() << ") " << endl;

    fitsfile *fptr;
    int ReadStatus = 0;
    if (fits_open_data(&fptr, FileName, READONLY, &ReadStatus))
    {
        vtkErrorMacro("vtkFitsReader::RequestInformation (# " << ProcInfo->GetPartitionId() << ")\n"
                                                              << "ERROR IN CFITSIO! Error reading "
                                                              << FileName << ":\n");
        fits_report_error(stderr, ReadStatus);
        return 0;
    }

    // Get axis information
    int maxaxis = 3;
    int imgtype = 0;
    int naxis = 0;
    long naxes[maxaxis];

    if (fits_get_img_param(fptr, maxaxis, &imgtype, &naxis, naxes, &ReadStatus))
    {
        vtkErrorMacro("vtkFitsReader::RequestInformation (# " << ProcInfo->GetPartitionId() << ")\n"
                                                              << "ERROR IN CFITSIO! Error reading image param "
                                                              << FileName << ":\n");
        fits_report_error(stderr, ReadStatus);
        return 0;
    }

    if (fits_close_file(fptr, &ReadStatus))
    {
        vtkErrorMacro("vtkFitsReader::RequestInformation (# " << ProcInfo->GetPartitionId() << ")\n"
                                                              << "ERROR IN CFITSIO! Error closing "
                                                              << FileName << ":\n");
        fits_report_error(stderr, ReadStatus);
        // We should have axes information, so we do not abort (i.e. no return here)
    }

    // Calculate the DataExtent and set the Spacing and Origin
    double spacings[3] = {1.0};
    double origin[3] = {0.0};
    int dataExtent[6];
    for (unsigned int axii = 0; axii < naxis; axii++)
    {
        dataExtent[2 * axii] = 0;
        dataExtent[2 * axii + 1] = naxes[axii] - 1;
    }

    this->SetPointDataType(vtkDataSetAttributes::SCALARS);
    this->SetNumberOfComponents(1);
    this->SetDataType(VTK_FLOAT);
    this->SetDataScalarType(VTK_FLOAT);
    this->SetDataExtent(dataExtent);
    this->SetDataSpacing(spacings);
    this->SetDataOrigin(origin);

    vtkInformation *outInfo = outVec->GetInformationObject(0);
    outInfo->Set(CAN_PRODUCE_SUB_EXTENT(), 1);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), dataExtent, 6);



    if (ProcInfo->GetPartitionId() == 0)
    {
        cout << "# of processors: " << ProcInfo->GetNumberOfLocalPartitions() << endl;
        cout << "NumberOfThreads: "<< vtkSMPTools::GetEstimatedNumberOfThreads() << endl;
        cout << "FileName: " << FileName
             << "\nImgType: " << imgtype
             << "\nNAXIS: " << naxis
             << "\nNAXIS = [" << naxes[0] << ", " << naxes[1]
             << ", " << naxes[2] << "]"
             << "\nDataExtent = [" << dataExtent[0] << ", " << dataExtent[1] << ", " << dataExtent[2]
             << ", " << dataExtent[3] << ", " << dataExtent[4] << ", " << dataExtent[5] << "]" << endl;
    }

    return 1;
}

//------------------------------------------------------------------------------
int vtkFitsReader::RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *outVec)
{
    int dataExtent[6] = {0, -1, 0, -1, 0, -1};
    vtkInformation *outInfo = outVec->GetInformationObject(0);
    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), dataExtent);

    vtkImageData *outData = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    outData->SetExtent(dataExtent);

    vtkProcessModule *ProcInfo;
    if (ProcInfo->GetPartitionId() == 0)
    {
        cout << "RequestData # of processors: " << ProcInfo->GetNumberOfLocalPartitions() << "\n";
    }
    cout << "RequestData " << FileName << " (#" << ProcInfo->GetPartitionId() << ")\n"
         << "\tDataExtent = [" << dataExtent[0] << ", " << dataExtent[1] << ", " << dataExtent[2]
         << ", " << dataExtent[3] << ", " << dataExtent[4] << ", " << dataExtent[5] << "]" << endl;

    if (this->GetFileName() == nullptr)
    {
        vtkErrorMacro("vtkFitsReader::RequestData (# " << ProcInfo->GetPartitionId() << ") "
                                                       << "Either a FileName or FilePrefix must be specified.");
        return 0;
    }

    vtkImageData *data = this->AllocateOutputData(outData, outInfo);
    if (data == nullptr)
    {
        vtkErrorMacro("vtkFitsReader::RequestData (# " << ProcInfo->GetPartitionId() << ") "
                                                       << "ERROR: data not allocated");
        return 0;
    }

    int size = strlen(FileName) + 100;
    char fn[size];

    snprintf(fn, size, "%s[%d:%d, %d:%d, %d:%d]",
             FileName, dataExtent[0] + 1, dataExtent[1] + 1, dataExtent[2] + 1,
             dataExtent[3] + 1, dataExtent[4] + 1, dataExtent[5] + 1);

    cout << ProcInfo->GetPartitionId() << " is opening the FITS with the following string: " << fn << endl;

    fitsfile *fptr;
    int ReadStatus = 0;
    if (fits_open_data(&fptr, fn, READONLY, &ReadStatus))
    {
        vtkErrorMacro("vtkFitsReader::RequestData (# " << ProcInfo->GetPartitionId() << ")\n"
                                                       << "ERROR IN CFITSIO! Error reading "
                                                       << FileName << ":\n");
        fits_report_error(stderr, ReadStatus);
        return 0;
    }

    // Get axis information
    int maxaxis = 3;
    int imgtype = 0;
    int naxis = 0;
    long naxes[maxaxis];

    if (fits_get_img_param(fptr, maxaxis, &imgtype, &naxis, naxes, &ReadStatus))
    {
        vtkErrorMacro("vtkFitsReader::RequestData (# " << ProcInfo->GetPartitionId() << ")\n"
                                                       << "ERROR IN CFITSIO! Error reading image param "
                                                       << FileName << ":\n");
        fits_report_error(stderr, ReadStatus);
        return 0;
    }

    // Get Data Pointer
    data->GetPointData()->GetScalars()->SetName("FITSImage");
    void *ptr = nullptr;
    ptr = data->GetPointData()->GetScalars()->GetVoidPointer(0);
    this->ComputeDataIncrements();

    long fpixel[3] = {1, 1, 1};
    long long nels = naxes[0] * naxes[1] * naxes[2];
    if (fits_read_pix(fptr, TFLOAT, fpixel, nels, NULL, ptr, NULL, &ReadStatus))
    {
        vtkErrorMacro("vtkFitsReader::RequestData (# " << ProcInfo->GetPartitionId() << ")\n"
                                                       << "ERROR IN CFITSIO! Error reading pixels "
                                                       << FileName << ":\n");
        fits_report_error(stderr, ReadStatus);
        return 0;
    }

    /*
        if (fits_read_img(this->fptr, TFLOAT, start_position, dim, &nullptrval, ptr, &anynullptr, &this->ReadStatus))
        {
            fits_report_error(stderr, this->ReadStatus);
            vtkErrorMacro(<< "vtkFITSReader::ExecuteDataWithInformation: data is nullptr.");
            return -1;
        }
    */

    if (fits_close_file(fptr, &ReadStatus))
    {
        vtkErrorMacro("vtkFitsReader::RequestData (# " << ProcInfo->GetPartitionId() << ")\n"
                                                       << "ERROR IN CFITSIO! Error closing "
                                                       << FileName << ":\n");
        fits_report_error(stderr, ReadStatus);
        // We should have axes information, so we do not abort (i.e. no return here)
    }

    return 1;
}
