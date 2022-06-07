#include "vtkFitsReader.h"

#include <fitsio.h>

// #include <boost/algorithm/string.hpp>
// #include <boost/lexical_cast.hpp>

#include <vtkCommand.h>
#include <vtkErrorCode.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkProcessModule.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkStringArray.h>
#include "vtkTable.h"

#include <vtkInformationStringKey.h>

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

// vtkCxxRevisionMacro(vtkFitsReader, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkFitsReader);

vtkFitsReader::vtkFitsReader()
{
    this->FileName = NULL;
    this->SetNumberOfOutputPorts(2);

#ifndef NDEBUG
    this->DebugOn();
#endif
}

vtkFitsReader::~vtkFitsReader()
{
    this->SetFileName(0);
}

void vtkFitsReader::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}

int vtkFitsReader::CanReadFile(const char *fname)
{
    return 1;
}

int vtkFitsReader::ReadFITSHeader()
{
    vtkNew<vtkStringArray> hName;
    hName->SetName("Name");
    table->AddColumn(hName);

    vtkNew<vtkStringArray> hValue;
    hValue->SetName("Value");
    table->AddColumn(hValue);

    fitsfile *fptr;
    int status = 0;
    if (fits_open_data(&fptr, FileName, READONLY, &status))
    {
        fits_report_error(stderr, status);
        return 1;
    }

    // Get number of keys in header
    int nKeys = 0;
    if (fits_get_hdrspace(fptr, &nKeys, 0, &status))
    {
        fits_report_error(stderr, status);
        return 2;
    }

    // Get header keys and values
    char name[80];
    char value[80];

    table->SetNumberOfRows(static_cast<vtkIdType>(nKeys));

    for (int i = 1; i <= nKeys; ++i)
    {
        if (fits_read_keyn(fptr, i, name, value, 0, &status))
        {
            fits_report_error(stderr, status);
            return 3;
        }

        std::string sName(name);
        std::string sValue(value);

        table->SetValue(static_cast<vtkIdType>(i - 1), 0, vtkVariant(sName));
        table->SetValue(static_cast<vtkIdType>(i - 1), 1, vtkVariant(sValue));
    }

    fits_close_file(fptr, &status);
    return 0;
}

int vtkFitsReader::RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *outVec)
{
    vtkProcessModule *ProcInfo;
    int ProcId = ProcInfo->GetPartitionId();
    vtkDebugMacro(<< this->GetClassName() << " (" << ProcId << "): RequestInformation " << FileName);

    fitsfile *fptr;
    int ReadStatus = 0;
    if (fits_open_data(&fptr, FileName, READONLY, &ReadStatus))
    {
        vtkErrorMacro(<< this->GetClassName() << " (" << ProcId << ") [CFITSIO] Error fits_open_data");
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
        vtkErrorMacro(<< this->GetClassName() << " (" << ProcId << ") [CFITSIO] Error fits_get_img_param");
        fits_report_error(stderr, ReadStatus);
        return 0;
    }

    if (fits_close_file(fptr, &ReadStatus))
    {
        vtkErrorMacro(<< this->GetClassName() << " (" << ProcId << ") [CFITSIO] Error fits_close_file");
        fits_report_error(stderr, ReadStatus);
        // We should have axes information, so we do not abort (i.e. no return here)
    }

    // Calculate the DataExtent and set the Spacing and Origin
    double spacings[3] = {1.0};
    double origin[3] = {0.0};
    int dataExtent[6];
    for (int axii = 0; axii < naxis; ++axii)
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

    if (ProcId == 0)
    {

        ReadFITSHeader();
        vtkDebugMacro(<< this->GetClassName() << " (" << ProcId << ")\n# of processors: " << ProcInfo->GetNumberOfLocalPartitions()
                      << "\nFileName: " << FileName
                      << "\nImgType: " << imgtype
                      << "\nNAXIS: " << naxis
                      << "\nNAXIS = [" << naxes[0] << ", " << naxes[1] << ", " << naxes[2] << "]"
                      << "\nDataExtent = [" << dataExtent[0] << ", " << dataExtent[1] << ", " << dataExtent[2]
                      << ", " << dataExtent[3] << ", " << dataExtent[4] << ", " << dataExtent[5] << "]");
    }

    return 1;
}

int vtkFitsReader::RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *outVec)
{
    vtkProcessModule *ProcInfo;
    int ProcId = ProcInfo->GetPartitionId();

    if (this->GetFileName() == nullptr)
    {
        vtkErrorMacro(<< this->GetClassName() << " (" << ProcId << "): Either a FileName or FilePrefix must be specified.");
        return 0;
    }

    // Get Data Extent assigned to this process
    int dataExtent[6] = {0, -1, 0, -1, 0, -1};
    vtkInformation *outInfo = outVec->GetInformationObject(0);
    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), dataExtent);
    vtkDebugMacro(<< this->GetClassName() << " (" << ProcId << "): RequestData " << FileName
                  << "\n\tDataExtent = [" << dataExtent[0] << ", " << dataExtent[1] << ", " << dataExtent[2]
                  << ", " << dataExtent[3] << ", " << dataExtent[4] << ", " << dataExtent[5] << "]");

    vtkImageData *data = this->AllocateOutputData(outInfo->Get(vtkDataObject::DATA_OBJECT()), outInfo);
    if (data == nullptr)
    {
        vtkErrorMacro(<< this->GetClassName() << " (" << ProcId << "): Data not allocated.");
        return 0;
    }

    // Create the string to open the sub-region of the fits file
    int size = strlen(FileName) * 2;
    char fn[size];
    snprintf(fn, size, "%s[%d:%d, %d:%d, %d:%d]",
             FileName, dataExtent[0] + 1, dataExtent[1] + 1, dataExtent[2] + 1,
             dataExtent[3] + 1, dataExtent[4] + 1, dataExtent[5] + 1);
    vtkDebugMacro(<< this->GetClassName() << " (" << ProcId << "): RequestData is opening the FITS with the following string: " << fn);

    fitsfile *fptr;
    int ReadStatus = 0;
    if (fits_open_data(&fptr, fn, READONLY, &ReadStatus))
    {
        vtkErrorMacro(<< this->GetClassName() << " (" << ProcId << ") [CFITSIO] Error fits_open_data");
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
        vtkErrorMacro(<< this->GetClassName() << " (" << ProcId << ") [CFITSIO] Error fits_get_img_param");
        fits_report_error(stderr, ReadStatus);
        return 0;
    }

    // Set Data Extend and get Data Pointer
    data->SetExtent(dataExtent);
    data->GetPointData()->GetScalars()->SetName("FITSImage");
    float *ptr = static_cast<float *>(data->GetPointData()->GetScalars()->GetVoidPointer(0));
    this->ComputeDataIncrements();

    long fpixel[3] = {1, 1, 1};
    long long nels = naxes[0] * naxes[1] * naxes[2];
    if (fits_read_pix(fptr, TFLOAT, fpixel, nels, 0, ptr, 0, &ReadStatus))
    {
        vtkErrorMacro(<< this->GetClassName() << " (" << ProcId << ") [CFITSIO] Error fits_read_pix");
        fits_report_error(stderr, ReadStatus);
        return 0;
    }

    if (fits_close_file(fptr, &ReadStatus))
    {
        vtkErrorMacro(<< this->GetClassName() << " (" << ProcId << ") [CFITSIO] Error fits_close_file");
        fits_report_error(stderr, ReadStatus);
        // We should have axes information, so we do not abort (i.e. no return failure here)
    }

    if (ProcId == 0)
    {
        auto output = vtkTable::GetData(outVec, 1);
        output->DeepCopy(table);
    }

    return 1;
}

int vtkFitsReader::FillOutputPortInformation(int port, vtkInformation *info)
{
    switch (port)
    {
    case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
        break;
    case 1:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
        break;
    }

    return 1;
}