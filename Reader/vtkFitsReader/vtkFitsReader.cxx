#include "vtkFitsReader.h"

#include <fitsio.h>

// #include <boost/algorithm/string.hpp>
// #include <boost/lexical_cast.hpp>

#include <vtkCommand.h>
#include <vtkErrorCode.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationStringKey.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkProcessModule.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkStringArray.h>
#include <vtkVariantArray.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

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
    char name[FLEN_KEYWORD];
    char value[FLEN_VALUE];
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
        vtkDebugMacro(<< this->GetClassName() << " (" << ProcId << ")"
                      << "\n  # of processors: " << ProcInfo->GetNumberOfLocalPartitions()
                      << "\n  FileName: " << FileName
                      << "\n  ImgType: " << imgtype
                      << "\n  NAXIS: " << naxis
                      << "\n  NAXIS = [" << naxes[0] << ", " << naxes[1] << ", " << naxes[2] << "]"
                      << "\n  DataExtent = [" << dataExtent[0] << ", " << dataExtent[1] << ", " << dataExtent[2]
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
                  << " - DataExtent = [" << dataExtent[0] << ", " << dataExtent[1] << ", " << dataExtent[2]
                  << ", " << dataExtent[3] << ", " << dataExtent[4] << ", " << dataExtent[5] << "]");

    vtkImageData *data = this->AllocateOutputData(outInfo->Get(vtkDataObject::DATA_OBJECT()), outInfo);
    if (data == nullptr)
    {
        vtkErrorMacro(<< this->GetClassName() << " (" << ProcId << "): Data not allocated.");
        return 0;
    }

    fitsfile *fptr;
    int ReadStatus = 0;
    if (fits_open_data(&fptr, FileName, READONLY, &ReadStatus))
    {
        vtkErrorMacro(<< this->GetClassName() << " (" << ProcId << ") [CFITSIO] Error fits_open_data");
        fits_report_error(stderr, ReadStatus);
        return 0;
    }

    // Set Data Extend and get Data Pointer
    data->SetExtent(dataExtent);
    data->GetPointData()->GetScalars()->SetName("FITSImage");
    float *ptr = static_cast<float *>(data->GetPointData()->GetScalars()->GetVoidPointer(0));
    this->ComputeDataIncrements();

    // Read the extent from the FITS file
    long fP[] = {dataExtent[0] + 1, dataExtent[2] + 1, dataExtent[4] + 1, 1};
    long lP[] = {dataExtent[1] + 1, dataExtent[3] + 1, dataExtent[5] + 1, 1};
    long inc[] = {1, 1, 1, 1};
    float nulval = 1e-30;
    int anynul = 0;
    if (fits_read_subset(fptr, TFLOAT, fP, lP, inc, &nulval, ptr, &anynul, &ReadStatus))
    {
        vtkErrorMacro(<< this->GetClassName() << " (" << ProcId << ") [CFITSIO] Error fits_read_subset");
        fits_report_error(stderr, ReadStatus);
        return 0;
    }

    // Get FITS Statistics
    long dimX = dataExtent[1] - dataExtent[0] + 1;
    long dimY = dataExtent[3] - dataExtent[2] + 1;
    long dimZ = dataExtent[5] - dataExtent[4] + 1;
    double mean = 0;
    double rms = 0;
    long goodpix = 0;
    if (fits_img_stats_float(ptr, dimX * dimY, dimZ, 1, nulval, &goodpix, 0, 0, 0, &rms, 0, 0, 0, 0, &ReadStatus))
    {
        vtkErrorMacro(<< this->GetClassName() << "[CFITSIO] Error fits_img_stats_float");
        fits_report_error(stderr, ReadStatus);
    }
    vtkDebugMacro(<< this->GetClassName() << " (" << ProcId << ") FITS Stats"
                  << "\n  anynul = " << anynul
                  << "\n  RMS: " << rms << "\n  THRESHOLD: " << (3.0 * rms)
                  << "\n  GOOD PIXELS: " << goodpix);

    if (fits_close_file(fptr, &ReadStatus))
    {
        vtkErrorMacro(<< this->GetClassName() << " (" << ProcId << ") [CFITSIO] Error fits_close_file");
        fits_report_error(stderr, ReadStatus);
        // We should have read the data, so we do not abort (i.e. no return failure here)
    }

    vtkTable *output = vtkTable::GetData(outVec, 1);

    if (ProcInfo->GetNumberOfLocalPartitions() == 1)
    {
        // We have the RMS value, so we put it in the output table.
        output->DeepCopy(table);

        vtkNew<vtkVariantArray> rmsRow;
        rmsRow->InsertNextValue(vtkVariant(std::string("RMS")));
        rmsRow->InsertNextValue(vtkVariant(rms));
        output->InsertNextRow(rmsRow);
    }
    else
    {
        // We have partial RMS values, so we put the MeanSquare values in the output tables

        if (ProcId == 0)
        {
            // Proc #0 outputs the entire FITS Header and the number of partial values
            output->DeepCopy(table);
            vtkNew<vtkVariantArray> numberOfValues;
            numberOfValues->InsertNextValue(vtkVariant(std::string("MSn")));
            numberOfValues->InsertNextValue(vtkVariant(ProcInfo->GetNumberOfLocalPartitions()));
            output->InsertNextRow(numberOfValues);
        }
        else
        {
            // Others provide just the partial MeanSquare, but we have to define the number of columns
            vtkNew<vtkStringArray> hName;
            hName->SetName("Name");
            output->AddColumn(hName);

            vtkNew<vtkStringArray> hValue;
            hValue->SetName("Value");
            output->AddColumn(hValue);
        }

        double meanSquare = rms * rms;
        vtkNew<vtkVariantArray> msRow;
        msRow->InsertNextValue(vtkVariant(std::string("MS" + std::to_string(ProcId))));
        msRow->InsertNextValue(vtkVariant(meanSquare));
        output->InsertNextRow(msRow);
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