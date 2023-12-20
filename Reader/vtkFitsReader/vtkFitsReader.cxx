#include "vtkFitsReader.h"

#include <fitsio.h>

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
#include <vtkTable.h>
#include <vtkVariantArray.h>

#include <algorithm>
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
    this->ImgType = imageType::EMPTY;

#ifndef NDEBUG
    this->DebugOn();
#endif
}

int vtkFitsReader::ReadPVSliceData(int ProcID, vtkInformationVector *outVec)
{
    // Get Data Extent assigned to this process
    int dataExtent[6];
    vtkInformation *outInfo = outVec->GetInformationObject(0);
    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), dataExtent);

    int xDim = std::ceil(std::sqrt(std::pow(this->PVEnd.first - this->PVStart.first, 2) + std::pow(this->PVEnd.second - this->PVStart.second, 2)));
    int yDim = this->PVZSubset.second - this->PVZSubset.first + 1;
    dataExtent[0] = dataExtent[2] = dataExtent[4] = dataExtent[5] = 0;
    dataExtent[1] = xDim;
    dataExtent[3] = yDim;
    long nels = xDim * yDim * 1;

    vtkDebugMacro(<< "(#" << ProcID << ") Creating PV Slice from " << FileName << " with DataExtent = [" << 0 << ", "
                  << xDim << ", " << 0 << ", " << yDim << ", " << 0 << ", " << 1 << "]");

    float *ptr;
    try
    {
        ptr = new float[nels];
    }
    catch (const std::bad_alloc &e)
    {
        vtkErrorMacro(<< "(#" << ProcID << ") Data not allocated.");
        return 0;
    }

    fitsfile *fptr;
    int ReadStatus = 0;
    if (fits_open_data(&fptr, FileName, READONLY, &ReadStatus))
    {
        vtkErrorMacro(<< "(#" << ProcID << ") [CFITSIO] Error fits_open_data");
        fits_report_error(stderr, ReadStatus);
        return 0;
    }

    auto botLeft = std::make_pair(std::min(PVStart.first, this->PVEnd.first), std::min(PVStart.second, this->PVEnd.second));
    auto topRight = std::make_pair(std::max(PVStart.first, this->PVEnd.first), std::max(PVStart.second, this->PVEnd.second));
    int columnNumber = topRight.first - botLeft.first;

    double xPixW = ((PVEnd.first - PVStart.first) * 1.0f) / xDim;
    double yPixW = ((PVEnd.second - PVStart.second) * 1.0f) / xDim;
    int sliceNels = (topRight.second - botLeft.second) * columnNumber;

    for (int rowCoord = 0; rowCoord < yDim; ++rowCoord){
        vtkDebugMacro(<< "Entering loop with index " << rowCoord);
        float *slicePtr;
        try
        {
            vtkDebugMacro(<< "Allocating memory for slicePtr at index " << rowCoord);
            slicePtr = new float[sliceNels];
        }
        catch (const std::bad_alloc &e)
        {
            vtkErrorMacro(<< "(#" << ProcID << ") Data not allocated.");
            return 0;
        }

        // Read the extent from the FITS file
        int sliceindex = rowCoord + this->PVZSubset.first + 1;
        long firstPixel[] = {botLeft.first + 1, botLeft.second + 1, sliceindex, 1};
        long lastPixel[] = {topRight.first + 1, topRight.second + 1, sliceindex , 1};
        long increment[] = {1, 1, 1, 1};
        float nulval = 1e-30;
        int anynul = 0;
        if (fits_read_subset(fptr, TFLOAT, firstPixel, lastPixel, increment, &nulval, slicePtr, &anynul, &ReadStatus))
        {
            vtkErrorMacro(<< "(#" << ProcID << ") [CFITSIO] Error fits_read_subset\n" <<
                          "First pixel: [" << firstPixel[0] << ", " << firstPixel[1] << ", " << firstPixel[2] << "]\n" <<
                          "Last pixel: [" << lastPixel[0] << ", " << lastPixel[1] << ", " << lastPixel[2] << "]\n" <<
                          "sliceNels: " << sliceNels);
            fits_report_error(stderr, ReadStatus);
            return 0;
        }
        vtkDebugMacro(<< "Successfully read subset for slice at index " << sliceindex);

        for (int columnCoord = 0; columnCoord < xDim; columnCoord++){
            auto idx = std::make_pair(PVStart.first + (xPixW * rowCoord), PVStart.second + (yPixW * rowCoord));
            double d00, d01, d10, d11;
            d00 = slicePtr[Convert2DIndexToLinear(std::floor(idx.first), std::floor(idx.second), columnNumber)];
            d01 = slicePtr[Convert2DIndexToLinear(std::floor(idx.first), std::floor(idx.second) + 1, columnNumber)];
            d10 = slicePtr[Convert2DIndexToLinear(std::floor(idx.first) + 1, std::floor(idx.second), columnNumber)];
            d11 = slicePtr[Convert2DIndexToLinear(std::floor(idx.first) + 1, std::floor(idx.second) + 1, columnNumber)];
            auto x = idx.first - std::floor(idx.first);
            auto y = idx.second - std::floor(idx.second);
            float pixIJ = interpolate(d00, d01, d10, d11, x, y);

            //We're creating 2D image, so z-axis coord is always 0. Writing to 0th component.
            ptr[Convert2DIndexToLinear(rowCoord, columnCoord, yDim)] = pixIJ;
        }

        vtkDebugMacro(<< "Loop at index " << rowCoord << " completed");
//        delete [] slicePtr;
    }


    if (fits_close_file(fptr, &ReadStatus))
    {
        vtkErrorMacro(<< "(#" << ProcID << ") [CFITSIO] Error fits_close_file");
        fits_report_error(stderr, ReadStatus);
        // We should have read the data, so we do not abort (i.e. no return failure here)
    }

    outInfo->Set(CAN_PRODUCE_SUB_EXTENT(), 1);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), dataExtent, 6);

    vtkNew<vtkFloatArray> scalars;
    scalars->SetName("FITSImage");
    scalars->SetNumberOfComponents(1);
    scalars->SetVoidArray(ptr, nels, 0, vtkAbstractArray::VTK_DATA_ARRAY_DELETE);

    vtkImageData *data = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    data->SetExtent(dataExtent);
    data->SetOrigin(0.0, 0.0, 0.0);
    data->SetSpacing(ScaleFactor, ScaleFactor, ScaleFactor);
    data->GetPointData()->SetScalars(scalars);


//    vtkTable *output = vtkTable::GetData(outVec, 1);

//    if (ProcInfo->GetNumberOfLocalPartitions() == 1)
//    {
//        // We have the RMS value, so we put it in the output table.
//        output->DeepCopy(table);

//        vtkNew<vtkVariantArray> rmsRow;
//        rmsRow->InsertNextValue(vtkVariant(std::string("RMS")));
//        rmsRow->InsertNextValue(vtkVariant(rms));
//        output->InsertNextRow(rmsRow);
//    }
//    else
//    {
//        // We have partial RMS values, so we put the MeanSquare values in the output tables

//        if (ProcID == 0)
//        {
//            // Proc #0 outputs the entire FITS Header and the number of partial values
//            output->DeepCopy(table);
//            vtkNew<vtkVariantArray> numberOfValues;
//            numberOfValues->InsertNextValue(vtkVariant(std::string("MSn")));
//            numberOfValues->InsertNextValue(vtkVariant(ProcInfo->GetNumberOfLocalPartitions()));
//            output->InsertNextRow(numberOfValues);
//        }
//        else
//        {
//            // Others provide just the partial MeanSquare, but we have to define the number of columns
//            vtkNew<vtkStringArray> hName;
//            hName->SetName("Name");
//            output->AddColumn(hName);

//            vtkNew<vtkStringArray> hValue;
//            hValue->SetName("Value");
//            output->AddColumn(hValue);
//        }

//        double meanSquare = rms * rms;
//        vtkNew<vtkVariantArray> msRow;
//        msRow->InsertNextValue(vtkVariant(std::string("MS" + std::to_string(ProcID))));
//        msRow->InsertNextValue(vtkVariant(meanSquare));
//        output->InsertNextRow(msRow);
//    }

    return 1;
}

vtkFitsReader::~vtkFitsReader()
{
    this->SetFileName(0);
}

void vtkFitsReader::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
    os << indent << "ImageType: " << (imageType)ImgType;
    os << indent << "ReadSubExtent: " << std::boolalpha << ReadSubExtent;
    if ((imageType)ImgType == imageType::FITS2DIMAGE)
        os << indent << indent << "SubExtent: " << SubExtent[0] << " " << SubExtent[1] << " " << SubExtent[2] << " "
           << SubExtent[3];
    else if ((imageType)ImgType == imageType::FITS3DIMAGE)
        os << indent << indent << "SubExtent: " << SubExtent[0] << " " << SubExtent[1] << " " << SubExtent[2] << " "
           << SubExtent[3] << " " << SubExtent[4] << " " << SubExtent[5];
    os << indent << "AutoScale: " << std::boolalpha << AutoScale;
    os << indent << "ScaleFactor: " << ScaleFactor;
}

int vtkFitsReader::CanReadFile(const char *fname)
{
    return 1;
}

int vtkFitsReader::ReadFITSHeader()
{
    table->Initialize();

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
    vtkDebugMacro(<< "(#" << ProcId << ") RequestInformation " << FileName);

    fitsfile *fptr;
    int ReadStatus = 0;
    if (fits_open_data(&fptr, FileName, READONLY, &ReadStatus))
    {
        vtkErrorMacro(<< "(#" << ProcId << ") [CFITSIO] Error fits_open_data");
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
        vtkErrorMacro(<< "(#" << ProcId << ") [CFITSIO] Error fits_get_img_param");
        fits_report_error(stderr, ReadStatus);
        return 0;
    }

    if (fits_close_file(fptr, &ReadStatus))
    {
        vtkErrorMacro(<< "(#" << ProcId << ") [CFITSIO] Error fits_close_file");
        fits_report_error(stderr, ReadStatus);
        // We should have axes information, so we do not abort (i.e. no return here)
    }

    // Check if the image is 2D or 3D
    if (naxis == 2)
    {
        this->ImgType = imageType::FITS2DIMAGE;
        naxes[2] = 1;
    }
    else if (naxis >= 3)
    {
        this->ImgType = imageType::FITS3DIMAGE;
    }
    else
    {
        vtkErrorMacro(<< "(#" << ProcId << ") [FITS Header] Error: header proclaims erroneous number of axes.");
        return 0;
    }

    int dataExtent[6];
    // Calculate and adjust DataExtent
    if (this->ImgType == imageType::FITS2DIMAGE)
    {
        dataExtent[0] = dataExtent[2] = dataExtent[4] = dataExtent[5] = 0;
        dataExtent[1] = static_cast<int>(naxes[0] - 1);
        dataExtent[3] = static_cast<int>(naxes[1] - 1);
    }
    else // if (ImgType == imageType::FITS3DIMAGE) always true
    {
        dataExtent[0] = dataExtent[2] = dataExtent[4] = 0;
        dataExtent[1] = static_cast<int>(naxes[0] - 1);
        dataExtent[3] = static_cast<int>(naxes[1] - 1);
        dataExtent[5] = static_cast<int>(naxes[2] - 1);
    }

    if (ReadSubExtent)
    {
        for (int i = 0; i < 6; ++i)
        {
            if (SubExtent[i] != -1)
            {
                dataExtent[i] = SubExtent[i];
            }
        }

        if (ProcId == 0)
        {
            vtkDebugMacro(<< "(#" << ProcId << ")\nReadSubExtent enabled [" << SubExtent[0] << ", " << SubExtent[1]
                          << ", " << SubExtent[2] << ", " << SubExtent[3] << ", " << SubExtent[4] << ", "
                          << SubExtent[5] << "]"
                          << "\nActual SubExtent [" << dataExtent[0] << ", " << dataExtent[1] << ", " << dataExtent[2]
                          << ", " << dataExtent[3] << ", " << dataExtent[4] << ", " << dataExtent[5] << "]");
        }
    }

    // If AutoScale is enabled, calculate ScaleFactor
    if (AutoScale)
    {
        long dimX = dataExtent[1] - dataExtent[0] + 1;
        long dimY = dataExtent[3] - dataExtent[2] + 1;
        long nels = dimX * dimY;

        if (this->ImgType == imageType::FITS3DIMAGE)
        {
            long dimZ = dataExtent[5] - dataExtent[4] + 1;
            nels *= dimZ;
        }
        size_t size = sizeof(float) * nels;
        size_t maxSize = static_cast<unsigned long>(CubeMaxSize) * 1024UL * 1024UL;

        if (size > maxSize)
        {
            int factor;
            if (this->ImgType == imageType::FITS2DIMAGE)
                factor = ceil(sqrt(1.0 * size / maxSize));
            else // if (ImgType == imageType::FITS3DIMAGE) always true
                factor = ceil(cbrt(1.0 * size / maxSize));
            SetScaleFactor(factor);
        }
    }

    if (ScaleFactor > 1)
    {
        for (int i = 0; i < 3; ++i)
        {
            dataExtent[2 * i] /= ScaleFactor;
            dataExtent[2 * i + 1] /= ScaleFactor;
        }
    }

    vtkInformation *outInfo = outVec->GetInformationObject(0);
    outInfo->Set(CAN_PRODUCE_SUB_EXTENT(), 1);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), dataExtent, 6);

    if (ProcId == 0)
    {
        ReadFITSHeader();
        vtkDebugMacro(<< "(#" << ProcId << ") FITS Info"
                      << "\n  # of processors: " << ProcInfo->GetNumberOfLocalPartitions()
                      << "\n  FileName: " << FileName << "\n  ImgType: " << imgtype << "\n  NAXIS: " << naxis
                      << "\n  NAXIS = [" << naxes[0] << ", " << naxes[1] << ", " << naxes[2] << "]"
                      << "\n  AutoScale: " << AutoScale << "\n  ScaleFactor: " << ScaleFactor << "\n  DataExtent = ["
                      << dataExtent[0] << ", " << dataExtent[1] << ", " << dataExtent[2] << ", " << dataExtent[3]
                      << ", " << dataExtent[4] << ", " << dataExtent[5] << "]");
    }

    return 1;
}

int vtkFitsReader::RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *outVec)
{
    vtkProcessModule *ProcInfo;
    int ProcId = ProcInfo->GetPartitionId();

    if (this->GetFileName() == nullptr)
    {
        vtkErrorMacro(<< "(#" << ProcId << ") Either a FileName or FilePrefix must be specified.");
        return 0;
    }

    if (this->GetReadAsPVSlice())
    {
        this->ImgType = imageType::FITS2DIMAGE;
        return this->ReadPVSliceData(ProcId, outVec);
    }

    // Get Data Extent assigned to this process
    int dataExtent[6];
    vtkInformation *outInfo = outVec->GetInformationObject(0);
    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), dataExtent);

    long dimX = dataExtent[1] - dataExtent[0] + 1;
    long dimY = dataExtent[3] - dataExtent[2] + 1;
    long dimZ = dataExtent[5] - dataExtent[4] + 1;
    long nels = dimX * dimY * dimZ;

    vtkDebugMacro(<< "(#" << ProcId << ") RequestData " << FileName << " - DataExtent = [" << dataExtent[0] << ", "
                  << dataExtent[1] << ", " << dataExtent[2] << ", " << dataExtent[3] << ", " << dataExtent[4] << ", "
                  << dataExtent[5] << "]");

    float *ptr;
    try
    {
        ptr = new float[nels];
    }
    catch (const std::bad_alloc &e)
    {
        vtkErrorMacro(<< "(#" << ProcId << ") Data not allocated.");
        return 0;
    }

    fitsfile *fptr;
    int ReadStatus = 0;
    if (fits_open_data(&fptr, FileName, READONLY, &ReadStatus))
    {
        vtkErrorMacro(<< "(#" << ProcId << ") [CFITSIO] Error fits_open_data");
        fits_report_error(stderr, ReadStatus);
        return 0;
    }

    // Read the extent from the FITS file
    long fP[] = {ScaleFactor * dataExtent[0] + 1, ScaleFactor * dataExtent[2] + 1, ScaleFactor * dataExtent[4] + 1, 1};
    long lP[] = {ScaleFactor * dataExtent[1] + 1, ScaleFactor * dataExtent[3] + 1, ScaleFactor * dataExtent[5] + 1, 1};
    long inc[] = {ScaleFactor, ScaleFactor, ScaleFactor, 1};
    float nulval = 1e-30;
    int anynul = 0;
    if (fits_read_subset(fptr, TFLOAT, fP, lP, inc, &nulval, ptr, &anynul, &ReadStatus))
    {
        vtkErrorMacro(<< "(#" << ProcId << ") [CFITSIO] Error fits_read_subset");
        fits_report_error(stderr, ReadStatus);
        return 0;
    }

    // Get FITS Statistics
    double mean = 0;
    double rms = 0;
    long goodpix = 0;
    if (fits_img_stats_float(ptr, dimX * dimY, dimZ, 1, nulval, &goodpix, 0, 0, 0, &rms, 0, 0, 0, 0, &ReadStatus))
    {
        vtkErrorMacro(<< "[CFITSIO] Error fits_img_stats_float");
        fits_report_error(stderr, ReadStatus);
    }

    vtkDebugMacro(<< "(#" << ProcId << ") FITS Stats"
                  << "\n  anynul = " << anynul << "\n  RMS: " << rms << "\n  THRESHOLD: " << (3.0 * rms)
                  << "\n  GOOD PIXELS: " << goodpix);

    if (fits_close_file(fptr, &ReadStatus))
    {
        vtkErrorMacro(<< "(#" << ProcId << ") [CFITSIO] Error fits_close_file");
        fits_report_error(stderr, ReadStatus);
        // We should have read the data, so we do not abort (i.e. no return failure here)
    }

    vtkNew<vtkFloatArray> scalars;
    scalars->SetName("FITSImage");
    scalars->SetNumberOfComponents(1);
    scalars->SetVoidArray(ptr, nels, 0, vtkAbstractArray::VTK_DATA_ARRAY_DELETE);

    vtkImageData *data = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    data->SetExtent(dataExtent);
    data->SetOrigin(0.0, 0.0, 0.0);
    data->SetSpacing(ScaleFactor, ScaleFactor, ScaleFactor);
    data->GetPointData()->SetScalars(scalars);

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

void vtkFitsReader::SetStartPoint(int newStartPointX, int newStartPointY)
{
    PVStart = std::make_pair(newStartPointX, newStartPointY);
}

void vtkFitsReader::SetEndPoint(int newEndPointX, int newEndPointY)
{
    PVEnd = std::make_pair(newEndPointX, newEndPointY);
}

void vtkFitsReader::SetZBounds(int newZMin, int newZMax)
{
    PVZSubset = std::make_pair(newZMin, newZMax);
}
