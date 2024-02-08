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

    double tempCDELT3, tempCRVAL3, tempCPIX3;

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
        if (sName == "CDELT3")
            tempCDELT3 = std::stod(sValue);
        if (sName == "CRPIX3")
            tempCPIX3 = std::stod(sValue);
        if (sName == "CRVAL3")
            tempCRVAL3 = std::stod(sValue);
    }

    this->CDELT3 = tempCDELT3;
    this->initSlice = tempCRVAL3 - (tempCDELT3 * (tempCPIX3 - 1));


    fits_close_file(fptr, &status);
    return 0;
}

int vtkFitsReader::RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *outVec)
{
    vtkProcessModule *ProcInfo;
    int ProcId = ProcInfo->GetPartitionId();
    vtkDebugMacro(<< "(#" << ProcId << ") RequestInformation " << FileName);

    if (this->GetFileName() == nullptr)
    {
        vtkErrorMacro(<< "(#" << ProcId << ") Either a FileName or FilePrefix must be specified.");
        return 0;
    }
    if (ProcId == 0)
    {
        switch (this->GetReadType()) {
        case readerType::RAW:
            vtkDebugMacro(<< " [ReqInfoReadType] Reading file as RAW:");
            break;
        case readerType::MOMENTMAP:
            vtkDebugMacro(<< " [ReqInfoReadType] Reading file as RAW:");
            return MomentMapRequestInfo(ProcId, outVec);
            break;
        default:
            break;
        }
    }
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
        ImgType = imageType::FITS2DIMAGE;
    }
    else if (naxis >= 3)
    {
        ImgType = imageType::FITS3DIMAGE;
    }
    else
    {
        vtkErrorMacro(<< "(#" << ProcId << ") [FITS Header] Error: header proclaims erroneous number of axes.");
        return 0;
    }

    int dataExtent[6];
    // Calculate and adjust DataExtent
    if (ImgType == imageType::FITS2DIMAGE)
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

        if (ImgType == imageType::FITS3DIMAGE)
        {
            long dimZ = dataExtent[5] - dataExtent[4] + 1;
            nels *= dimZ;
        }
        size_t size = sizeof(float) * nels;
        size_t maxSize = static_cast<unsigned long>(CubeMaxSize) * 1024UL * 1024UL;

        if (size > maxSize)
        {
            int factor;
            if (ImgType == imageType::FITS2DIMAGE)
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
                      << "\n  FileName: " << FileName << "\n  Reading as: Raw\n  ImgType: " << imgtype << "\n  NAXIS: " << naxis
                      << "\n  NAXIS = [" << naxes[0] << ", " << naxes[1] << ", " << naxes[2] << "]"
                      << "\n  AutoScale: " << AutoScale << "\n  ScaleFactor: " << ScaleFactor << "\n  DataExtent = ["
                      << dataExtent[0] << ", " << dataExtent[1] << ", " << dataExtent[2] << ", " << dataExtent[3]
                      << ", " << dataExtent[4] << ", " << dataExtent[5] << "]");
    }

    return 1;
}

int vtkFitsReader::MomentMapRequestInfo(int ProcId, vtkInformationVector *outVec)
{
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

    ImgType = imageType::FITS2DIMAGE;

    int dataExtent[6];
    // Calculate and adjust DataExtent
    dataExtent[0] = dataExtent[2] = dataExtent[4] = dataExtent[5] = 0;
    dataExtent[1] = static_cast<int>(naxes[0] - 1);
    dataExtent[3] = static_cast<int>(naxes[1] - 1);

    vtkInformation *outInfo = outVec->GetInformationObject(0);
    outInfo->Set(CAN_PRODUCE_SUB_EXTENT(), 0);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), dataExtent, 6);

    ReadFITSHeader();
    vtkDebugMacro(<< "(#" << ProcId << ") FITS Info"
                  << "\n  FileName: " << FileName << "\n  Reading as: MomentMap order " << this->MomentOrder << "\n  ImgType: " << imgtype << "\n  NAXIS: " << naxis
                  << "\n  NAXIS = [" << naxes[0] << ", " << naxes[1] << ", " << naxes[2] << "]"
                  << "\n  AutoScale: " << AutoScale << "\n  ScaleFactor: " << ScaleFactor << "\n  DataExtent = ["
                  << dataExtent[0] << ", " << dataExtent[1] << ", " << dataExtent[2] << ", " << dataExtent[3]
                  << ", " << dataExtent[4] << ", " << dataExtent[5] << "]");

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
    if (ProcId == 0)
    {
        switch (this->GetReadType()) {
        case readerType::RAW:
            vtkDebugMacro(<< " [ReadType] Reading file as RAW:");
            break;
        case readerType::MOMENTMAP:
            vtkDebugMacro(<< " [ReadType] Reading file as Moment Map:");
            return MomentMapRequestData(ProcId, outVec);
            break;
        default:
            break;
        }
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

int vtkFitsReader::MomentMapRequestData(int ProcId, vtkInformationVector *outVec)
{
    fitsfile *fptr;
    int ReadStatus = 0;
    if (fits_open_data(&fptr, FileName, READONLY, &ReadStatus))
    {
        vtkErrorMacro(<< "(#" << ProcId << ") [CFITSIO] Error fits_open_data");
        fits_report_error(stderr, ReadStatus);
        return 0;
    }

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

    ImgType = imageType::FITS2DIMAGE;

    int dataExtent[6];
    vtkInformation *outInfo = outVec->GetInformationObject(0);

    // Calculate and adjust DataExtent
    dataExtent[0] = dataExtent[2] = dataExtent[4] = dataExtent[5] = 0;
    dataExtent[1] = static_cast<int>(naxes[0] - 1);
    dataExtent[3] = static_cast<int>(naxes[1] - 1);

    outInfo->Set(CAN_PRODUCE_SUB_EXTENT(), 0);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), dataExtent, 6);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), dataExtent, 6);

    auto scalars = this->CalculateMoment(this->GetMomentOrder());

    vtkImageData *data = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    data->SetExtent(dataExtent);
    data->SetOrigin(0.0, 0.0, 0.0);
    data->GetPointData()->SetScalars(scalars);
    std::string name = scalars->GetName();
    vtkDebugMacro(<< "Returning vtkImageData with scalars named " << name << ".");
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

vtkNew<vtkFloatArray> vtkFitsReader::CalculateMoment(int order)
{
    vtkDebugMacro(<< " Calculating moment map " << order);
    ReadFITSHeader();
    fitsfile *fptr;
    int status = 0;
    if (fits_open_file(&fptr, this->GetFileName(), READONLY, &status))
        vtkErrorMacro(<< "[CFITSIO] Error: " << status << "!");

    int nfound = 0;
    long naxes[3];
    if (fits_read_keys_lng(fptr, "NAXIS", 1, 3, naxes, &nfound, &status))
        vtkErrorMacro(<< "[CFITSIO] Error: " << status << "!");

    long buffsize = naxes[0] * naxes[1];
    float *buffer = new float[buffsize];
    float *scalars = new float[buffsize];
    std::fill_n(scalars, buffsize, 0);

    int anynull = 0;
    float fpixel = 1, nullval = 0;

    double vdelt = std::abs(CDELT3);

    switch (order) {
    case 0: {
        // the integrated value of the spectrum
        for (int slice = 0; slice < naxes[2]; ++slice) {
            if (fits_read_img(fptr, TFLOAT, fpixel, buffsize, &nullval, buffer, &anynull, &status))
                vtkErrorMacro(<< "[CFITSIO] Error: " << status << "!");

            for (long i = 0; i < buffsize; ++i) {
                if (std::isfinite(buffer[i])) {
                    scalars[i] += buffer[i] * vdelt;
                }
            }

            fpixel += buffsize;
        }

        break;
    }
    case 1: {
        // the intensity weighted coordinate
        float m0[buffsize];
        for (int slice = 0; slice < naxes[2]; ++slice) {
            if (fits_read_img(fptr, TFLOAT, fpixel, buffsize, &nullval, buffer, &anynull, &status))
                vtkErrorMacro(<< "[CFITSIO] Error: " << status << "!");

            for (long i = 0; i < buffsize; ++i) {
                if (std::isfinite(buffer[i])) {
                    m0[i] += buffer[i] * vdelt;
                }
            }

            fpixel += buffsize;
        }

        fpixel = 1;

        for (int slice = 0; slice < naxes[2]; ++slice) {
            double velocityValue = initSlice + CDELT3 * (slice);
            if (fits_read_img(fptr, TFLOAT, fpixel, buffsize, &nullval, buffer, &anynull, &status))
                vtkErrorMacro(<< "[CFITSIO] Error: " << status << "!");

            for (long i = 0; i < buffsize; ++i) {
                auto val = m0[i];
                if (val != 0 && std::isfinite(buffer[i])) {
                    scalars[i] = scalars[i] + (buffer[i] * velocityValue * vdelt) / val;
                }
            }

            fpixel += buffsize;
        }

        break;
    }
    case 2: {
        // the intensity weighted dispersion of the coordinate

        // the intensity weighted coordinate
        float m0[buffsize];
        for (int slice = 0; slice < naxes[2]; ++slice) {
            if (fits_read_img(fptr, TFLOAT, fpixel, buffsize, &nullval, buffer, &anynull, &status))
                vtkErrorMacro(<< "[CFITSIO] Error: " << status << "!");

            for (long i = 0; i < buffsize; ++i) {
                if (std::isfinite(buffer[i])) {
                    m0[i] += buffer[i] * vdelt;
                }
            }

            fpixel += buffsize;
        }

        fpixel = 1;

        // Moment 1
        float m1[buffsize];
        std::fill_n(m1, buffsize, 0);
        for (int slice = 0; slice < naxes[2]; ++slice) {
            double velocityValue = initSlice + CDELT3 * (slice);

            if (fits_read_img(fptr, TFLOAT, fpixel, buffsize, &nullval, buffer, &anynull, &status))
                vtkErrorMacro(<< "[CFITSIO] Error: " << status << "!");

            for (long i = 0; i < buffsize; ++i) {
                auto val = m0[i];
                if (val != 0 && std::isfinite(buffer[i])) {
                    m1[i] = m1[i] + (buffer[i] * velocityValue * vdelt) / val;
                }
            }

            fpixel += buffsize;
        }

        // Moment 2
        fpixel = 1;
        for (int slice = 0; slice < naxes[2]; ++slice) {
            double velocityValue = initSlice + CDELT3 * (slice);

            if (fits_read_img(fptr, TFLOAT, fpixel, buffsize, &nullval, buffer, &anynull, &status))
                vtkErrorMacro(<< "[CFITSIO] Error: " << status << "!");

            for (long i = 0; i < buffsize; ++i) {
                auto val = m0[i];
                if (val != 0 && std::isfinite(buffer[i])) {
                    scalars[i] = scalars[i]
                                 + (buffer[i] * std::pow(velocityValue - m1[i], 2) * vdelt)
                                       / val;
                }
            }

            fpixel += buffsize;
        }

        break;
    }
    case 6: {
        // root mean square of the spectrum (noise map)
        for (int slice = 0; slice < naxes[2]; ++slice) {
            if (fits_read_img(fptr, TFLOAT, fpixel, buffsize, &nullval, buffer, &anynull, &status))
                vtkErrorMacro(<< "[CFITSIO] Error: " << status << "!");

            for (long i = 0; i < buffsize; ++i) {
                if (std::isfinite(buffer[i])) {
                    scalars[i] += buffer[i] * buffer[i];
                }
            }

            fpixel += buffsize;
        }

        for (int i = 0; i < buffsize; ++i) {
            scalars[i] = std::sqrt(scalars[i] / buffsize);
        }

        break;
    }
    case 8: {
        // maximum value of the spectrum (peak map)
        for (int slice = 0; slice < naxes[2]; ++slice) {
            if (fits_read_img(fptr, TFLOAT, fpixel, buffsize, &nullval, buffer, &anynull, &status))
                vtkErrorMacro(<< "[CFITSIO] Error: " << status << "!");

            for (long i = 0; i < buffsize; ++i) {
                if (std::isfinite(buffer[i])) {
                    scalars[i] = fmax(scalars[i], buffer[i]);
                }
            }

            fpixel += buffsize;
        }

        break;
    }
    case 10: {
        // the minimum value of the spectrum
        for (int slice = 0; slice < naxes[2]; ++slice) {
            if (fits_read_img(fptr, TFLOAT, fpixel, buffsize, &nullval, buffer, &anynull, &status))
                vtkErrorMacro(<< "[CFITSIO] Error: " << status << "!");

            for (long i = 0; i < buffsize; ++i) {
                if (std::isfinite(buffer[i])) {
                    scalars[i] = fmin(scalars[i], buffer[i]);
                }
            }

            fpixel += buffsize;
        }

        break;
    }
    }

    delete[] buffer;

    float datamin, datamax;
    double mean, rms;
    if (fits_img_stats_float(scalars, naxes[0], naxes[1], 0, 0, 0, &datamin, &datamax,
                             &mean, &rms, 0, 0, 0, 0, &status)) {
        fits_report_error(stderr, status);
    }

    if (fits_close_file(fptr, &status)) {
        vtkErrorMacro(<< "[CFITSIO] Error: " << status << "!");
    }

    vtkNew<vtkFloatArray> values;
    std::string name = "FITSImage" + std::to_string(order);
    values->SetName(name.c_str());
    values->SetNumberOfComponents(1);
    values->SetVoidArray(scalars, buffsize, 0, vtkAbstractArray::VTK_DATA_ARRAY_DELETE);
    return values;
}
