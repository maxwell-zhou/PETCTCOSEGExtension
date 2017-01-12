#ifndef IMGIO_H
#define IMGIO_H

#include "stdlib.h"
#include "stdio.h"
#include <string>
#include "ImageType.h" 

#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkExtractImageFilter.h>

/*! \class IMGIO
*  \brief The image IO class
*  Include several static member function used to read and write the 3D/4D image
*/

using namespace std;
class IMGIO
{
public:
	inline IMGIO(){};
	inline ~IMGIO(){};

	template <typename Timage, typename TImagePointer> 
		static TImagePointer LoadImg(string filename);
	template <typename Timage, typename TImagePointer>
		static bool  WriteImg(TImagePointer Img, std::string filename); 
	template <typename Timage, typename TImagePointer, typename SizeType, typename IndexType, typename PointType, typename SpacingType> 
		static void ImgInitialization(TImagePointer& InputImage, SizeType InputImageSize, 
		IndexType InputImageIndex, PointType InputImageStart, SpacingType InputImageSpace); 
	template <typename Timage4D, typename Timage3D, typename Timage4DPointer, typename Timage3DPointer>
		static Timage3DPointer Extract3DFrom4DImg(Timage4DPointer Img, int phase);
	template <typename Timage, typename TimagePointer, typename TimageSize, typename TimageIndex>
		static TimagePointer ExtractROI(TimagePointer Img, TimageSize ROI_Size, TimageIndex ROI_Start);
	

	/*! \fn template <typename Timage, typename TImagePointer> static TImagePointer LoadImg(string filename)
	*  \brief An image Loader.
	*  \param filename the name of the image needed to be loaded.
	*  \return An ITK image.
	*/

	/*! \fn static template <typename Timage, typename TImagePointer> static bool  WriteImg(TImagePointer Img, std::string filename)
	*  \brief An image Writer.
	*  \param Img the image needed to be writen.
	*  \param filename the name of the image needed to be write.
	*/

	/*! \fn template <typename Timage, typename TImagePointer, typename SizeType, typename IndexType, typename PointType, typename SpacingType> static void ImgInitialization(TImagePointer& InputImage, SizeType InputImageSize, IndexType InputImageIndex, PointType InputImageStart, SpacingType InputImageSpace)
	*  \brief Initialize an input image.
	*/

	/*! \fn template <typename Timage4D, typename Timage3D, typename Timage4DPointer, typename Timage3DPointer> static Timage3DPointer Extract3DFrom4DImg(Timage4DPointer Img, int phase)
	*  \brief Extract a single phase 3D image from a 4D image
	*  \param Img The input 4D image.
	*  \param phase The phase number of the target 3D image.
	*  \return The 3D phase image
	*/

	/*! \fn 	template <typename Timage, typename TimagePointer, typename TimageSize, typename TimageIndex> static TimagePointer ExtractROI(TimagePointer Img, TimageSize ROI_Size, TimageIndex ROI_Start)
	*  \brief Extract a ROI from image
	*  \param[in] Img The input image.
	*  \param[in] ROI_Size The size of the ROI image.
	*  \param[in] ROI_Start The start index of the ROI image.
	*  \return The ROI image
	*/

	
};


template <typename Timage, typename TImagePointer> 
TImagePointer IMGIO::LoadImg(string filename)
{
	
	typedef itk::ImageFileReader<Timage> ReaderType;
	typename ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename);
	reader->Update();
	return reader->GetOutput();
}

template <typename Timage, typename TImagePointer>
bool IMGIO::WriteImg(TImagePointer Img, std::string filename)
{
	typename itk::ImageFileWriter< Timage >::Pointer Imgwriter = 
		itk::ImageFileWriter< Timage >::New();
	Imgwriter->SetFileName( filename);
	Imgwriter->SetInput(Img);
	Imgwriter->Update();
	return true;
}

template <typename Timage, typename TImagePointer, typename SizeType, typename IndexType, typename PointType, typename SpacingType> 
void IMGIO::ImgInitialization(TImagePointer& InputImage, SizeType InputImageSize, 
								  IndexType InputImageIndex, PointType InputImageStart, SpacingType InputImageSpace)
{
	//Calculate the cost image for resampled Image
	InputImage = Timage::New();
	typename Timage::RegionType Region;
	Region.SetSize(InputImageSize);
	Region.SetIndex(InputImageIndex);
	InputImage->SetLargestPossibleRegion( Region );
	InputImage->SetBufferedRegion( Region );
	InputImage->SetRequestedRegion( Region );
	InputImage->Allocate();
	InputImage->SetOrigin(InputImageStart);
	InputImage->SetSpacing(InputImageSpace);
	typename itk::ImageRegionIterator< Timage>  InputImageIt(InputImage, InputImage->GetRequestedRegion());
	/*
	for(InputImageIt.GoToBegin();!InputImageIt.IsAtEnd();++InputImageIt)
	{
		InputImageIt.Set(0);
	}
	*/
	return;
}

template <typename Timage4D, typename Timage3D, typename Timage4DPointer, typename Timage3DPointer>
Timage3DPointer IMGIO::Extract3DFrom4DImg(Timage4DPointer Img, int phase)
{
	typedef itk::ExtractImageFilter< Timage4D, Timage3D> ExtractFilterType;
	typename ExtractFilterType::Pointer ExtractFilter = ExtractFilterType::New();
	typename Timage4D::SizeType desiredsize = Img->GetLargestPossibleRegion().GetSize();
	desiredsize[3] = 0;
	typename Timage4D::IndexType desiredstart = Img->GetLargestPossibleRegion().GetIndex();
	desiredstart[3] = phase;

	typename Timage4D::RegionType desiredRegion;
	desiredRegion.SetSize(desiredsize);
	desiredRegion.SetIndex(desiredstart);
	ExtractFilter->SetInput(Img);
	ExtractFilter->SetExtractionRegion( desiredRegion );
	ExtractFilter->Update();
	return ExtractFilter->GetOutput();
}

template <typename Timage, typename TimagePointer, typename TimageSize, typename TimageIndex>
TimagePointer IMGIO::ExtractROI(TimagePointer Img, TimageSize ROI_Size, TimageIndex ROI_Start)
{
	typedef itk::ExtractImageFilter< Timage, Timage> ExtractFilterType;
	typename ExtractFilterType::Pointer ExtractFilter = ExtractFilterType::New();

	typename Timage::RegionType desiredRegion;
	desiredRegion.SetSize(ROI_Size);
	desiredRegion.SetIndex(ROI_Start);
	ExtractFilter->SetInput(Img);
	ExtractFilter->SetExtractionRegion( desiredRegion );
	ExtractFilter->Update();
	return ExtractFilter->GetOutput();
}


extern IMGIO ImageIO;



#endif
