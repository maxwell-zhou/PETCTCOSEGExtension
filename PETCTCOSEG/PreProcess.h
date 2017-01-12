#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkCastImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkAntiAliasBinaryImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkStatisticsImageFilter.h"
#include "ImageType.h"


using namespace std;

int round(float a) {
return int(a + 0.5);
}

template < class ImageType >
typename ImageType::Pointer scaleImage(typename ImageType::Pointer inputImage)
{
	float fMax,fMin,fScale;
	typedef itk::MinimumMaximumImageCalculator<ImageType> CalculatorType;
	typename CalculatorType::Pointer calculator=CalculatorType::New();
	calculator->SetImage(inputImage);
	calculator->SetRegion(inputImage->GetLargestPossibleRegion());
	calculator->Compute();
	fMax=static_cast<float>(calculator->GetMaximum());
	fMin=static_cast<float>(calculator->GetMinimum());
	fScale=255.0/(fMax-fMin);
	cout<<"Max:"<<fMax<<endl;
	cout<<"Min:"<<fMin<<endl;
	cout<<"Scale:"<<fScale<<endl;
	typedef itk::ImageRegionIterator<ImageType> IteratorType;
	IteratorType inputIt(inputImage, inputImage->GetLargestPossibleRegion());
	for (inputIt.GoToBegin();!inputIt.IsAtEnd();++inputIt){
		//inputIt.Set(255.0-(inputIt.Get()-fMin)*fScale);
		inputIt.Set( 255 - (inputIt.Get()-fMin)*fScale);
	}
	cout<<"Scaling done!"<<endl;
	return inputImage;

}

template < typename TInputImageType, typename TOutputImageType >
typename TOutputImageType::Pointer scaleCastImage( typename TInputImageType::Pointer inputImage, float Scale )
{
	typedef itk::CastImageFilter< TInputImageType, TOutputImageType > CastType;
	typename CastType::Pointer caster = CastType::New();
	caster->SetInput( inputImage );
	caster->Update();
	typename TOutputImageType::Pointer castImage = caster->GetOutput();
	
	float fMax,fMin,fScale;
	typedef itk::MinimumMaximumImageCalculator< TOutputImageType > CalculatorType;
	typename CalculatorType::Pointer calculator = CalculatorType::New();
	calculator->SetImage( castImage );
	calculator->SetRegion( castImage->GetLargestPossibleRegion() );
	calculator->Compute();
	fMax = static_cast<float>(calculator->GetMaximum());
	fMin = static_cast<float>(calculator->GetMinimum());
	fScale = Scale / ( fMax - fMin );
	cout << "Max:" << fMax << endl;
	cout << "Min:" << fMin << endl;
	cout << "Scale:" << fScale << endl;
	
	typedef itk::ImageRegionIterator< TOutputImageType > IteratorType;
	IteratorType castIt( castImage, castImage->GetLargestPossibleRegion() );
	for ( castIt.GoToBegin(); !castIt.IsAtEnd(); ++castIt)
	{		
		castIt.Set( ( castIt.Get() - fMin ) * fScale );
	}
	cout << "Scaling done!" << endl;

	return castImage;

}

template < class ImageType >
typename ImageType::Pointer SmoothImage( typename ImageType::Pointer inputImage )
{
	typedef itk::GradientAnisotropicDiffusionImageFilter< ImageType, ImageType > SmoothFilterType;
	typename SmoothFilterType::Pointer smoothFilter = SmoothFilterType::New();
	smoothFilter->SetInput( inputImage );
	smoothFilter->SetNumberOfIterations( 5 );
	smoothFilter->SetTimeStep( 0.125 );
	smoothFilter->SetConductanceParameter( 9 );
	smoothFilter->Update();
	return smoothFilter->GetOutput();
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

template < typename TInputImageType >
typename TInputImageType::Pointer antiAlias( typename TInputImageType::Pointer inputImage )
{
	
	typedef itk::AntiAliasBinaryImageFilter< TInputImageType, ImageType3DFLOAT > AntiAliasFilterType;
	typename AntiAliasFilterType::Pointer antiAliasFilter = AntiAliasFilterType::New();
	antiAliasFilter->SetInput( inputImage );
	antiAliasFilter->SetMaximumRMSError( 0.001 );
	antiAliasFilter->Update();
	//return antiAliasFilter->GetOutput();

	ImageType3DFLOAT::Pointer internalImage = antiAliasFilter->GetOutput();
	typedef itk::ImageRegionIterator< ImageType3DFLOAT > IteratorType;
	IteratorType internalIt( internalImage, internalImage->GetLargestPossibleRegion() );
	for ( internalIt.GoToBegin(); !internalIt.IsAtEnd(); ++internalIt)
	{
	    if ( internalIt.Get() >= 0 )
			internalIt.Set( 255 );
		else
			internalIt.Set( 0 );
	}
	typedef itk::CastImageFilter< ImageType3DFLOAT, TInputImageType> CastImageType;
	typename CastImageType::Pointer castFilter=CastImageType::New();
	castFilter->SetInput( internalImage );
	castFilter->Update();
	cout<<"AntiAliasing done!"<<endl;
	return castFilter->GetOutput();
	
}
template <typename TInputImageType, typename TObImageType >
typename TInputImageType::Pointer ComputePETRegionCost( typename TInputImageType::Pointer inputImage, typename TObImageType::Pointer obImage, float upThres, float lowThres)
{

	ImageType3DFLOAT::Pointer costImage;

	
	typedef itk::CastImageFilter< TInputImageType, ImageType3DFLOAT > CastFLOATType;
	typename CastFLOATType::Pointer casterFLOAT = CastFLOATType::New();
	casterFLOAT->SetInput( inputImage );
	casterFLOAT->Update();
	costImage = casterFLOAT->GetOutput();

	//cout<<"Scaling cost..."<<endl;
	float fMax, fMin, fScale;

	typedef itk::MinimumMaximumImageCalculator<ImageType3DFLOAT> CalculatorType;
	typename CalculatorType::Pointer cal = CalculatorType::New();
	cal->SetImage( costImage );
	cal->SetRegion( costImage->GetLargestPossibleRegion() );
	cal->Compute();
	fMax = static_cast<float>(cal->GetMaximum());
	fMin = static_cast<float>(cal->GetMinimum());
	fScale = 1.0 / (fMax - fMin);
	
	////////////////////////////////
	typedef itk::StatisticsImageFilter<ImageType3DFLOAT> StatisticsFilterType;
	typename StatisticsFilterType::Pointer statisticsFilter = StatisticsFilterType::New();
	statisticsFilter->SetInput( costImage );
	statisticsFilter->Update();
	statisticsFilter->GetMean();
	float fmean, fmax, fmin, fstd, fcof;
    fmean = statisticsFilter->GetMean();
	fmax = statisticsFilter->GetMaximum();
	fmin = statisticsFilter->GetMinimum();
	fstd = sqrt(statisticsFilter->GetVariance());
	fcof = (fmean+3*fstd-fmin)/(fmax-fmin);
    
	cout << "The computed mean is " << statisticsFilter->GetMean() << endl;
	cout << "The computed maximum is " << statisticsFilter->GetMaximum() << endl;
	cout << "The computed variance is " << statisticsFilter->GetVariance() << endl;
	cout << "The coefficient is " << fcof<< endl;
	
	typedef itk::ImageRegionIterator< TObImageType > IteratorObType;
	IteratorObType obIt( obImage, obImage->GetLargestPossibleRegion() );
	//cout << "The maximum pet value is " << fMax << endl;
	IteratorType3DFLOAT costIt( costImage, costImage->GetLargestPossibleRegion() );
	float petMax = 0;
	for ( obIt.GoToBegin(), costIt.GoToBegin(); !costIt.IsAtEnd(); ++costIt, ++obIt)
	{
		if ( obIt.Get() == 1)
		{
			if (costIt.Get() > petMax)
			petMax = costIt.Get();
		}
		costIt.Set( ( costIt.Get() - fMin ) * fScale );
	}
	
	//float b = 0.4 * (petMax/fMax);
	float b = upThres;
	float p = 3;
	//float a = 1000 / fMax;
	float a = lowThres;
	//float sigma_1 = 0.5, sigma_2 = 0.5, sigma_3 = 0.5; 
	float alpha = 1;
	float tmpCost;
	cout << "The value of a is " << a << endl;
	for ( costIt.GoToBegin(); !costIt.IsAtEnd(); ++costIt)
	{

		tmpCost = costIt.Get();
        
		if (tmpCost > b)
		costIt.Set(255);
		else if (tmpCost < a)
		costIt.Set(0);
		else
		{
			/*
			tmpCost = (tmpCost - a)/(b - a) - sigma_2;
			if (tmpCost > 0)
				tmpCost = pow(sigma_1, (p-1)/p) * pow(tmpCost, 1/p) + sigma_3;
			else if (tmpCost == 0)
				tmpCost = sigma_3;
			else
			tmpCost = -pow(sigma_1, (p-1)/p) * pow( -tmpCost, 1/p ) + sigma_3;
			*/
			
	///// replaced by Junjie to try!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
			//tmpCost = (tmpCost - a)/(b - a) - 0.5;
			tmpCost = (tmpCost - a)/(b - a) - 0.0;
			tmpCost = 1 / ( 1 + exp( -tmpCost/alpha ) ) * 255;
			costIt.Set( tmpCost );
		}
	}
	
	//fScale = 255.0 / (fMax - fMin);
	
	typedef itk::CastImageFilter<ImageType3DFLOAT, TInputImageType> CastImageType;
	typename CastImageType::Pointer castFilter = CastImageType::New();
	castFilter->SetInput( costImage );
	castFilter->Update();
	return castFilter->GetOutput();
	
}

template <typename TInputImageType, typename TObImageType >
typename TInputImageType::Pointer ComputeRegionCost( typename TInputImageType::Pointer inputImage, typename TObImageType::Pointer obImage)
{
	typedef itk::ImageRegionIterator< TInputImageType > IteratorInputType;
	typedef itk::ImageRegionIterator< TObImageType > IteratorObType;

	IteratorInputType inputIt( inputImage, inputImage->GetLargestPossibleRegion() );
	IteratorObType obIt( obImage, obImage->GetLargestPossibleRegion() );

	int sum = 0, sumSq = 0, numPt = 0;
	for ( inputIt.GoToBegin(),obIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++obIt)
	{
		if ( obIt.Get() == 1)
		{
			sum = sum + inputIt.Get();
			sumSq = sumSq + inputIt.Get() * inputIt.Get();
			numPt++;
		}
			
	}

	float mean, var;
	mean = static_cast<float> (sum) / numPt;
	var =  1.0 / numPt  * ( sumSq ) - mean * mean;
	std::cout << "mean is " << mean << std::endl;
	std::cout << "var is " << var << std::endl;

	ImageType3DFLOAT::Pointer costImage;

	typedef itk::CastImageFilter< TInputImageType, ImageType3DFLOAT > CastFLOATType;
	typename CastFLOATType::Pointer casterFLOAT = CastFLOATType::New();
	casterFLOAT->SetInput( inputImage );
	casterFLOAT->Update();
	costImage = casterFLOAT->GetOutput();

	IteratorType3DFLOAT costIt( costImage, costImage->GetLargestPossibleRegion() );

	for ( inputIt.GoToBegin(),costIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++costIt)
	{
	
		float cost = exp ( -( inputIt.Get() - mean ) * ( inputIt.Get() - mean ) / ( 2 * var ) );  
		//float cost = ( inputIt.Get() - mean ) * ( inputIt.Get() - mean ) / ( 2 * var ) ;
		costIt.Set( cost );
	}

	//cout<<"Scaling cost..."<<endl;
	float fMax, fMin, fScale;

	typedef itk::MinimumMaximumImageCalculator<ImageType3DFLOAT> CalculatorType;
	
	typename CalculatorType::Pointer cal = CalculatorType::New();
	cal->SetImage( costImage );
	cal->SetRegion( costImage->GetLargestPossibleRegion() );
	cal->Compute();
	fMax = static_cast<float>(cal->GetMaximum());
	fMin = static_cast<float>(cal->GetMinimum());
	fScale = 255.0 / (fMax - fMin);
	

	for ( costIt.GoToBegin(); !costIt.IsAtEnd(); ++costIt)
	{
		costIt.Set( ( costIt.Get() - fMin ) * fScale );
	}

	typedef itk::CastImageFilter<ImageType3DFLOAT, TInputImageType> CastImageType;
	typename CastImageType::Pointer castFilter = CastImageType::New();
	castFilter->SetInput( costImage );
	castFilter->Update();
	return castFilter->GetOutput();
	
}

template < typename TInputImageType >
typename TInputImageType::Pointer ConnectThres( typename TInputImageType::Pointer inputImage, typename TInputImageType::IndexType& index )
{
	typedef itk::ConnectedThresholdImageFilter< TInputImageType, TInputImageType > ConnectedFilterType;
    typename ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
	connectedThreshold->SetInput( inputImage );
	connectedThreshold->SetLower(  254  );
    connectedThreshold->SetUpper(  255  );
	connectedThreshold->SetReplaceValue( 255 );
	connectedThreshold->SetSeed( index );
	connectedThreshold->Update();
	return connectedThreshold->GetOutput();
	
}
template < typename TInputImageType >
typename TInputImageType::Pointer MorpSmooth( typename TInputImageType::Pointer inputImage, typename TInputImageType::IndexType& index, typename TInputImageType::PixelType TubeRadius1, typename TInputImageType::PixelType TubeRadius2, int flagNoConnected )
{
    typename TInputImageType::Pointer tmpImage;
	if (flagNoConnected == 1)
	tmpImage = inputImage;
	else
	{
		typedef itk::ConnectedThresholdImageFilter< TInputImageType, TInputImageType > ConnectedFilterType;
		typename ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
		//connectedThreshold->SetInput( erodeImageFilter->GetOutput() );
		connectedThreshold->SetInput( inputImage );
		connectedThreshold->SetLower(  254  );
		connectedThreshold->SetUpper(  255  );
		connectedThreshold->SetReplaceValue( 255 );
		connectedThreshold->SetSeed( index );
		connectedThreshold->Update();
		tmpImage = connectedThreshold->GetOutput();
	}

	typedef itk::BinaryBallStructuringElement< float, 3> StructuringElementType;

    StructuringElementType structuringElement1;
    structuringElement1.SetRadius( TubeRadius1 );
    structuringElement1.CreateStructuringElement();

	typedef itk::BinaryErodeImageFilter< TInputImageType, TInputImageType, StructuringElementType > ErodeImageFilterType;
    typename ErodeImageFilterType::Pointer erodeImageFilter = ErodeImageFilterType::New();
	//erodeImageFilter->SetInput( inputImage );
	erodeImageFilter->SetInput( tmpImage );
	erodeImageFilter->SetKernel( structuringElement1 );
	erodeImageFilter->SetErodeValue( 255 );
	erodeImageFilter->Update();

    /*
	//temporary test
	index[0]=42; index[1]=44; index[2]=44;
    typename ConnectedFilterType::Pointer connectedThreshold2 = ConnectedFilterType::New();
	connectedThreshold2->SetInput( erodeImageFilter->GetOutput() );
	//connectedThreshold->SetInput( inputImage );
	connectedThreshold2->SetLower(  254  );
    connectedThreshold2->SetUpper(  255  );
	connectedThreshold2->SetReplaceValue( 255 );
	connectedThreshold2->SetSeed( index );
	connectedThreshold2->Update();
	*/


    StructuringElementType structuringElement2;
    structuringElement2.SetRadius( TubeRadius2 );
    structuringElement2.CreateStructuringElement();

	typedef itk::BinaryDilateImageFilter< TInputImageType, TInputImageType, StructuringElementType > DilateImageFilterType;
    typename DilateImageFilterType::Pointer dilateImageFilter = DilateImageFilterType::New();
	//dilateImageFilter->SetInput( connectedThreshold2->GetOutput() );
	dilateImageFilter->SetInput( erodeImageFilter->GetOutput() );
	dilateImageFilter->SetKernel( structuringElement2 );
	dilateImageFilter->SetDilateValue( 255 );
	dilateImageFilter->Update();

	return dilateImageFilter->GetOutput();
	
}