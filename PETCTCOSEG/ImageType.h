
#ifndef IMGAGETYPE_H
#define IMGAGETYPE_H


#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkRGBPixel.h"

typedef itk::Image<short, 2 > ImageType2D;
typedef itk::Image<unsigned char, 2 > ImageType2DCHAR;
typedef itk::Image<float, 2 > ImageType2DFLOAT;

//Color image
typedef itk::RGBPixel< unsigned char >    PixelType;
typedef itk::Image< PixelType, 2 >   ImageType2DCOLOR;
typedef itk::Image< PixelType, 3 >   ImageType3DCOLOR;

//
typedef itk::Image<short, 3 > ImageType3D;
typedef itk::Image<unsigned char, 3 > ImageType3DCHAR;
typedef itk::Image<float, 3 > ImageType3DFLOAT;
typedef itk::Image<double, 3 > ImageType3DDOUBLE;
typedef itk::Image<short, 4 > ImageType4D;
typedef itk::Image<unsigned char, 4 > ImageType4DCHAR;
typedef itk::Image<float, 4 > ImageType4DFLOAT;
typedef itk::Image<double, 4 > ImageType4DDOUBLE;

typedef itk::ImageRegionIterator< ImageType2D> IteratorType2D;
typedef itk::ImageRegionIterator< ImageType2DCHAR> IteratorType2DCHAR;
typedef itk::ImageRegionIterator< ImageType2DFLOAT> IteratorType2DFLOAT;
typedef itk::ImageRegionIterator< ImageType2DCOLOR> IteratorType2DCOLOR;

typedef itk::ImageRegionIterator< ImageType3D> IteratorType3D;
typedef itk::ImageRegionIterator< ImageType3DCHAR> IteratorType3DCHAR;
typedef itk::ImageRegionIterator< ImageType3DFLOAT> IteratorType3DFLOAT;
typedef itk::ImageRegionIterator< ImageType4D> IteratorType4D;
typedef itk::ImageRegionIterator< ImageType4DCHAR> IteratorType4DCHAR;
typedef itk::ImageRegionIterator< ImageType4DFLOAT> IteratorType4DFLOAT;

typedef itk::Point<float, 3> PointType3DFloat;
typedef itk::Vector<float, 3> VectorType3DFloat;
#endif	