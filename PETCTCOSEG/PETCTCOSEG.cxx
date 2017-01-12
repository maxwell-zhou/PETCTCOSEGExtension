#include "stdlib.h"
#include "stdio.h"
#include <itkImageFileReader.h>
#include "ImageType.h"
#include "ImgIO.h"
#include "PreProcess.h"
#include "time.h"
#include "optnet_vce_lib/optnet/_base/array_ref.hxx"
#include "optnet_vce_lib/optnet_graphcut/optnet_gs_gt_multi_dir.hxx"

#include "itkPluginUtilities.h"

#include "PETCTCOSEGCLP.h"

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
#define MAXCOST 10000 

using namespace std;
using namespace optnet;
 

int main(int argc, char* argv[]){
	
	PARSE_ARGS;	
	
	const char * inputCTFile = inputVolume_CT.c_str();
	const char * inputPETFile = inputVolume_PET.c_str();
	int flagMultiSeeds = flag_MultiSeeds;
	int withContext = with_Context;
	int contextCoef = Context_Coef;
	float upThres = up_Thres;
	float lowThres = low_Thres;
	const char * seedOb = inputVolume_OBJ.c_str();
	const char * seedBg = inputVolume_BKG.c_str();	
	
	int numSurf_graphcut = 2;
	
	typedef ImageType3DFLOAT InputImageType;
	typedef ImageType3DCHAR OutputImageType;
	typedef ImageType3DCHAR SeedImageType;
	typedef ImageType3DFLOAT InternalImageType;
	
	InputImageType::Pointer originCTImage, originPETImage;
	originCTImage = ImageIO.LoadImg< InputImageType, InputImageType::Pointer>( inputCTFile );
	originPETImage = ImageIO.LoadImg< InputImageType, InputImageType::Pointer>( inputPETFile );
	
	SeedImageType::Pointer seedImage[2];
	seedImage[0] = ImageIO.LoadImg< SeedImageType, SeedImageType::Pointer>( seedOb );
	seedImage[1] = ImageIO.LoadImg< SeedImageType, SeedImageType::Pointer>( seedBg );
	
	// Configuration configs;
	// ConfigReader config_reader;
	// config_reader.readConfiguration(argv[1],configs);

	

	// string inputCTFile = configs.input_path_name + "InputROI.hdr";
	// string inputPETFile = configs.input_path_name + "PETReg.hdr";
	// int flagMultiSeeds = configs.flag_multi_seeds;
	// cout << "Processing " << inputCTFile << endl;

	
	// int withContext = configs.flag_context;
	//int withContext = 0;
	// int contextCoef = configs.context_coef;
	
	// float upThres = configs.up_thres;
	// float lowThres = configs.low_thres;	
	
	// string seedOb = configs.input_seed_name + "ob.hdr";
	// string seedBg = configs.input_seed_name + "bg.hdr";	

	//Finish image loading
	//Calculate the running time;
	clock_t t1, t2;
	t1 = clock();
	//

	InternalImageType::Pointer scaleCTImage, scalePETImage, smoothImage;
	scaleCTImage = scaleCastImage< InputImageType, InternalImageType >( originCTImage, 255 );
	scalePETImage = scaleCastImage< InputImageType, InternalImageType >( originPETImage, 255 );
	// ImageIO.WriteImg< InternalImageType >( scalePETImage,configs.input_path_name + "scalePETImage.hdr" );
	// ImageIO.WriteImg< InternalImageType >( scaleCTImage,configs.input_path_name + "scaleCTImage.hdr" );
	
	InternalImageType::SizeType imgSize = scaleCTImage->GetLargestPossibleRegion().GetSize();
	cout << "The image size is" << imgSize <<endl; 
		
    //////////////////////////////////	
	typedef optnet_gs_gt_multi_dir<int, long, net_f_xy> OptNet;

    InternalImageType::SizeType CostImgSize;
	CostImgSize[0] = imgSize[0];
	CostImgSize[1] = imgSize[1];
	CostImgSize[2] = imgSize[2];

    OptNet::cost_array_type cost_ob, cost_bg, cost_neigh, cost_context;
	//cost_gs.create( CostImgSize[0], CostImgSize[1],CostImgSize[2], numSurf_graphsearch );
    cost_ob.create( CostImgSize[0], CostImgSize[1],CostImgSize[2], numSurf_graphcut );
	cost_bg.create( CostImgSize[0], CostImgSize[1],CostImgSize[2], numSurf_graphcut );
	cost_neigh.create( CostImgSize[0], CostImgSize[1],CostImgSize[2], numSurf_graphcut );
	cost_context.create( CostImgSize[0], CostImgSize[1],CostImgSize[2], 2 );

	InternalImageType::IndexType index3D;

	InternalImageType::Pointer costCTRegionImage, costPETRegionImage;

	costPETRegionImage = ComputePETRegionCost<InputImageType, SeedImageType>( originPETImage, seedImage[0], upThres, lowThres);
	costCTRegionImage = ComputeRegionCost<InternalImageType, SeedImageType>( scaleCTImage, seedImage[0]);
	
	
	// ImageIO.WriteImg< InternalImageType >( costPETRegionImage,"costPETRegion.hdr" );
	// ImageIO.WriteImg< InternalImageType >( costCTRegionImage,"costCTRegion.hdr" );
   
    //Assign cost	
	cout << "Assigning cost..."<<endl;

	typedef itk::ImageRegionIterator< InternalImageType > IteratorInternalType;
    typedef itk::ImageRegionIterator< SeedImageType > IteratorSeedType;

	IteratorInternalType scaleCTIt( scaleCTImage, scaleCTImage->GetLargestPossibleRegion() );
	IteratorInternalType scalePETIt( scalePETImage, scalePETImage->GetLargestPossibleRegion() );
	IteratorInternalType costCTRegion_It( costCTRegionImage, costCTRegionImage->GetLargestPossibleRegion() );
	IteratorInternalType costPETRegion_It( costPETRegionImage, costPETRegionImage->GetLargestPossibleRegion() );
	IteratorSeedType seedObIt( seedImage[0], seedImage[0]->GetLargestPossibleRegion() );
	IteratorSeedType seedBgIt( seedImage[1], seedImage[1]->GetLargestPossibleRegion() );
	
	scaleCTIt.GoToBegin();
	scalePETIt.GoToBegin();
	costCTRegion_It.GoToBegin();
	costPETRegion_It.GoToBegin();
	seedObIt.GoToBegin();
	seedBgIt.GoToBegin();
    
	int context_cost;
    
	for( index3D[2] = 0;index3D[2] < static_cast<int> (CostImgSize[2]); ++index3D[2] )
	   for ( index3D[1] = 0; index3D[1] < static_cast<int> (CostImgSize[1]); ++index3D[1] )
           for( index3D[0] = 0; index3D[0] < static_cast<int> (CostImgSize[0]); ++index3D[0] )
			  {       
				
				  cost_ob( index3D[0], index3D[1], index3D[2], 0 ) = costCTRegion_It.Get() * 1;
				  cost_ob( index3D[0], index3D[1], index3D[2], 1 ) = costPETRegion_It.Get() * 10;				  
				  
				  cost_bg( index3D[0], index3D[1], index3D[2], 0 ) = (255 - costCTRegion_It.Get()) * 1;				  
				  cost_bg( index3D[0], index3D[1], index3D[2], 1 ) = (255 - costPETRegion_It.Get()) * 10;
				
				  cost_ob( index3D[0], index3D[1], index3D[2], 0 ) = seedObIt.Get() * MAXCOST + cost_ob( index3D[0], index3D[1], index3D[2], 0 ) ;
				  cost_ob( index3D[0], index3D[1], index3D[2], 1 ) = seedObIt.Get() * MAXCOST + cost_ob( index3D[0], index3D[1], index3D[2], 1 ) ;
				  cost_bg( index3D[0], index3D[1], index3D[2], 0 ) = ( 1 - seedBgIt.Get() ) * MAXCOST + cost_bg( index3D[0], index3D[1], index3D[2], 0 );
				  cost_bg( index3D[0], index3D[1], index3D[2], 1 ) = ( 1 - seedBgIt.Get() ) * MAXCOST + cost_bg( index3D[0], index3D[1], index3D[2], 1 );		
						
				  cost_neigh( index3D[0], index3D[1], index3D[2], 0 ) = scaleCTIt.Get();
				  cost_neigh( index3D[0], index3D[1], index3D[2], 1 ) = scalePETIt.Get();
				  
				  context_cost = 0.1*( 255 - abs(costCTRegion_It.Get() - costPETRegion_It.Get()) ) * contextCoef + 250;
				  costCTRegion_It.Set( 255 - abs(costCTRegion_It.Get() - costPETRegion_It.Get()) ); //Output Context term for test.
				  cost_context( index3D[0], index3D[1], index3D[2], 0 ) = context_cost;
				  cost_context( index3D[0], index3D[1], index3D[2], 1 ) = context_cost;
				  	
				  ++scaleCTIt;
				  ++scalePETIt;
				  ++costCTRegion_It;
				  ++costPETRegion_It;
				  ++seedObIt;
				  ++seedBgIt;
				  
			  }
	
	//ImageIO.WriteImg< InternalImageType >( costCTRegionImage,"contextCost.hdr" );
	cout << "Finish cost image assignment" << endl;

    //
    OptNet::net_type resImage( CostImgSize[0], CostImgSize[1], CostImgSize[2], numSurf_graphcut);
    OptNet optnet_graphcut;
	cout << "Create the graph " << endl;
    optnet_graphcut.create( CostImgSize[0], CostImgSize[1],CostImgSize[2], 0, numSurf_graphcut  );
	
	
	
	/////////////////////////////////////////////////////
	/////////////////////////////////////////////////////
    optnet_graphcut.set_ob_cost( cost_ob );
    optnet_graphcut.set_bg_cost( cost_bg );
    optnet_graphcut.set_neigh_cost( cost_neigh );
	long neigh_coef[2];
	neigh_coef[0] = 10000;
	neigh_coef[1] = 1;
	long *p = neigh_coef;
	optnet_graphcut.set_neigh_coef(p);

	if ( withContext == 1 )
	{
		//Set context between graphs
		OptNet::inter_cutcut_type cut_context;
		cut_context.k[0] = 0;
		cut_context.k[1] = 1;
		cut_context.cost_context_cut = &cost_context;
		optnet_graphcut.set_cutcut_relation( cut_context );
	}
		
    optnet_graphcut.solve_all ( resImage, NULL);
    cout<<"solve the graph"<<endl;
    cost_ob.clear();
    cost_bg.clear();
    cost_neigh.clear();
	cost_context.clear();
	
	////////////////////////////////////////////////////////////////
	t2 = clock();
	float diff = ((float)t2 - (float)t1) / CLOCKS_PER_SEC;
	cout << "The running time is " << diff << endl;

    OutputImageType::Pointer resultImage[3];
    typedef itk::CastImageFilter< InternalImageType, OutputImageType > CastType;

	
	for( int i = 0; i < 3; i++ )
	{
	   CastType::Pointer caster = CastType::New();
	   caster->SetInput( scaleCTImage );
	   caster->Update();
	   resultImage[i] = caster->GetOutput();
	   resultImage[i]->FillBuffer( 0 );
	}

	for ( index3D[2] = 0; index3D[2] < static_cast<int>(CostImgSize[2]); ++index3D[2] )
        for ( index3D[1] = 0; index3D[1] < static_cast<int>(CostImgSize[1]); ++index3D[1] )
            for( index3D[0] = 0; index3D[0] < static_cast<int>(CostImgSize[0]); ++index3D[0] )
			{
				  int tmpValue = 0;
				  
				  for ( int i = 0; i < numSurf_graphcut; i++ )
				  {
					  tmpValue =  resImage( index3D[0], index3D[1], index3D[2], i ) * 255;
					  resultImage[i]->SetPixel(index3D, tmpValue );
				  }

				  
			}              

	string resultFileName0, resultFileName1;
	
	stringstream contextCoefString, upThresString, lowThresString;
    contextCoefString << contextCoef;
	upThresString << upThres;
	lowThresString << lowThres;
	

	// if ( withContext == 1 )
	// {
		// resultFileName0 = "SEG_CT_" + contextCoefString.str() + "_" + upThresString.str() + "_" + lowThresString.str() + ".hdr";
		// resultFileName1 = "SEG_PET_" + contextCoefString.str() + "_" + upThresString.str() + "_" + lowThresString.str() + ".hdr";
	// }
	// else
	// {
		// resultFileName0 = "SEG_CT_NC.hdr";
		// resultFileName1 = "SEG_PET_NC_" + upThresString.str() + "_" + lowThresString.str() + ".hdr";
	// }

 
    // ImageIO.WriteImg< OutputImageType >( resultImage[0], resultFileName0 );
	// ImageIO.WriteImg< OutputImageType >( resultImage[1], resultFileName1 );
	
	
	OutputImageType::IndexType seed;
	
	for ( seedObIt.GoToBegin(); !seedObIt.IsAtEnd(); ++seedObIt)
	{
		if ( seedObIt.Get() != 0)
		{
			seed = seedObIt.GetIndex();
			break;
		}
	}
    	
	OutputImageType::Pointer morpImage[2];
		
	morpImage[0] = MorpSmooth< OutputImageType >( resultImage[0], seed, 0, 0, flagMultiSeeds );
	morpImage[1] = MorpSmooth< OutputImageType >( resultImage[0], seed, 1, 1, flagMultiSeeds );


	string morpFileName[2];
/*
	if ( withContext == 1 )
	{
	   morpFileName[0] = "SEG_CT_MP0_" + contextCoefString.str() + "_" + upThresString.str() + "_" + lowThresString.str() + ".hdr";
	   morpFileName[1] = "SEG_CT_MP1_" + contextCoefString.str() + "_" + upThresString.str() + "_" + lowThresString.str() + ".hdr";
	}
	else
	{
	   morpFileName[0] = "SEG_CT_NC_MP0.hdr";
	   morpFileName[1] = "SEG_CT_NC_MP1.hdr";
	}
*/

//	add here
	morpFileName[0] = outputVolume_CTMP0.c_str();
	morpFileName[1] = outputVolume_CTMP1.c_str();
	for ( int i = 0; i < 2; i++)
	{
      
	   ImageIO.WriteImg< OutputImageType >( morpImage[i],morpFileName[i] );
	}
	

   return 0;
   
}
