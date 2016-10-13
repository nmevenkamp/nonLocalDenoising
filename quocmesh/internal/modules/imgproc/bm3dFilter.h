#ifndef __BM3DFILTER_H
#define __BM3DFILTER_H

#include "neighborhoodFilter.h"
#include "blockStructures.h"
#include "blockOperators.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace im {

class BM3D_SIMILARITYSEARCH_METHOD : public SIMILARITYSEARCH_METHOD {
public:
  static std::string getIdentifier ( const int SimilaritySearchMethod ) {
    return SIMILARITYSEARCH_METHOD::getIdentifier ( SimilaritySearchMethod );
  }
  
  static std::string toString ( const int SimilaritySearchMethod ) {
    return SIMILARITYSEARCH_METHOD::toString ( SimilaritySearchMethod );
  }
};

class BM3D_PROFILE {
public:
  static const int NP = 0;      // Normal Profile (default, balanced quality)
  static const int LC = 1;      // Low Complexity Profile (fast, lower quality)
  static const int HIGH = 2;    // High Profile (high quality, not documented in [1])
  
  static const int NUM = 3;
  
  static std::string getIdentifier ( const int Profile ) {
    if ( Profile == NP ) return "np";
    else if ( Profile == LC ) return "lc";
    else if ( Profile == HIGH ) return "high";
    else throw aol::Exception ( "Did not recognize profile!", __FILE__, __LINE__ );
  }
  
  static std::string toString ( const int Profile ) {
    if ( Profile == NP ) return "Normal";
    else if ( Profile == LC ) return "Low complexity";
    else if ( Profile == HIGH ) return "High quality";
    else throw aol::Exception ( "Did not recognize profile!", __FILE__, __LINE__ );
  }
};

class BM3D_ESTIMATION_METHOD {
public:
  static const int HT = 0;      // Parameters are set for hard-thresholding (HT) in 2D and 3D transform domain
  static const int WIENER = 1;  // Parameters are set for Wiener filtering in 3D transform domain
  
  static const int NUM = 2;
  
  static std::string getIdentifier ( const int EstimationMethod ) {
    if ( EstimationMethod == HT ) return "ht";
    if ( EstimationMethod == WIENER ) return "wie";
    else throw aol::Exception ( "Did not recognize estimation method!", __FILE__, __LINE__ );
  }
  
  static std::string toString ( const int EstimationMethod ) {
    if ( EstimationMethod == HT ) return "Hard-thresholding";
    if ( EstimationMethod == WIENER ) return "Wiener filtering";
    else throw aol::Exception ( "Did not recognize estimation method!", __FILE__, __LINE__ );
  }
};

template <typename _RealType, typename _PictureType, typename _BlockType, typename _BlockSizeType>
struct BaseBM3DOptions : public NHFilterOptions<_RealType, _PictureType, _BlockType, _BlockSizeType> {
  int profile;
  int method;
  int N2;
  int refStep;
  int denoiseInTransformDomain;
  bool restrictStackSizeToPowersOfTwo;
  int levelSetNumPointsPerSample, levelSetOverSamplingFactor;
  _RealType filterParameter;
  
  BaseBM3DOptions ( const std::string OutputDir = "", const bool Quiet = true )
    : NHFilterOptions<_RealType, _PictureType, _BlockType, _BlockSizeType> ( OutputDir, Quiet ),
      profile ( BM3D_PROFILE::NP ), method ( BM3D_ESTIMATION_METHOD::HT ), N2 ( 0 ), refStep ( 0 ), denoiseInTransformDomain ( true ),
      restrictStackSizeToPowersOfTwo ( true ),
      levelSetNumPointsPerSample ( 64 ), levelSetOverSamplingFactor ( 1 ),
      filterParameter ( 0 ) { }
  
  BaseBM3DOptions ( const BaseBM3DOptions<_RealType, _PictureType, _BlockType, _BlockSizeType> &Options )
    : NHFilterOptions<_RealType, _PictureType, _BlockType, _BlockSizeType> ( Options ),
      profile ( Options.profile ), method ( Options.method ), N2 ( Options.N2 ), refStep ( Options.refStep ), denoiseInTransformDomain ( Options.denoiseInTransformDomain ),
      restrictStackSizeToPowersOfTwo ( Options.restrictStackSizeToPowersOfTwo ),
      levelSetNumPointsPerSample ( Options.levelSetNumPointsPerSample ), levelSetOverSamplingFactor ( Options.levelSetOverSamplingFactor ),
      filterParameter ( Options.filterParameter ) { }
};

template <typename _RealType, qc::Dimension Dim = qc::QC_2D, typename _PictureType = qc::ScalarArray<_RealType, Dim>,
          typename _BlockType = qc::ScalarArray<_RealType, qc::QC_2D>, typename _BlockSizeType = qc::GridSize<qc::QC_2D> >
struct BM3DOptions;

template <typename _RealType, qc::Dimension Dim>
struct BM3DOptions<_RealType, Dim, qc::ScalarArray<_RealType, Dim>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D> >
  : public BaseBM3DOptions<_RealType, qc::ScalarArray<_RealType, Dim>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D> > {
  typedef _RealType RealType;
  typedef qc::ScalarArray<_RealType, Dim> PictureType;
  typedef qc::ScalarArray<_RealType, qc::QC_2D> BlockType;
  typedef qc::GridSize<qc::QC_2D> BlockSizeType;
  
  BM3DOptions ( const std::string OutputDir = "", const bool Quiet = true )
    : BaseBM3DOptions<RealType, PictureType, BlockType, BlockSizeType> ( OutputDir, Quiet ) {
    this->blockSize.setAll ( 8 );
  }
  
  BM3DOptions ( const BM3DOptions<RealType, Dim, PictureType, BlockType, BlockSizeType> &Options )
    : BaseBM3DOptions<RealType, PictureType, BlockType, BlockSizeType> ( Options ) { }
};

template<typename _RealType>
struct BM3DOptions<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_3D>, aol::Vector<_RealType>, int>
  : public BaseBM3DOptions<_RealType, qc::ScalarArray<_RealType, qc::QC_3D>, aol::Vector<_RealType>, int> {
  typedef _RealType RealType;
  const static qc::Dimension Dim = qc::QC_2D;
  typedef qc::ScalarArray<_RealType, qc::QC_3D> PictureType;
  typedef aol::Vector<_RealType> BlockType;
  typedef int BlockSizeType;
    
  BM3DOptions ( const std::string OutputDir = "", const bool Quiet = true )
    : BaseBM3DOptions<RealType, PictureType, BlockType, BlockSizeType> ( OutputDir, Quiet ) { }
  
  BM3DOptions ( const BM3DOptions<RealType, Dim, PictureType, BlockType, BlockSizeType> &Options )
    : BaseBM3DOptions<RealType, PictureType, BlockType, BlockSizeType> ( Options ) { }
};

template<typename _RealType>
struct BM3DOptions<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_3D>, qc::ScalarArray<_RealType, qc::QC_3D>, qc::GridSize<qc::QC_3D> >
: public BaseBM3DOptions<_RealType, qc::ScalarArray<_RealType, qc::QC_3D>, qc::ScalarArray<_RealType, qc::QC_3D>, qc::GridSize<qc::QC_3D> > {
  typedef _RealType RealType;
  const static qc::Dimension Dim = qc::QC_2D;
  typedef qc::ScalarArray<_RealType, qc::QC_3D> PictureType;
  typedef qc::ScalarArray<_RealType, qc::QC_3D> BlockType;
  typedef qc::GridSize<qc::QC_3D> BlockSizeType;
  
  BM3DOptions ( const std::string OutputDir = "", const bool Quiet = true )
    : BaseBM3DOptions<RealType, PictureType, BlockType, BlockSizeType> ( OutputDir, Quiet ) {
    this->blockSize[1] = 8;
    this->blockSize[2] = 8;
  }
  
  BM3DOptions ( const BM3DOptions<RealType, Dim, PictureType, BlockType, BlockSizeType> &Options )
    : BaseBM3DOptions<RealType, PictureType, BlockType, BlockSizeType> ( Options ) { }
};


struct BaseBM3DSharedDataStructure {
#ifdef _OPENMP
  const int numThreadsMax;
  
  BaseBM3DSharedDataStructure ( ) : numThreadsMax ( omp_get_max_threads ( ) ) { }
#endif
};

template <typename RealType, qc::Dimension Dim, typename PictureType>
struct BM3DSharedDataStructure : public BaseBM3DSharedDataStructure {
  PictureType eBuff, wBuff;
  qc::ScalarArray<unsigned int, Dim> numAggregates;
  const int numBufferEntries, numAggregatesEntries;
#ifdef _OPENMP
  aol::RandomAccessContainer<PictureType> eBuffsThreads, wBuffsThreads;
  aol::RandomAccessContainer<qc::ScalarArray<unsigned int, Dim> > numAggregatesThreads;
#endif
  
  BM3DSharedDataStructure ( const PictureType &Input, const qc::CoordType &NX )
    : eBuff ( Input, aol::STRUCT_COPY ), wBuff ( Input, aol::STRUCT_COPY ), numAggregates ( NX ),
      numBufferEntries ( eBuff.size ( ) ), numAggregatesEntries ( numAggregates.size ( ) )
#ifdef _OPENMP
      , eBuffsThreads ( numThreadsMax, eBuff ), wBuffsThreads ( numThreadsMax, wBuff ), numAggregatesThreads ( numThreadsMax, numAggregates )
#endif
    { }
};


template <typename RealType>
struct BaseBM3DPrivateDataStructure {
  aol::RandomAccessContainer<qc::CoordType> Sx;
  aol::Vector<RealType> blockDistances;
  int Nz;
  RealType weight;
#ifdef _OPENMP
  const int numThreads;
  const int thrIdx;
#endif
  
  BaseBM3DPrivateDataStructure ( const int NumNodes )
    : Sx ( NumNodes ), blockDistances ( NumNodes ), Nz ( 0 ), weight ( 0 )
#ifdef _OPENMP
      , numThreads ( omp_get_num_threads ( ) ), thrIdx ( omp_get_thread_num ( ) )
#endif
      { }
};

template <typename RealType, typename PictureType, typename BlockType, typename BlockSizeType, typename BlockStackType, typename BlockCollectionType>
struct BM3DPrivateDataStructure : public BaseBM3DPrivateDataStructure<RealType> {
  BlockCollectionType blocksInitial, blocksEstimate;
  BlockCollectionType *blocksMethod;
  HTBlockStackDenoisingOp<RealType, BlockType, BlockSizeType> htOp;
  WienerBlockStackDenoisingOp<RealType, BlockType, BlockSizeType> wienerOp;
  BlockStackType initialBlockStack, estimateBlockStack, denoisedBlockStack;
  
  BM3DPrivateDataStructure ( const int Method, const BlockSizeType &BlockSize, const int NumNodes, const int N2,
                             const bool DenoiseInTransformDomain,
                             const std::string &TransformXY, const std::string &TransformZ, const RealType HTThreshold, const RealType StdDev )
    : BaseBM3DPrivateDataStructure<RealType> ( NumNodes ),
      blocksMethod ( Method == BM3D_ESTIMATION_METHOD::HT ? &blocksInitial : &blocksEstimate ),
      initialBlockStack ( N2, BlockSize ), estimateBlockStack ( N2, BlockSize ), denoisedBlockStack ( N2, BlockSize ) {
    if ( DenoiseInTransformDomain ) {
      htOp.setTransforms ( BlockSize, TransformXY, N2, TransformZ );
      htOp.setThreshold ( HTThreshold );
      wienerOp.setTransforms ( BlockSize, TransformXY, N2, TransformZ );
      wienerOp.setNoiseStdDev ( StdDev );
    }
  }
};


template <typename _RealType, qc::Dimension Dim, typename _PictureType, typename _BlockType, typename _BlockSizeType>
class BM3DFilterTrait {
public:
  typedef BM3DOptions<_RealType, Dim, _PictureType, _BlockType, _BlockSizeType> OptionsType;
  typedef BlockStack<_RealType, _BlockType, _BlockSizeType> BlockStackType;
  typedef BlockCollection<_RealType, Dim, _PictureType, _BlockType, _BlockSizeType> BlockCollectionType;
  typedef BM3D_SIMILARITYSEARCH_METHOD SimilaritySearchMethodType;
  
  typedef BM3DSharedDataStructure<_RealType, Dim, _PictureType> SharedDataStructureType;
  typedef BM3DPrivateDataStructure<_RealType, _PictureType, _BlockType, _BlockSizeType, BlockStackType, BlockCollectionType> PrivateDataStructureType;
};


//  The following class BM3DFilter has been implemented by
//
//    Niklas Mevenkamp, email: mevenkamp@aices.rwth-aachen.de
//
//  based on the description of the BM3D algorithm as published in:
//
//      K. Dabov, A. Foi, V. Katkovnik, and K. Egiazarian, "Image Denoising
//      by Sparse 3D Transform-Domain Collaborative Filtering,"
//      IEEE Transactions on Image Processing, vol. 16, no. 8, August, 2007.
//      preprint at http://www.cs.tut.fi/~foi/GCF-BM3D
//
//  AUTHORS:
//      Kostadin Dabov, email: dabov _at_ cs.tut.fi
/**
 * \brief Abstract implementation of the BM3D filter (needs specialization depending on dimension of input image, node space, blocks)
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, qc::Dimension Dim, typename _PictureType, typename _BlockType, typename _BlockSizeType, typename FilterTrait>
class BaseBM3DFilter : public CollaborativeNeighborhoodFilter<_RealType, Dim, _PictureType, _BlockType, _BlockSizeType, FilterTrait> {
public:
  typedef _RealType RealType;
  typedef _PictureType PictureType;
  typedef _BlockType BlockType;
  typedef _BlockSizeType BlockSizeType;
  
  typedef qc::MultiArray<RealType, Dim, 3> ColoredPictureType;
 
  typedef typename FilterTrait::OptionsType OptionsType;
  typedef typename FilterTrait::BlockStackType BlockStackType;
  typedef typename FilterTrait::BlockCollectionType BlockCollectionType;
  typedef typename FilterTrait::SimilaritySearchMethodType SimilaritySearchMethodType;
  
  typedef typename FilterTrait::SharedDataStructureType SharedDataStructureType;
  typedef typename FilterTrait::PrivateDataStructureType PrivateDataStructureType;
protected:
  int _N2;                                 // Maximum number of similar blocks for 3D stacks (must be a power of 2, if wavelet denoising is used along 3rd dimension)
  RealType _tauMatch;                      // Threshold for the block distance (d-distance)
  std::string _transformXY, _transformZ;   // ID of the block and stack transform
  RealType _htThreshold;                   // Parameter of hard-thresholding operator
  RealType _filterParameter, _filterParameterSqr, _stdDevSqr; // Parameters used in case of NLM-like denoising (no transform domain)
public:
  BaseBM3DFilter ( ) {
    this->_blockAnchor.set ( 0, 0 );
  }
  
  void apply ( const OptionsType &Options, const PictureType &Arg, PictureType &Dest ) {
    this->_input = &Arg;
    this->_estimate = &Dest;

    this->setCtrlCHandler();
    
    // Preprocess
    OptionsType options ( Options );
    initializeOnce ( options );
    this->preprocess ( options );
    
    // Denoise
    options.method = BM3D_ESTIMATION_METHOD::HT;
    if ( !this->wantsInterrupt ( ) ) denoise ( options );
    if ( Options.denoiseInTransformDomain ) {
      options.method = BM3D_ESTIMATION_METHOD::WIENER;
      if ( !this->wantsInterrupt ( ) ) denoise ( options );
    }

    // Postprocess
    this->postprocess ( options );
    
    // (Optional) console output
    if ( !Options.quietMode
      && options.groundTruth != NULL
      && (*options.groundTruth).getNumX ( ) == this->_estimate->getNumX ( ) && (*options.groundTruth).getNumY ( ) == this->_estimate->getNumY ( ) ) {
      if ( options.noiseType == NOISE_TYPE::GAUSSIAN ) std::cerr << "PSNR = " << aol::PSNR<RealType> ( *options.groundTruth, *(this->_estimate), 255 ) << " dB" << std::endl;
      else if ( options.noiseType == NOISE_TYPE::POISSON ) std::cerr << "PSNR = " << aol::PSNR<RealType> ( *options.groundTruth, *(this->_estimate) ) << " dB" << std::endl;
    }
    
    this->unsetCtrlCHandler();
  }
  
  void apply ( const OptionsType &Options, const qc::CoordType &XRef, const PictureType &Arg, BlockType &Dest ) {
    OptionsType options ( Options );
    options.method = BM3D_ESTIMATION_METHOD::HT;
    BlockCollectionType blocksInitial, blocksEstimate;
    prepareGetBlockStack ( options, Arg, blocksInitial, blocksEstimate );
    
    // Process the specified reference block
    HTBlockStackDenoisingOp<RealType, BlockType, BlockSizeType> htOp ( this->_blockSize, _transformXY, _N2, _transformZ, _htThreshold );
    WienerBlockStackDenoisingOp<RealType, BlockType, BlockSizeType> wienerOp ( this->_blockSize, _transformXY, _N2, _transformZ, this->_stdDev );
    aol::RandomAccessContainer<qc::CoordType> Sx ( this->_numNodes );
    aol::Vector<RealType> blockDistances ( this->_numNodes );
    BlockStackType initialBlockStack ( _N2, this->_blockSize ), estimateBlockStack ( _N2, this->_blockSize ), denoisedBlockStack ( _N2, this->_blockSize );
    int Nz = 0;
    RealType weight = 0;
    
    blockMatching ( Options, blocksInitial, XRef, Sx, blockDistances, Nz );
    blockStacking ( Options, blocksInitial, blocksEstimate, Sx, Nz, initialBlockStack, estimateBlockStack );
    blockStackDenoising ( Options, htOp, wienerOp, blockDistances, Nz, initialBlockStack, estimateBlockStack, denoisedBlockStack, weight );
    
    options.method = BM3D_ESTIMATION_METHOD::WIENER;
    initializeIteration ( options );
    for ( int z=0; z<Nz ; ++z ) estimateBlockStack[z] = denoisedBlockStack[z];
    blockStackDenoising ( Options, htOp, wienerOp, blockDistances, Nz, initialBlockStack, estimateBlockStack, denoisedBlockStack, weight );
    
    Dest.reallocate ( this->_blockSize );
    Dest = denoisedBlockStack[0];
  }
  
  void setSimilaritySearchGrid ( const OptionsType &Options, const PictureType &Arg, const qc::CoordType &XRef,
                                 aol::RandomAccessContainer<qc::CoordType> &CoordsFinal, aol::RandomAccessContainer<qc::CoordType> &CoordsTmp ) {
    // Preprocess
    OptionsType options ( Options );
    options.method = BM3D_ESTIMATION_METHOD::HT;
    BlockCollectionType blocksInitial, blocksEstimate;
    prepareGetBlockStack ( options, Arg, blocksInitial, blocksEstimate );
    
    // Perform block-matching and in output PNG mark all corners of iterated blocks in the search region red
    RealType dist;
    for ( SimilaritySearchIterator<RealType, Dim>* searchIt = this->getSimilaritySearchIterator ( Options, XRef ); searchIt->notAtEnd ( ) ; ++(*searchIt) ) {
      this->setBlockDistance ( blocksInitial, dist, XRef, *(*searchIt) );
      
      if ( searchIt->isCurFinal ( ) ) CoordsFinal.pushBack ( *(*searchIt) );
      else CoordsTmp.pushBack ( *(*searchIt) );
      
      searchIt->update ( dist );
    }
  }
  
  void prepareGetBlockStack ( OptionsType &Options, const PictureType &Arg,
                              BlockCollectionType &BlocksInitial, BlockCollectionType &BlocksEstimate ) {
    this->_input = &Arg;
    initializeOnce ( Options );
    this->preprocess ( Options );
    initializeIteration ( Options );
    fillBlocksFromPreviousEstimate ( Options, BlocksInitial, BlocksEstimate );
  }
  
  void getBlockStack ( const OptionsType &Options, const qc::CoordType &XRef,
                       BlockCollectionType &Blocks,
                       BlockStackType &BlockStack, aol::RandomAccessContainer<qc::CoordType > &Sx, aol::Vector<RealType> &BlockDistances,
                       int &Nz ) const {
    blockMatching ( Options, Blocks, XRef, Sx, BlockDistances, Nz );
    blockStacking ( Options, Blocks, Sx, Nz, BlockStack );
  }
  
  void getWeightedMeanBlock ( const BlockStackType &BlockStack, const aol::Vector<RealType> &BlockDistances, const int Nz, BlockType &Mean ) const {
    Mean.setZero ( );
    RealType normFactor = 0;
    for ( int z=0; z<Nz ; ++z ) {
      const RealType weight = exp ( -aol::Max<RealType> ( BlockDistances[z] / this->_blockNumPixels - 2 * _stdDevSqr, 0.0 ) / _filterParameterSqr  );
      Mean.addMultiple ( BlockStack[z], weight );
      normFactor += weight;
    }
    Mean /= normFactor;
  }
  
  static std::string getMethodIdentifier ( const OptionsType &Options ) {
    return getMethodPrefix ( Options ) + "-" + name ( );
  }
  
  static std::string getMethodPrefix ( const OptionsType &Options ) {
    std::string methodPrefix = BM3D_PROFILE::getIdentifier ( Options.profile );
    methodPrefix += "-" + NeighborhoodFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::getMethodPrefix ( Options );
    return methodPrefix;
  }
  
  static std::string name ( ) { return "BM3D"; }
protected:
  void preprocess ( OptionsType &Options ) {
    if ( this->_input != NULL ) {
      CollaborativeNeighborhoodFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::preprocess ( Options );
      
      // For BM3D with Poisson Maximum-Likelihood ratios, the Anscombe transformed image is still required to fill the block stacks
      if ( Options.noiseType == NOISE_TYPE::POISSON && Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::MAXIMUMLIKELIHOOD ) {
        AnscombeForward<RealType> anscombeFW;
        this->applyVSTForwardAndScale ( *(this->_input), this->_preprocessedInput, anscombeFW, true );
      }
    } else throw aol::Exception ( "Input points to NULL", __FILE__, __LINE__ );
  }
  
  virtual void setBlockDistanceFunction ( const OptionsType &Options ) {
    if ( Options.method == BM3D_ESTIMATION_METHOD::HT ) {
      if ( Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::MAXIMUMLIKELIHOOD )
        this->_setBlockDistance = &BaseBM3DFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::setBlockDistancePoissonLikelihoodRatio;
      else
        this->_setBlockDistance = &BaseBM3DFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::setBlockDistanceL2NormSqr;
    } else if ( Options.method == BM3D_ESTIMATION_METHOD::WIENER ) {
      this->_setBlockDistance = &BaseBM3DFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::setBlockDistanceL2NormSqr;
    }
    
    setBlockDistanceFunctionNormFactor ( Options );
  }
  
  virtual void setBlockDistanceFunctionNormFactor ( const OptionsType &Options ) {
    if ( Options.method == BM3D_ESTIMATION_METHOD::HT ) {
      if ( Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::MAXIMUMLIKELIHOOD ) {
        this->_setBlockDistanceNormFactor = &BaseBM3DFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::setBlockDistanceNormFactorL2NormSqr;
        std::cerr << "NeighborhoodFilter: WARNING, L2NormSqr norm factor is used, despite NoiseType being POISSON!" << std::endl;
      } else
        this->_setBlockDistanceNormFactor = &BaseBM3DFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::setBlockDistanceNormFactorL2NormSqr;
    } else if ( Options.method == BM3D_ESTIMATION_METHOD::WIENER ) {
      this->_setBlockDistanceNormFactor = &BaseBM3DFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::setBlockDistanceNormFactorL2NormSqr;
    }
  }
    
  virtual void initializeOnce ( OptionsType &Options ) {
    CollaborativeNeighborhoodFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::initialize ( Options );
    this->_X0.setZero ( );
  }
  
  virtual void initializeIteration ( const OptionsType &Options ) {
    setBlockDistanceFunction ( Options );
    setParameters ( Options );
  }
  
  virtual SharedDataStructureType* allocateAndReturnSharedDataPtr ( const OptionsType &/*Options*/ ) const {
    if ( this->_input != NULL ) {
      return new SharedDataStructureType ( *(this->_input), this->_NX );
    } else throw aol::Exception ( "Input points to NULL", __FILE__, __LINE__ );
  }
  
  virtual PrivateDataStructureType* allocateAndReturnPrivateDataPtr ( const OptionsType &Options ) const {
    return new PrivateDataStructureType ( Options.method, this->_blockSize, this->_numNodes, _N2, Options.denoiseInTransformDomain, _transformXY, _transformZ, _htThreshold, this->_stdDev );
  }

  virtual void denoise ( const OptionsType &Options ) {
    initializeIteration ( Options );
    
    // Assemble an iterable container of all reference coordinates
    aol::RandomAccessContainer<qc::CoordType> refCoords;
    this->getReferenceCoordinates ( Options, refCoords );
    const int numRefCoords = refCoords.size ( );
    if ( Options.progressBar != NULL ) {
      if ( Options.method == BM3D_ESTIMATION_METHOD::HT ) Options.progressBar->setText ( "BM3D: Computing initial estimate (step 2/5)" );
      else if ( Options.method == BM3D_ESTIMATION_METHOD::WIENER ) Options.progressBar->setText ( "BM3D: Computing final estimate (step 5/5)" );
      Options.progressBar->start ( numRefCoords, 1 );
    }
    
    aol::DeleteFlagPointer<SharedDataStructureType> sharedDataPtr ( allocateAndReturnSharedDataPtr ( Options ), true );
#ifdef _OPENMP
    #pragma omp parallel
    {
#endif
      aol::DeleteFlagPointer<PrivateDataStructureType> privateDataPtr ( allocateAndReturnPrivateDataPtr ( Options ), true );
      fillBlocksFromPreviousEstimate ( Options, *privateDataPtr );
#ifdef _OPENMP
      #pragma omp for
#endif
      for ( int k=0; k<numRefCoords ; ++k ) {
        blockMatching ( Options, *privateDataPtr, refCoords[k] );
        blockStacking ( Options, *privateDataPtr );
        blockStackDenoising ( Options, *privateDataPtr );
        aggregate ( Options, *privateDataPtr, *sharedDataPtr );
#ifdef _OPENMP
        if ( Options.progressBar != NULL && privateDataPtr->thrIdx == 0 ) for ( int thr=0; thr<privateDataPtr->numThreads ; ++thr ) (*Options.progressBar)++;
#else
        if ( Options.progressBar != NULL ) (*Options.progressBar)++;
#endif
      }
      processDataStructures ( Options, *privateDataPtr, *sharedDataPtr );
#ifdef _OPENMP
      mergeSharedData ( *privateDataPtr, *sharedDataPtr );
    }
#endif
    if ( Options.progressBar != NULL ) Options.progressBar->finish ( );
    setEstimateFromBuffers ( Options, *sharedDataPtr );
    saveCurrentEstimate ( Options, *sharedDataPtr );
  }
public:
  virtual void blockMatching ( const OptionsType &Options,
                               const BlockCollectionType &Blocks,
                               const qc::CoordType &XRef, aol::RandomAccessContainer<qc::CoordType > &Sx, aol::Vector<RealType> &BlockDistances, int &Nz,
                               const PictureType &NumAggregates = PictureType ( ) ) const {
    // Calculate distances to all blocks in the search region
    int z = 0;
    RealType dist;
    for ( SimilaritySearchIterator<RealType, Dim>* searchIt = this->getSimilaritySearchIterator ( Options, XRef ); searchIt->notAtEnd ( ) ; ++(*searchIt) ) {
      this->setBlockDistance ( Blocks, dist, XRef, *(*searchIt) );
      if ( searchIt->isCurFinal ( ) && dist <= _tauMatch ) {
        BlockDistances[z] = dist;
        Sx[z] = *(*searchIt);
        ++z;
      }
      searchIt->update ( dist );
    }
    Nz = z;
    
    blockMatchingPostprocessing ( Options, Sx, BlockDistances, Nz, NumAggregates );
  }
  
  virtual void blockMatchingPostprocessing ( const OptionsType &Options,
                                             aol::RandomAccessContainer<qc::CoordType > &Sx, aol::Vector<RealType> &BlockDistances, int &Nz,
                                             const PictureType &/*NumAggregates*/ = PictureType ( ),
                                             const int /*NumNonLocalWindows*/ = 0 ) const {
    // If Nz > _N2 or Nz is not a power of 2 drop coordinates with highest distances until |Sx| = Nz = _N2 or |Sx| = Nz is a power of 2
    short finalNz;
    if ( Options.restrictStackSizeToPowersOfTwo ) finalNz = aol::Min<int> ( aol::Pow ( 2, floor ( log2 ( Nz ) ) ), _N2 );
    else finalNz = aol::Min<int> ( Nz, _N2 );
    if ( Nz > finalNz ){
      qc::CoordType tmpPos;
      RealType tmpDist;
      int zMin;
      for ( int z=0; z<finalNz ; ++z ) {
        zMin = z;
        for ( int z2=z+1; z2<Nz ; ++z2 )
          if ( BlockDistances[z2] < BlockDistances[zMin] ) zMin = z2;
        tmpPos = Sx[z];
        tmpDist = BlockDistances[z];
        Sx[z] = Sx[zMin];
        BlockDistances[z] = BlockDistances[zMin];
        Sx[zMin] = tmpPos;
        BlockDistances[zMin] = tmpDist;
      }
      Nz = finalNz;
    }
  }
  
  virtual void blockMatching ( const OptionsType &Options, PrivateDataStructureType &PrivateData, const qc::CoordType &XRef ) const {
    blockMatching ( Options, *PrivateData.blocksMethod, XRef, PrivateData.Sx, PrivateData.blockDistances, PrivateData.Nz );
  }
  
  virtual void blockMatchingPostprocessing ( const OptionsType &Options, PrivateDataStructureType &PrivateData ) const {
    blockMatchingPostprocessing ( Options, PrivateData.Sx, PrivateData.blockDistances, PrivateData.Nz );
  }
  
  virtual void setBlockDistanceL2NormSqrWithoutCenterPixel ( const BlockCollectionType &Blocks, RealType &Dist, const qc::CoordType &XRef, const qc::CoordType &X ) const {
    Dist = 0.0;
    for ( int k=0; k<(this->_blockNumPixels-1)/2; ++k )
      Dist += aol::Sqr<RealType> ( Blocks.get ( XRef, k ) - Blocks.get ( X, k ) );
    for ( int k=(this->_blockNumPixels-1)/2+1; k<this->_blockNumPixels; ++k )
      Dist += aol::Sqr<RealType> ( Blocks.get ( XRef, k ) - Blocks.get ( X, k ) );
  }
protected:
  virtual void fillBlocksFromPreviousEstimate ( const OptionsType &Options, BlockCollectionType &BlocksInitial, BlockCollectionType &BlocksEstimate ) const {
    if ( this->_input != NULL && this->_estimate != NULL ) {
      if ( Options.noiseType == NOISE_TYPE::POISSON && Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::MAXIMUMLIKELIHOOD && Options.method == BM3D_ESTIMATION_METHOD::HT )
        BlocksInitial.initialize ( *(this->_input), this->_blockSize, this->_blockAnchor );
      else {
        BlocksInitial.initialize ( this->_preprocessedInput, this->_blockSize, this->_blockAnchor );
        if ( Options.method == BM3D_ESTIMATION_METHOD::WIENER )
          BlocksEstimate.initialize ( *(this->_estimate), this->_blockSize, this->_blockAnchor );
      }
    } else throw aol::Exception ( "Input and/or estimate point to NULL", __FILE__, __LINE__ );
  }
  
  virtual void fillBlocksFromPreviousEstimate ( const OptionsType &Options, PrivateDataStructureType &PrivateData ) const {
    fillBlocksFromPreviousEstimate ( Options, PrivateData.blocksInitial, PrivateData.blocksEstimate );
  }
  
  virtual void blockStacking ( const OptionsType &/*Options*/,
                               const BlockCollectionType &Blocks,
                               const aol::RandomAccessContainer<qc::CoordType > &Sx, const short Nz,
                               BlockStackType &BlockStack ) const {
    for ( int z=0; z<Nz ; ++z ) Blocks.getBlock ( BlockStack[z], Sx[z] );
  }
  
  virtual void blockStacking ( const OptionsType &Options,
                               const BlockCollectionType &BlocksInitial,
                               const BlockCollectionType &BlocksEstimate,
                               const aol::RandomAccessContainer<qc::CoordType > &Sx, const short Nz,
                               BlockStackType &Initial, BlockStackType &Estimate ) const {
    blockStacking ( Options, BlocksInitial, Sx, Nz, Initial );
    if ( Options.method == BM3D_ESTIMATION_METHOD::WIENER )
      blockStacking ( Options, BlocksEstimate, Sx, Nz, Estimate );
  }
  
  virtual void blockStacking ( const OptionsType &Options, PrivateDataStructureType &PrivateData ) const {
    blockStacking ( Options, PrivateData.blocksInitial, PrivateData.blocksEstimate, PrivateData.Sx, PrivateData.Nz, PrivateData.initialBlockStack, PrivateData.estimateBlockStack );
  }
  
  virtual void blockStackDenoising ( const OptionsType &Options,
                                     HTBlockStackDenoisingOp<RealType, BlockType, BlockSizeType> &HTOp,
                                     WienerBlockStackDenoisingOp<RealType, BlockType, BlockSizeType> &WienerOp,
                                     const aol::Vector<RealType> &BlockDistances,
                                     const short Nz,
                                     const BlockStackType &Initial, const BlockStackType &Estimate,
                                     BlockStackType &Denoised, RealType &Weight ) const {
    if ( Options.denoiseInTransformDomain ) {
      if ( Options.method == BM3D_ESTIMATION_METHOD::HT )
        HTOp.apply ( Nz, Initial, Denoised, Weight );
      else
        WienerOp.apply ( Nz, Initial, Estimate, Denoised, Weight );
    } else {
      getWeightedMeanBlock ( Initial, BlockDistances, Nz, Denoised[0] );
      Weight = 1;
    }
  }
  
  virtual void blockStackDenoising ( const OptionsType &Options, PrivateDataStructureType &PrivateData ) const {
    blockStackDenoising ( Options, PrivateData.htOp, PrivateData.wienerOp, PrivateData.blockDistances, PrivateData.Nz,
                          PrivateData.initialBlockStack, PrivateData.estimateBlockStack, PrivateData.denoisedBlockStack, PrivateData.weight );
  }
  
  virtual void aggregate ( const OptionsType &Options, const aol::RandomAccessContainer<qc::CoordType > &Sx, const short Nz,
                           const BlockStackType &BlockStack, const RealType Weight,
                           PictureType &EBuff, PictureType &WBuff, qc::ScalarArray<unsigned int, Dim> &NumAggregates ) const = 0;
  
  virtual void aggregate ( const OptionsType &Options, PrivateDataStructureType &PrivateData, SharedDataStructureType &SharedData ) const {
#ifdef _OPENMP
    aggregate ( Options, PrivateData.Sx, PrivateData.Nz, PrivateData.denoisedBlockStack, PrivateData.weight,
                SharedData.eBuffsThreads[PrivateData.thrIdx], SharedData.wBuffsThreads[PrivateData.thrIdx], SharedData.numAggregatesThreads[PrivateData.thrIdx] );
#else
    aggregate ( Options, PrivateData.Sx, PrivateData.Nz, PrivateData.denoisedBlockStack, PrivateData.weight,
                SharedData.eBuff, SharedData.wBuff, SharedData.numAggregates );
#endif
  }
  
  virtual void processDataStructures ( const OptionsType &/*Options*/, const PrivateDataStructureType &/*PrivateData*/, const SharedDataStructureType &/*SharedData*/ ) { }
  
#ifdef _OPENMP
  virtual void mergeSharedData ( const PrivateDataStructureType &PrivateData, SharedDataStructureType &SharedData ) {
#pragma omp for
    for ( int k=0; k<SharedData.numBufferEntries ; ++k ) {
      for ( int thr=0; thr<PrivateData.numThreads ; ++thr ) {
        SharedData.eBuff[k] += SharedData.eBuffsThreads[thr][k];
        SharedData.wBuff[k] += SharedData.wBuffsThreads[thr][k];
      }
    }
#pragma omp for
    for ( int k=0; k<SharedData.numAggregatesEntries ; ++k )
      for ( int thr=0; thr<PrivateData.numThreads ; ++thr )
        SharedData.numAggregates[k] += SharedData.numAggregatesThreads[thr][k];
  }
#endif
  
  virtual void setEstimateFromBuffers ( const OptionsType &/*Options*/,
                                        const PictureType &EBuff, const PictureType &WBuff ) {
    this->_estimate->reallocate ( EBuff );
    for ( int k=0; k<this->_estimate->size ( ) ; ++k )
      (*this->_estimate)[k] = ( WBuff[k] != 0 ) ? EBuff[k] / WBuff[k] : 0;
  }
  
  virtual void setEstimateFromBuffers ( const OptionsType &Options, const SharedDataStructureType &SharedData ) {
    setEstimateFromBuffers ( Options, SharedData.eBuff, SharedData.wBuff );
  }
  
  virtual void saveCurrentEstimate ( const OptionsType &Options, const SharedDataStructureType &SharedData ) {
    if ( Options.outputDir != "" ) {
      this->_estimate->save ( aol::strprintf ( "%s/estimate_%s%s", Options.outputDir.c_str ( ), BM3D_ESTIMATION_METHOD::getIdentifier ( Options.method ).c_str ( ),
                                              qc::getDefaultArraySuffix ( PictureType::Dim ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
      SharedData.numAggregates.save ( aol::strprintf ( "%s/numAggregates_%s%s", Options.outputDir.c_str ( ), BM3D_ESTIMATION_METHOD::getIdentifier ( Options.method ).c_str ( ),
                                                       qc::getDefaultArraySuffix ( Dim ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    }
  }
  
  virtual void setParameters ( const OptionsType &Options ) {
    // Initially set all parameters for Normal Profile
    this->_searchWindowSize.setAll ( 39 );
    if ( Options.method == BM3D_ESTIMATION_METHOD::HT ) {
      _N2 = 16;
      _tauMatch = 3000;
    } else if ( Options.method == BM3D_ESTIMATION_METHOD::WIENER ) {
      _N2 = 32;
      _tauMatch = 400;
    }
    
    // Adjust parameters according to actually chosen profile
    if ( Options.profile == BM3D_PROFILE::LC ) {
      this->_searchWindowSize.setAll ( 25 );
      if ( BM3D_ESTIMATION_METHOD::WIENER ) _N2 = 16;
    }
    
    if ( Options.searchWindowSize.getNumX ( ) > 0 && Options.searchWindowSize.getNumY ( ) > 0 ) {
      if ( Options.searchWindowSize.getNumX ( ) % 2 == 1 && Options.searchWindowSize.getNumY ( ) % 2 == 1 ) this->_searchWindowSize = Options.searchWindowSize;
      else throw aol::Exception ( "Specified search window size must be odd!", __FILE__, __LINE__ );
    }
    
    if ( Options.method == BM3D_ESTIMATION_METHOD::HT && Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::MAXIMUMLIKELIHOOD ) {
      _tauMatch = -log(0.55); // TODO
      _tauMatch *= this->_blockNumPixels;
    } else
      _tauMatch *= this->_blockNumPixels / aol::Sqr<RealType> ( 255.0 );
    
    if ( Options.N2 > 0 ) _N2 = Options.N2;
    
    if ( Options.denoiseInTransformDomain ) {
      // Set parameters for block stack filtering
      _transformZ = "haar";
      if ( Options.method == BM3D_ESTIMATION_METHOD::HT ) {
        _transformXY = "bior1.5";
        _htThreshold = ( ( Options.profile == BM3D_PROFILE::HIGH ) ? 2.5 : 2.7 ) * this->_stdDev;
      } else if ( Options.method == BM3D_ESTIMATION_METHOD::WIENER ) {
        _transformXY = "dct";
      }
    } else {
      // Set parameters for NLM-like denoising
      _filterParameter = ( Options.filterParameter > 0 ) ? Options.filterParameter : 0.40 * this->_stdDev;
      _filterParameterSqr = aol::Sqr<RealType> ( _filterParameter );
      _stdDevSqr = aol::Sqr<RealType> ( this->_stdDev );
    }
  }
};

  
  
template <qc::Dimension Dim, typename BlockSizeType>
class BaseHasOverlapChecker : public aol::VectorFeasibilityChecker {
protected:
  qc::BitArray<Dim> _slots;
  BlockSizeType _blockSize;
public:
  BaseHasOverlapChecker ( ) : _blockSize ( static_cast<short> ( 0 ) ) { }
  
  BaseHasOverlapChecker ( const qc::CoordType NX, const BlockSizeType &BlockSize ) : _slots ( qc::GridSize<Dim> ( NX ) ), _blockSize ( BlockSize ) { }
  
  virtual bool isFeasible ( const aol::Vector<int> &X ) const {
    qc::CoordType x;
    for ( int d=0; d<Dim ; ++d ) x[d] = X[d];
    return !hasOverlap ( x );
  }
  
  virtual bool isFeasible ( const qc::CoordType &X ) const {
    return !hasOverlap ( X );
  }
  
  void resetOverlaps ( ) {
    _slots.setAll ( false );
  }
  
  virtual void addOverlap ( const qc::CoordType &X ) = 0;
  
  virtual bool hasOverlap ( const qc::CoordType &X ) const = 0;
  
  void setSlots ( const qc::BitArray<Dim> &Slots ) {
    _slots = Slots;
  }
  
  const qc::BitArray<Dim>& getSlotsConstReference ( ) const {
    return _slots;
  }
  
  void reallocate ( const qc::CoordType NX, const BlockSizeType &BlockSize ) {
    _slots.reallocate ( qc::GridSize<Dim> ( NX ) );
    _blockSize = BlockSize;
  }
  
  int numFeasible ( ) {
    int res = 0;
    qc::CoordType lower, upper;
    upper[0] = _slots.getNumX ( )-_blockSize.getNumX ( )+1;
    upper[1] = _slots.getNumY ( )-_blockSize.getNumY ( )+1;
    if ( Dim == qc::QC_3D ) upper[2] = _slots.getNumZ ( );
    for ( qc::RectangularIterator<Dim> it ( lower, upper ); it.notAtEnd ( ) ; ++it )
      if ( isFeasible ( *it ) ) ++res;
    return res;
  }
};

template <qc::Dimension Dim, typename BlockSizeType> class HasOverlapChecker;

template<qc::Dimension Dim>
class HasOverlapChecker<Dim, qc::GridSize<qc::QC_2D> > : public BaseHasOverlapChecker<Dim, qc::GridSize<qc::QC_2D> > {
public:
  HasOverlapChecker ( const qc::CoordType NX, const qc::GridSize<qc::QC_2D> &BlockSize )
  : BaseHasOverlapChecker<Dim, qc::GridSize<qc::QC_2D> > ( NX, BlockSize ) { }
  
  virtual void addOverlap ( const qc::CoordType &X ) {
    for ( int yBlock=0; yBlock<this->_blockSize.getNumY ( ); ++yBlock ) {
      const int xBlock = ( this->_blockSize.getNumX ( ) - 1 ) / 2;
      qc::CoordType x ( X );
      x[0] += xBlock;
      x[1] += yBlock;
      this->_slots.set ( x, true );
    }
  }
  
  virtual bool hasOverlap ( const qc::CoordType &X ) const {
    bool res = false;
    for ( int yBlock=0; yBlock<this->_blockSize.getNumY ( ); ++yBlock ) {
      const int xBlock = ( this->_blockSize.getNumX ( ) - 1 ) / 2;
      qc::CoordType x ( X );
      x[0] += xBlock;
      x[1] += yBlock;
      if ( this->_slots.get ( x ) ) res = true;
    }
    return res;
  }
};
  
  
  
  
/**
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, qc::Dimension Dim = qc::QC_2D, typename _PictureType = qc::ScalarArray<_RealType, qc::QC_2D>,
          typename _BlockType = qc::ScalarArray<_RealType, qc::QC_2D>, typename _BlockSizeType = qc::GridSize<qc::QC_2D>,
          typename FilterTrait = BM3DFilterTrait<_RealType, Dim, _PictureType, _BlockType, _BlockSizeType> >
class BM3DFilter;

  
/**
 * \brief Implements the BM3D filter for 2D scalar images (or stacks of 2D scalar images)
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, qc::Dimension Dim, typename FilterTrait>
class BM3DFilter<_RealType, Dim, qc::ScalarArray<_RealType, Dim>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D>, FilterTrait>
  : public BaseBM3DFilter<_RealType, Dim, qc::ScalarArray<_RealType, Dim>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D>, FilterTrait> {
public:
  typedef _RealType RealType;
  typedef qc::ScalarArray<_RealType, Dim> PictureType;
  typedef qc::ScalarArray<_RealType, qc::QC_2D> BlockType;
  typedef qc::GridSize<qc::QC_2D> BlockSizeType;
  
  typedef typename FilterTrait::OptionsType OptionsType;
  typedef BlockStack<RealType, BlockType, BlockSizeType> BlockStackType;
  typedef typename FilterTrait::BlockCollectionType BlockCollectionType;
    
  typedef qc::MultiArray<RealType, Dim, 3> ColoredPictureType;
protected:
  RealType _beta; // Parameter for 2D Kaiser window
  BlockType _kaiserWindow2D;
public:
  BM3DFilter ( ) : BaseBM3DFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait> ( ) { }
    
  void initializeOnce ( OptionsType &Options ) {
    BaseBM3DFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::initializeOnce ( Options );
    
    if ( Dim == qc::QC_3D ) {
      this->_NX[2] = this->_input->getNumZ ( );
      this->_XEnd[2] = this->_NX[2];
    }
  }
    
  void getLevelSetSamples ( const OptionsType &Options, const PictureType &Input, aol::VectorContainer<aol::Vector<RealType> > &Samples ) {
    aol::Vector<RealType> deltas;
    getLevelSetSamples ( Options, Input, Samples, deltas );
  }
  
  virtual void getLevelSetSamples ( const OptionsType &Options, const PictureType &Input, aol::VectorContainer<aol::Vector<RealType> > &Samples,
                                    aol::Vector<RealType> &/*Deltas*/ ) {
    OptionsType options ( Options );
    options.method = BM3D_ESTIMATION_METHOD::HT;
    options.N2 = options.levelSetNumPointsPerSample;
    options.similaritySearchMethod = BM3D_SIMILARITYSEARCH_METHOD::GLOBAL;
    options.restrictStackSizeToPowersOfTwo = false;
    options.noiseType = NOISE_TYPE::GAUSSIAN;
    options.sigma = aol::NumberTrait<RealType>::Inf;
    
    // 1. extract level-set coordinates from preprocessed input image via block-matching using L2-distance and static stack size
    aol::VectorContainer<aol::Vector<int> > samplesCoords;
    getLevelSetSampleCoords ( options, Input, samplesCoords );
    getLevelSetSamples ( Input, Samples, samplesCoords );
  }
  
  void getLevelSetSamples ( const PictureType &RawInput, aol::VectorContainer<aol::Vector<RealType> > &Samples, const aol::VectorContainer<aol::Vector<int> > &SampleCoords ) const {
    Samples.clear ( );
    for ( int i=0; i<SampleCoords.size ( ) ; ++i ) {
      aol::Vector<RealType> sample;
      for ( int j=0; j<SampleCoords[i].size ( ) ; ++j ) sample.pushBack ( RawInput[SampleCoords[i][j]] );
      Samples.pushBack ( sample );
    }
  }
  
  virtual void getLevelSetSampleCoords ( OptionsType &Options, const PictureType &Input, aol::VectorContainer<aol::Vector<int> > &SampleCoords ) {
    SampleCoords.clear ( );
    BlockCollectionType blocksInitial, blocksEstimate;
    this->prepareGetBlockStack ( Options, Input, blocksInitial, blocksEstimate );
    
    this->_tauMatch = aol::NumberTrait<RealType>::Inf;
    this->_setBlockDistance = static_cast<typename BaseNeighborhoodFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::BlockDistanceFunctionType> ( &BaseBM3DFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::setBlockDistanceL2NormSqrWithoutCenterPixel );
    
    aol::VectorContainer<aol::Vector<int> > sampleCoords;
    aol::Vector<RealType> sampleBlockDistancesSums;
    qc::CoordType xRef;
    BlockStackType blockStack ( Options.N2, this->_blockSize );
    aol::Vector<RealType> blockDistances ( this->_numNodes );
    aol::RandomAccessContainer<qc::CoordType > Sx ( this->_numNodes );
    int Nz = 0;
    aol::Vector<RealType> meanValues;
    HasOverlapChecker<Dim, BlockSizeType> hasOverlapChecker ( this->_NX, Options.blockSize );
    qc::FastILexMapper<Dim> mapper ( qc::GridSize<Dim> ( this->_NX ) );
    
    if ( Options.outputDir != "" ) aol::makeDirectory ( aol::strprintf ( "%s/noiseAnalysis", Options.outputDir.c_str ( ) ).c_str ( ), false );
    aol::ProgressBar<> progressBar ( "Analyzing noise", std::cerr );
    if ( !Options.quietMode ) progressBar.start ( Options.levelSetNumSamples * Options.levelSetOverSamplingFactor );
    int numReferences = 0;
    while ( numReferences < Options.levelSetNumSamples * Options.levelSetOverSamplingFactor ) {
      // 1. Choose free reference pixel that is nearest to the center of the largest interval between two sampling points
      if ( numReferences <= 1 ) {
        PictureType u ( Input, aol::STRUCT_COPY );
        u.setAll ( aol::NumberTrait<RealType>::Inf );
        if ( numReferences == 1 )
          u *= -1;
        qc::CoordType x0 ( this->_X0 );
        for ( int d=0; d<2 ; ++d )
          x0[d] = ( this->_blockSize[d] - 1 ) / 2;
        for ( qc::RectangularIterator<Dim> it ( x0, this->_XEnd ); it.notAtEnd ( ) ; ++it )
          u.set ( *it, Input.get ( *it ) );
        std::pair<int, RealType> indVal;
        if ( numReferences == 0 )
          indVal = u.getMinIndexAndValue ( );
        else
          indVal = u.getMaxIndexAndValue ( );
        for ( int d=0; d<Dim ; ++d ) xRef[d] = mapper.splitGlobalIndex ( indVal.first )[d];
        for ( int d=0; d<2 ; ++d )
          xRef[d] -= ( this->_blockSize[d] - 1 ) / 2;
      } else {
        meanValues.sortValues ( );
        RealType maxIntervalSize = 0;
        int maxIntervalIdx = 0;
        for ( int i=1; i<meanValues.size ( )-1 ; ++i ) {
          if ( meanValues[i+1] - meanValues[i] > maxIntervalSize ) {
            maxIntervalSize = meanValues[i+1] - meanValues[i];
            maxIntervalIdx = i;
          }
        }
        const RealType newMean = 0.5 * ( meanValues[maxIntervalIdx+1] + meanValues[maxIntervalIdx] );
        PictureType diffToNewMean ( Input, aol::STRUCT_COPY );
        diffToNewMean.setAll ( aol::NumberTrait<RealType>::Inf );
        for ( qc::RectangularIterator<Dim> it ( this->_X0, this->_XEnd ); it.notAtEnd ( ) ; ++it )
          if ( !hasOverlapChecker.hasOverlap ( *it ) ) {
            diffToNewMean.set ( *it, 0.0 );
            for ( int dx=-1; dx<=1 ; ++dx ) {
              qc::CoordType x ( *it );
              x[0] += ( this->_blockSize.getNumX ( ) - 1 ) / 2 + dx;
              diffToNewMean.add ( *it, Input.get ( x ) );
            }
            diffToNewMean.set ( *it, aol::Abs<RealType> ( diffToNewMean.get ( *it ) / 3.0 - newMean ) );
          }
        std::pair<int, RealType> minIndVal = diffToNewMean.getMinIndexAndValue ( );
        for ( int d=0; d<Dim ; ++d ) xRef[d] = mapper.splitGlobalIndex ( minIndVal.first )[d];
      }
      hasOverlapChecker.addOverlap ( xRef );
      
      // 2. Perform block-matching and add coordinates of the centers of the matched blocks as level-set samples
      this->getBlockStack ( Options, xRef, blocksInitial, blockStack, Sx, blockDistances, Nz );
      RealType blockDistancesSum = 0;
      for ( int z=0; z<Nz ; ++z ) blockDistancesSum += blockDistances[z];
      if ( blockDistancesSum < aol::NumberTrait<RealType>::Inf ) {
        const int kBlock = ( this->_blockNumPixels - 1 ) / 2;
        aol::Vector<int> sample ( Nz );
        RealType meanVal = 0;
        for ( int zStack=0; zStack<Nz ; ++zStack ) {
          qc::CoordType x ( Sx[zStack] );
          x[0] += ( this->_blockSize[0] - 1 ) / 2;
          x[1] += ( this->_blockSize[1] - 1 ) / 2;
          sample[zStack] = mapper.getGlobalIndex ( x );
          meanVal += blocksInitial.get ( Sx[zStack], kBlock );
        }
        
        sampleCoords.pushBack ( sample );
        sampleBlockDistancesSums.pushBack ( blockDistancesSum );
        meanValues.pushBack ( meanVal / static_cast<RealType> ( Nz ) );
        
        if ( Options.outputDir != "" && numReferences < 10 ) {
          qc::ScalarArray<RealType, qc::QC_2D> blockStackImg ( this->_blockSize.getNumX ( ), blockStack.getNumBlocks ( ) );
          for ( int xBlock=0; xBlock<this->_blockSize.getNumX ( ) ; ++xBlock )
            for ( int zStack=0; zStack<blockStack.getNumBlocks ( ) ; ++zStack )
              blockStackImg.set ( xBlock, zStack, blockStack[zStack].get ( xBlock, ( this->_blockSize[1] - 1 ) / 2 ) );
          blockStackImg.save ( aol::strprintf ( "%s/noiseAnalysis/blockStack_%d%s", Options.outputDir.c_str ( ), numReferences, qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
          blockDistances.saveASCII ( aol::strprintf ( "%s/noiseAnalysis/blockDistances_%d.csv", Options.outputDir.c_str ( ), numReferences ).c_str ( ), 12 );
        }
        
        if ( numReferences == 1 && Dim == qc::QC_2D ) {
          // Save the level-set with the largest mean value
          qc::MultiArray<RealType, qc::QC_2D, 3> inputLevelSet ( Input.getNumX ( ), Input.getNumY ( ) );
          for ( int d=0; d<3 ; ++d ) inputLevelSet[d] = Input;
          inputLevelSet.scaleValuesTo01 ( );
          for ( int i=0; i<sampleCoords[0].size ( ) ; ++i ) {
            inputLevelSet[0][sampleCoords[1][i]] = 1;
            inputLevelSet[1][sampleCoords[1][i]] = 0;
            inputLevelSet[2][sampleCoords[1][i]] = 0;
          }
          inputLevelSet.setOverflowHandlingToCurrentValueRange ( );
          inputLevelSet.savePNG ( aol::strprintf ( "%s/noiseAnalysis/levelSet.png", Options.outputDir.c_str ( ) ).c_str ( ) );
        }
        
        ++numReferences;
        if ( !Options.quietMode ) progressBar++;
      }
    }
    if ( !Options.quietMode ) progressBar.finish ( );
    
    if ( Options.outputDir != "" )
      hasOverlapChecker.getSlotsConstReference ( ).save ( aol::strprintf ( "%s/noiseAnalysis/slots.pgm", Options.outputDir.c_str ( ) ).c_str ( ) );
    
    // Keep only the levelSetNumSamples samples with the least sum of block distances
    if ( Options.levelSetOverSamplingFactor > 1 ) {
      aol::Vector<RealType> sampleBlockDistancesSumsSorted ( sampleBlockDistancesSums );
      sampleBlockDistancesSumsSorted.sortValues ( );
      const RealType threshold = sampleBlockDistancesSumsSorted[aol::Min<int> ( sampleBlockDistancesSumsSorted.size ( ) - 1, Options.levelSetNumSamples )];
      for ( int i=0; i<sampleBlockDistancesSums.size ( ) ; ++i ) {
        if ( sampleBlockDistancesSums[i] < threshold )
          SampleCoords.pushBack ( sampleCoords[i] );
      }
    } else {
      for ( int i=0; i<sampleCoords.size ( ) ; ++i )
        SampleCoords.pushBack ( sampleCoords[i] );
    }
  }
protected:
  void preprocess ( OptionsType &Options ) {
    BaseBM3DFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::preprocess ( Options );
    
    if ( Options.outputDir != "" ) this->_preprocessedInput.save ( aol::strprintf ( "%s/preprocessedInput%s", Options.outputDir.c_str ( ), qc::getDefaultArraySuffix ( Dim ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
  }
  
  using BaseBM3DFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::aggregate;
  void aggregate ( const OptionsType &Options, const aol::RandomAccessContainer<qc::CoordType > &Sx, const short Nz,
                   const BlockStackType &BlockStack, const RealType Weight,
                   PictureType &EBuff, PictureType &WBuff, qc::ScalarArray<unsigned int, Dim> &NumAggregates ) const {
    qc::FastILexMapper<Dim> mapper ( EBuff );
    qc::CoordType X;
    RealType weight;
    int flatIdx;
    const int nz = ( Options.denoiseInTransformDomain ) ? Nz : 1;
    for ( short yBlock=0; yBlock<this->_blockSize.getNumY ( ) ; ++yBlock ) {
      for ( short xBlock=0; xBlock<this->_blockSize.getNumX ( ) ; ++xBlock ) {
        weight = Weight * _kaiserWindow2D.get ( xBlock, yBlock );
        for ( short z=0; z<nz ; ++z ) {
          X = Sx[z];
          X[0] += xBlock;
          X[1] += yBlock;
          flatIdx = mapper.getGlobalIndex ( X );
          EBuff[flatIdx] += weight * BlockStack[z].get ( xBlock, yBlock );
          WBuff[flatIdx] += weight;
          ++NumAggregates[flatIdx];
        }
      }
    }
  }
    
  void setParameters ( const OptionsType &Options ) {
    this->_blockSize = Options.blockSize;
    this->_blockNumPixels = this->_blockSize.getNumberOfNodes ( );
    this->_refStep = 3;
    _beta = 2.0;
    
    if ( Options.profile == BM3D_PROFILE::LC ) {
      if ( Options.method == BM3D_ESTIMATION_METHOD::HT ) this->_refStep = 6;
      else if ( Options.method == BM3D_ESTIMATION_METHOD::WIENER ) this->_refStep = 5;
    }
    
    if ( Options.profile == BM3D_PROFILE::HIGH ) {
      this->_refStep = 2;
      if ( Options.method == BM3D_ESTIMATION_METHOD::HT )
        _beta = 2.5;
      else if ( Options.method == BM3D_ESTIMATION_METHOD::WIENER )
        _beta = 1.5;
    }
    
    if ( Options.refStep > 0 ) this->_refStep = Options.refStep;
    
    _kaiserWindow2D.reallocate ( this->_blockSize );
    if ( Options.denoiseInTransformDomain ) setKaiserWindow ( _kaiserWindow2D, _beta );
    else _kaiserWindow2D.setAll ( 1 );
    
    this->_XEnd[0] = this->_NX[0] - this->_blockSize.getNumX ( ) + 1;
    this->_XEnd[1] = this->_NX[1] - this->_blockSize.getNumY ( ) + 1;
    this->_numNodes = this->_XEnd[0] * this->_XEnd[1];
    if ( Dim == qc::QC_3D )
      this->_numNodes *= this->_NX[2];
    
    BaseBM3DFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::setParameters ( Options );
  }
                            
  int numAggregatesSum ( const PictureType &NumAggregates, const qc::CoordType &X ) const {
    int res = 0;
    qc::CoordType x ( X );
    for ( int yBlock=0; yBlock<this->_blockSize.getNumY ( ) ; ++yBlock ) {
      for ( int xBlock=0; xBlock<this->_blockSize.getNumX ( ) ; ++xBlock ) {
        x[0] = X[0] + xBlock;
        x[1] = X[1] + yBlock;
        res += NumAggregates.get ( x );
      }
    }
    return res;
  }
};
    

  
/**
 * \brief Implements the BM3D filter for multi-spectral 2D data (blocks are restricted to individual spectra)
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, typename FilterTrait>
class BM3DFilter<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_3D>, aol::Vector<_RealType>, int, FilterTrait>
  : public BaseBM3DFilter<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_3D>, aol::Vector<_RealType>, int, FilterTrait> {
public:
  typedef _RealType RealType;
  const static qc::Dimension Dim = qc::QC_2D;
  typedef qc::ScalarArray<_RealType, qc::QC_3D> PictureType;
  typedef aol::Vector<_RealType> BlockType;
  typedef int BlockSizeType;
   
  typedef typename FilterTrait::OptionsType OptionsType;  
  typedef BlockStack<RealType, BlockType, BlockSizeType> BlockStackType;  
public:
  BM3DFilter ( ) : BaseBM3DFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait> ( ) { }
protected:
  using BaseBM3DFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::aggregate;
  void aggregate ( const OptionsType &Options, const aol::RandomAccessContainer<qc::CoordType > &Sx, const short Nz,
                   const BlockStackType &BlockStack, const RealType Weight,
                   PictureType &EBuff, PictureType &WBuff, qc::ScalarArray<unsigned int, Dim> &NumAggregates ) const {
    qc::FastILexMapper<qc::QC_3D> mapper ( EBuff );
    int flatIdx;
    const int nz = ( Options.denoiseInTransformDomain ) ? Nz : 1;
    for ( short k=0; k<this->_blockSize ; ++k ) {
      for ( short z=0; z<nz ; ++z ) {
        flatIdx = mapper.getGlobalIndex ( k, Sx[z][0], Sx[z][1] );
        EBuff[flatIdx] += Weight * BlockStack[z][k];
        WBuff[flatIdx] += Weight;
        NumAggregates.add ( Sx[z][0], Sx[z][1], 1 );
      }
    }
  }

  void setParameters ( const OptionsType &Options ) {
    this->_blockSize = this->_input->getNumX ( );
    this->_blockNumPixels = this->_blockSize;
    this->_refStep = 1;
    
    this->_XEnd = this->_NX;
    this->_numNodes = this->_NX[0] * this->_NX[1];
    
    BaseBM3DFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::setParameters ( Options );
  }
                            
  int numAggregatesSum ( const PictureType &NumAggregates, const qc::CoordType &X ) const {
    return NumAggregates.get ( X );
  }
};

  
/**
 * \brief Implements the BM3D filter for multi-dimensional 2D data (e.g. color or multi-spectral data)
 * 
 *        In contrast to the special multi-spectral implementation, the block size in the 2D plane can be set to any power of 2 (like in the BM3D filter for scalar 2D images).
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, typename FilterTrait>
class BM3DFilter<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_3D>, qc::ScalarArray<_RealType, qc::QC_3D>, qc::GridSize<qc::QC_3D>, FilterTrait>
: public BaseBM3DFilter<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_3D>, qc::ScalarArray<_RealType, qc::QC_3D>, qc::GridSize<qc::QC_3D>, FilterTrait> {
public:
  typedef _RealType RealType;
  const static qc::Dimension Dim = qc::QC_2D;
  typedef qc::ScalarArray<_RealType, qc::QC_3D> PictureType;
  typedef qc::ScalarArray<RealType, qc::QC_3D> BlockType;
  typedef qc::GridSize<qc::QC_3D> BlockSizeType;
  
  typedef typename FilterTrait::OptionsType OptionsType;
  typedef BlockStack<RealType, BlockType, BlockSizeType> BlockStackType;
public:
  BM3DFilter ( ) : BaseBM3DFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait> ( ) { }
protected:
  using BaseBM3DFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::aggregate;
  void aggregate ( const OptionsType &Options, const aol::RandomAccessContainer<qc::CoordType > &Sx, const short Nz,
                   const BlockStackType &BlockStack, const RealType Weight,
                   PictureType &EBuff, PictureType &WBuff, qc::ScalarArray<unsigned int, Dim> &NumAggregates ) const {
    qc::FastILexMapper<qc::QC_3D> mapper ( EBuff );
    int flatIdx;
    const int nz = ( Options.denoiseInTransformDomain ) ? Nz : 1;
    for ( int yBlock=0; yBlock<this->_blockSize.getNumZ ( ) ; ++yBlock ) {
      for ( int xBlock=0; xBlock<this->_blockSize.getNumY ( ) ; ++xBlock ) {
        for ( int x=0; x<this->_blockSize.getNumX ( ) ; ++x ) {
          for ( short z=0; z<nz ; ++z ) {
            flatIdx = mapper.getGlobalIndex ( x, Sx[z][0] + xBlock, Sx[z][1] + yBlock );
            EBuff[flatIdx] += Weight * BlockStack[z].get ( x, xBlock, yBlock );
            WBuff[flatIdx] += Weight;
            NumAggregates.add ( Sx[z][0] + xBlock, Sx[z][1] + yBlock, 1 );
          }
        }
      }
    }
  }
  
  void setParameters ( const OptionsType &Options ) {
    if ( this->_input != NULL ) {
      this->_blockSize[0] = this->_input->getNumX ( );
      this->_blockSize[1] = Options.blockSize.getNumY ( );
      this->_blockSize[2] = Options.blockSize.getNumZ ( );
      this->_blockNumPixels = this->_blockSize.getNumberOfNodes ( );
      this->_refStep = 3;
      
      if ( Options.refStep > 0 ) this->_refStep = Options.refStep;
      
      this->_XEnd[0] = this->_NX[0] - this->_blockSize.getNumY ( ) + 1;
      this->_XEnd[1] = this->_NX[1] - this->_blockSize.getNumZ ( ) + 1;
      this->_numNodes = this->_XEnd[0] * this->_XEnd[1];
      
      BaseBM3DFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::setParameters ( Options );
    } else throw aol::Exception ( "Input points to NULL", __FILE__, __LINE__ );
  }
  
  int numAggregatesSum ( const PictureType &NumAggregates, const qc::CoordType &X ) const {
    int res = 0;
    for ( int yBlock=0; yBlock<this->_blockSize.getNumZ ( ) ; ++yBlock )
      for ( int xBlock=0; xBlock<this->_blockSize.getNumY ( ) ; ++xBlock )
        res += NumAggregates.get ( X[0] + xBlock, X[1] + yBlock, X[2] );
    return res;
  }
};

} // namespace im

#endif