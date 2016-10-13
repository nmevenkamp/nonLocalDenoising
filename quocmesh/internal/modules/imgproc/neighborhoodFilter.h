#ifndef __NEIGHBORHOODFILTER_H
#define __NEIGHBORHOODFILTER_H

//// standard
//#define _USE_MATH_DEFINES
//#define _OPEN_SYS
//#include <climits>
//#include <dirent.h>
//// Under VC++, dirent.h defines macros conflicting with std code.
//#ifdef _MSC_VER
//#ifdef max
//#undef max
//#endif
//#ifdef min
//#undef min
//#endif
//#endif

// quocmesh
#include <aol.h>
#include <ctrlCCatcher.h>
#include <multiArray.h>
#include <op.h>
#include <scalarArray.h>
#include <qmException.h>
#include <bzipiostream.h>
#include "blockStructures.h"
#include "statistics.h"
#include "anscombe.h"
#include "similaritySearchIterators.h"
#include "referenceBlockIterators.h"
#include <noiseAnalysis.h>


namespace im {


template <typename FilterType, typename OptionsType>
class BM3DFilterInterface;


template <typename T>
int resolveIdentifier ( const std::string &Identifier ) {
  int identifierIdx = -1;
  for ( int i=0; i<T::NUM ; ++i )
    if ( Identifier == T::getIdentifier ( i ) ) identifierIdx = i;
  if ( identifierIdx == -1 ) throw aol::Exception ( aol::strprintf ( "Invalid identifier '%s'!", Identifier.c_str() ).c_str(), __FILE__, __LINE__ );
  return identifierIdx;
}


class NOISE_TYPE {
public:
  static const int GAUSSIAN = 0;  // Remove additive Gaussian white noise (AGWN)
  static const int POISSON  = 1;  // Remove Poisson noise
  static const int CMOS     = 2;  // Remove mixed scaled Poisson and additive Gaussian noise (proposed in the literature for digital images from CMOS/CCD sensors)
  
  static const int NUM = 3;
  
  static std::string getIdentifier ( const int NoiseType ) {
    if ( NoiseType == GAUSSIAN ) return "G";
    else if ( NoiseType == POISSON ) return "P";
    else if ( NoiseType == CMOS ) return "CMOS";
    else throw aol::Exception ( "Did not recognize noise type!", __FILE__, __LINE__ );
  }
  
  static std::string toString ( const int NoiseType ) {
    if ( NoiseType == GAUSSIAN ) return "Gaussian noise";
    else if ( NoiseType == POISSON ) return "Poisson noise";
    else if ( NoiseType == CMOS ) return "CMOS (mixed scaled Poisson + Gaussian)";
    else throw aol::Exception ( "Did not recognize noise type!", __FILE__, __LINE__ );
  }
};

class POISSON_NOISE_ADAPTATION {
public:
  static const int ANSCOMBE           = 0;  // Three step procedure: 1. Anscombe transform, 2. Remove Gaussian noise, 3. Inverse Anscombe transform
  static const int MAXIMUMLIKELIHOOD  = 1;  // Directly apply denoising procedures on Poisson noise statistics using maximum-likelihood ratios as intensity distance measure
  
  static const int NUM = 2;
  
  static std::string getIdentifier ( const int PoissonNoiseAdaptation ) {
    if ( PoissonNoiseAdaptation == ANSCOMBE ) return "A";
    else if ( PoissonNoiseAdaptation == MAXIMUMLIKELIHOOD ) return "ML";
    else throw aol::Exception ( "Did not recognize Poisson noise adaptation!", __FILE__, __LINE__ );
  }
  
  static std::string toString ( const int PoissonNoiseAdaptation ) {
    if ( PoissonNoiseAdaptation == ANSCOMBE ) return "Anscombe transform";
    else if ( PoissonNoiseAdaptation == MAXIMUMLIKELIHOOD ) return "Maximum-likelihood ratios";
    else throw aol::Exception ( "Did not recognize Poisson noise adaptation!", __FILE__, __LINE__ );
  }
};

class SIMILARITYSEARCH_METHOD {
public:
  static const int LOCAL  = 0;  // A local square is used as a search window
  static const int GLOBAL = 1;  // Similarity search is performed on the whole image
  
  static const int NUM = 2;
  
  static std::string getIdentifier ( const int SimilaritySearchMethod ) {
    if ( SimilaritySearchMethod == LOCAL ) return "l";
    else if ( SimilaritySearchMethod == GLOBAL ) return "gl";
    else throw aol::Exception ( "NH similarity search method: Did not recognize similarity search method!", __FILE__, __LINE__ );
  }
  
  static std::string toString ( const int SimilaritySearchMethod ) {
    if ( SimilaritySearchMethod == LOCAL ) return "Local";
    else if ( SimilaritySearchMethod == GLOBAL ) return "Global";
    else throw aol::Exception ( "NH similarity search method: Did not recognize similarity search method!", __FILE__, __LINE__ );
  }
};

struct REFERENCEITERATION_METHOD {
  static const int GLOBAL = 0;  // The whole image is denoised (possibly except some boundary pixels)
  static const int PIXELS = 1;  // Only a specified set of pixels is denoised
  
  static const int NUM = 2;
};


template<typename _RealType, typename _PictureType, typename _BlockType, typename _BlockSizeType>
struct NHFilterOptions {
  _PictureType *groundTruth;
  std::string outputDir;
  bool quietMode;
  aol::ProgressBar<> *progressBar;
  _BlockSizeType blockSize;
  qc::GridSize<qc::QC_2D> searchWindowSize;
  int noiseType;
  _RealType alpha, sigma, mu;
  int similaritySearchMethod;
  _RealType maxVal, smoothingCoefficient;
  int levelSetNumSamples;
  int poissonNoiseAdaptation;
  bool removeHotPixels;
  
  
  NHFilterOptions ( const std::string OutputDir = "", const bool Quiet = true )
    : groundTruth ( NULL ),
      outputDir ( OutputDir ), quietMode ( Quiet ), progressBar ( NULL ),
      blockSize ( static_cast<short> ( 0 ) ), searchWindowSize ( static_cast<short> ( 0 ) ), noiseType ( NOISE_TYPE::GAUSSIAN ),
      alpha ( 0.0 ), sigma ( 0.0 ), mu ( 0.0 ),
      similaritySearchMethod ( SIMILARITYSEARCH_METHOD::LOCAL ),
      maxVal ( 255.0 ), smoothingCoefficient ( 1.0 ),
      levelSetNumSamples ( 100 ),
      poissonNoiseAdaptation ( POISSON_NOISE_ADAPTATION::ANSCOMBE ),
      removeHotPixels ( false ) { }
  
  NHFilterOptions ( const NHFilterOptions<_RealType, _PictureType, _BlockType, _BlockSizeType> &Options )
    : groundTruth ( Options.groundTruth ),
      outputDir ( Options.outputDir ), quietMode ( Options.quietMode ), progressBar ( Options.progressBar ),
      blockSize ( Options.blockSize ), searchWindowSize ( Options.searchWindowSize ), noiseType ( Options.noiseType ),
      alpha ( Options.alpha ), sigma ( Options.sigma ), mu ( Options.mu ),
      similaritySearchMethod ( Options.similaritySearchMethod ),
      maxVal ( Options.maxVal ), smoothingCoefficient ( Options.smoothingCoefficient ),
      levelSetNumSamples ( Options.levelSetNumSamples ),
      poissonNoiseAdaptation ( Options.poissonNoiseAdaptation ),
      removeHotPixels ( Options.removeHotPixels ) { }
};


template <typename _RealType, qc::Dimension Dim, typename _PictureType, typename _BlockType, typename _BlockSizeType>
class NHFilterTrait {
public:
  typedef NHFilterOptions<_RealType, _PictureType, _BlockType, _BlockSizeType> OptionsType;
  typedef BlockCollection<_RealType, Dim, _PictureType, _BlockType, _BlockSizeType> BlockCollectionType;
  typedef SIMILARITYSEARCH_METHOD SimilaritySearchMethodType;
};


template <typename _BlockSizeType> struct getBlockSizeStr;

template<>
struct getBlockSizeStr<qc::GridSize<qc::QC_2D> > {
  static std::string doGetBlockSizeStr ( const qc::GridSize<qc::QC_2D> &BlockSize ) {
    if ( BlockSize.getNumX ( ) == BlockSize.getNumY ( ) )
      return aol::strprintf ( "%d", BlockSize.getNumX ( ) );
    else
      return aol::strprintf ( "(%d,%d)", BlockSize.getNumX ( ), BlockSize.getNumY ( ) );
  }
};

template<>
struct getBlockSizeStr<int> {
  static std::string doGetBlockSizeStr ( const int &/*BlockSize*/ ) {
    return "";
  }
};

template<>
struct getBlockSizeStr<qc::GridSize<qc::QC_3D> > {
  static std::string doGetBlockSizeStr ( const qc::GridSize<qc::QC_3D> &BlockSize ) {
    if ( BlockSize.getNumY ( ) == BlockSize.getNumZ ( ) )
      return aol::strprintf ( "%d", BlockSize.getNumY ( ) );
    else
      return aol::strprintf ( "(%d,%d)", BlockSize.getNumY ( ), BlockSize.getNumZ ( ) );
  }
};

  
/**
 * \brief Abstract implementation of a neighborhood filter (provides basic functionality such as noise transformation, managing of block distance functions, etc.)
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template<typename _RealType, qc::Dimension Dim, typename _PictureType, typename _BlockType, typename _BlockSizeType, typename FilterTrait>
class BaseNeighborhoodFilter {
public:
  typedef _RealType RealType;
  typedef _PictureType PictureType;
  typedef _BlockType BlockType;
  typedef _BlockSizeType BlockSizeType;
  typedef aol::FullMatrix<RealType> MatrixType;
  typedef qc::MultiArray<RealType, Dim, 3> ColoredPictureType;
  typedef typename FilterTrait::OptionsType OptionsType;
  typedef typename FilterTrait::BlockCollectionType BlockCollectionType;
  typedef typename FilterTrait::SimilaritySearchMethodType SimilaritySearchMethodType;
  typedef void ( BaseNeighborhoodFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::*BlockDistanceFunctionType )( const BlockCollectionType&, RealType&, const qc::CoordType&, const qc::CoordType& ) const;
  typedef void ( BaseNeighborhoodFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::*BlockDistanceNormFactorFunctionType )( const BlockCollectionType&, RealType&, const qc::CoordType& ) const;
protected:
  // Input
  const PictureType *_input;                                              // Noisy input image
  PictureType _preprocessedInput;                                         // Input image after pre-processing
  PictureType *_estimate;                                                 // Current estimate of the mean pixel intensities
  RealType _smoothingCoefficient;                                         // Amplification factor of noise standard deviation that can be used to adjust smoothing after application of VST
  RealType _stdDev;                                                       // Noise standard deviation in AGWN model (e.g. after Ancombe transformation in case of (mixed) Poisson(-Gaussian) noise)
  // Partial terminating
  bool _catchCtrlC;
  mutable sigfunc _previousCtrlCHandler;
  
  // Similarity search
  aol::Vec2<short> _blockAnchor;
  BlockSizeType _blockSize;                                               // Block Size (e.g. aol::Vec2<short> in case of 2D scalar blocks)
  int _blockNumPixels;                                                    // Block total linear size (including any extra dimensions))
  qc::CoordType _X0, _NX, _XEnd;                                          // Parameters defining the dimensions of the subregion of the image containing well-defined block anchors
  int _numNodes;                                                          // Parameters defining the dimensions of the subregion of the image containing well-defined block anchors
  qc::GridSize<qc::QC_2D> _searchWindowSize;                              // Size of the window used for similarity search (in case of a local similarity search)
  aol::MultiVector<short> _pixels;                                        // Array of reference pixels for which similarity search is performed (in case of a sparse reconstruction)
  
  // Blocks
  BlockCollectionType _blocksInitial;                                     // Block collection containing values from input image
  BlockDistanceFunctionType _setBlockDistance;                            // Pointer to currently used block distance function
  BlockDistanceNormFactorFunctionType _setBlockDistanceNormFactor;        // Pointer to currently used block distancen normalization factor function
  
  // Variables used for (generalized) Anscombe transformation
  RealType _alpha, _sigma, _mu;                                           // mixed Poisson-Gaussian noise parameters
  RealType _GATb;
  const RealType _scaleRange, _scaleShift;
  RealType _minTransformed, _maxTransformed;
  aol::VectorContainer<aol::Vector<RealType> > _samples;
public:
  BaseNeighborhoodFilter ( )
    : _catchCtrlC ( false ),
      _blockSize ( static_cast<short> ( 0 ) ),
      _searchWindowSize ( 0, 0 ),
      _scaleRange ( 0.7 ), _scaleShift ( 0.5 * ( 1.0 - _scaleRange ) ) { }
  
  void setCatchCtrlC ( bool catchCtrlC ) {
    _catchCtrlC = catchCtrlC;
  }
  
  virtual ~BaseNeighborhoodFilter ( ) { }
  
  virtual void apply ( const OptionsType &Options, const PictureType &Arg, PictureType &Dest ) = 0;
  
  virtual void setSimilaritySearchGrid ( const OptionsType &Options, const PictureType &Arg, const qc::CoordType &XRef,
                                         aol::RandomAccessContainer<qc::CoordType> &CoordsFinal, aol::RandomAccessContainer<qc::CoordType> &CoordsTMP ) = 0;
  
  static void setSimilaritySearchGridImage ( const PictureType &Input, const qc::CoordType &XRef,
                                             const aol::RandomAccessContainer<qc::CoordType> &CoordsFinal, const aol::RandomAccessContainer<qc::CoordType> &CoordsTmp,
                                             ColoredPictureType &Dest ) {
    Dest.reallocate ( Input );
    Dest.setZero ( );
    Dest[0] = Input; Dest[1] = Input; Dest[2] = Input;
    Dest.scaleValuesTo01 ( );
    for ( int i=0; i<CoordsTmp.size ( ) ; ++i ) Dest[2].set ( CoordsTmp[i], aol::Min<RealType> ( 1, Dest[2].get ( CoordsTmp[i] ) + 0.5 ) );
    for ( int i=0; i<CoordsFinal.size ( ) ; ++i ) {
      Dest[0].set ( CoordsFinal[i], 1 );
      Dest[1].set ( CoordsFinal[i], 0 );
      Dest[2].set ( CoordsFinal[i], 0 );
    }
    Dest[0].set ( XRef, 0 );
    Dest[1].set ( XRef, 1 );
    Dest[2].set ( XRef, 0 );
    Dest.setOverflowHandlingToCurrentValueRange ( );
  }
  
  virtual void getLevelSetSamples ( const OptionsType &/*Options*/, const PictureType &/*Input*/, aol::VectorContainer<aol::Vector<RealType> > &/*Samples*/ ) {
    throw aol::UnimplementedCodeException ( "Not implemented", __FILE__, __LINE__ );
  }
  
  static std::string getMethodIdentifier ( const OptionsType &Options ) {
    return getMethodPrefix ( Options ) + "-" + name ( );
  }
  
  static std::string getMethodPrefix ( const OptionsType &Options ) {
    std::stringstream ss;
    ss << NOISE_TYPE::getIdentifier ( Options.noiseType );
    ss << "-" << Options.smoothingCoefficient;
    if ( Options.noiseType == NOISE_TYPE::POISSON ) ss << "-" << POISSON_NOISE_ADAPTATION::getIdentifier ( Options.poissonNoiseAdaptation );
    ss << "-" << SimilaritySearchMethodType::getIdentifier ( Options.similaritySearchMethod );
    ss << "-" << getBlockSizeStr<BlockSizeType>::doGetBlockSizeStr ( Options.blockSize );
    return ss.str ( );
  }
  
  static const std::string name ( ) { return "NHF"; }
  
  int getNumNodes ( ) {
    return _numNodes ( );
  }
  
  BlockSizeType blockSize ( ) const {
    return _blockSize;
  }
  
  RealType alpha ( ) {
    return _alpha;
  }
  
  RealType sigma ( ) {
    return _sigma;
  }
  
  RealType mu ( ) {
    return _mu;
  }
protected:
  virtual void initialize ( OptionsType &Options ) {
    if ( _input != NULL ) {
      if ( !Options.quietMode && Options.progressBar == NULL ) Options.progressBar = new aol::ProgressBar<>;
      _preprocessedInput.reallocate ( *_input );
      
      // The following has to be done before everything else, since it applies similarity search with altered settings to get level sets
      if ( ( Options.noiseType == NOISE_TYPE::GAUSSIAN && Options.sigma <= 0 )
        || ( Options.noiseType == NOISE_TYPE::CMOS && ( Options.alpha <= 0 || Options.sigma < 0 ) ) ) {
        _samples.clear ( );
        getLevelSetSamples ( Options, *this->_input, _samples );
      }
      
      setBlockDistanceFunction ( Options );
    } else throw aol::Exception ( "Input points to NULL", __FILE__, __LINE__ );
  }
  
  virtual void preprocess ( OptionsType &Options ) {
    if ( _input != NULL ) {
      if ( Options.noiseType == NOISE_TYPE::GAUSSIAN ) {
        if ( Options.sigma > 0 )
          _sigma = Options.sigma;
        else {
          const std::string outputDirNoiseAnalyzer = ( Options.outputDir != "" ) ? aol::strprintf ( "%s/noiseAnalysis", Options.outputDir.c_str ( ) ) : "";
          if ( outputDirNoiseAnalyzer != "" ) aol::makeDirectory ( outputDirNoiseAnalyzer.c_str ( ), false );
          _sigma = 0;
          im::NoiseAnalyzer<RealType>::getGaussianVarianceParams ( _samples, _sigma, Options.quietMode, outputDirNoiseAnalyzer );
        }
        _preprocessedInput = *_input;
        _preprocessedInput /= static_cast<RealType> ( Options.maxVal );
        _stdDev = _sigma / static_cast<RealType> ( Options.maxVal );
      } else if ( Options.noiseType == NOISE_TYPE::POISSON ) {
        if ( Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::ANSCOMBE ) {
          AnscombeForward<RealType> anscombeFW;
          applyVSTForwardAndScale ( *_input, _preprocessedInput, anscombeFW, true );
        } else if ( Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::MAXIMUMLIKELIHOOD ) {
          _preprocessedInput = *_input;
        } else throw aol::Exception ( "Did not recognize Poisson noise adaptation!", __FILE__, __LINE__ );
      } else if ( Options.noiseType == NOISE_TYPE::CMOS ) {
        if ( Options.alpha > 0 && Options.sigma >= 0 ) {
          _alpha = Options.alpha;
          _sigma = Options.sigma;
          _mu = Options.mu;

          _GATb = aol::Sqr<RealType> ( _sigma ) - _alpha * _mu;
        } else {
          const std::string outputDirNoiseAnalyzer = ( Options.outputDir != "" ) ? aol::strprintf ( "%s/noiseAnalysis", Options.outputDir.c_str ( ) ) : "";
          if ( outputDirNoiseAnalyzer != "" ) aol::makeDirectory ( outputDirNoiseAnalyzer.c_str ( ), false );
          im::NoiseAnalyzer<RealType>::getCMOSVarianceParams ( _samples, _alpha, _GATb, Options.quietMode, outputDirNoiseAnalyzer );
          im::NoiseAnalyzer<RealType>::getCMOSPedestalParam ( _samples, _mu, _alpha, _GATb, Options.quietMode, outputDirNoiseAnalyzer );
          im::NoiseAnalyzer<RealType>::getCMOSGaussianStdDevParam ( _sigma, _alpha, _GATb, _mu, Options.quietMode, outputDirNoiseAnalyzer );
        }
        GeneralizedAnscombeOp<RealType>::clampParameters ( _alpha, _GATb, _input->getMinValue ( ), Options.quietMode );
        GeneralizedAnscombeForward<RealType> gatFW ( _alpha, _GATb );
        applyVSTForwardAndScale ( *_input, _preprocessedInput, gatFW, true );
      } else throw aol::Exception ( "Did not recognize noise model!", __FILE__, __LINE__ );
      
      if ( Options.noiseType == NOISE_TYPE::POISSON || Options.noiseType == NOISE_TYPE::CMOS ) {
        _smoothingCoefficient = Options.smoothingCoefficient;
        _stdDev *= _smoothingCoefficient;
      }
      
      if ( Options.removeHotPixels ) {
        if ( Options.noiseType == NOISE_TYPE::POISSON && POISSON_NOISE_ADAPTATION::MAXIMUMLIKELIHOOD )
          std::cerr << "NeighborhoodFilter: removing hot pixels not supported for Poisson noise without VST! Skipping this step." << std::endl;
        
        removeHotPixels ( Options );
      }
    } else throw aol::Exception ( "Input points to NULL", __FILE__, __LINE__ );
  }
  
  virtual void postprocess ( const OptionsType &Options ) {
    if ( _estimate != NULL ) {
      if ( Options.outputDir != "" && ( Options.noiseType == NOISE_TYPE::POISSON || Options.noiseType == NOISE_TYPE::CMOS ) )
        saveMethodNoise ( _preprocessedInput, *_estimate, aol::strprintf ( "%s/methodNoise_VST%s", Options.outputDir.c_str ( ), getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ) );
      
      if ( Options.noiseType == NOISE_TYPE::GAUSSIAN ) {
        (*_estimate) *= Options.maxVal; // Scale estimate back to [0,255]
      } else if ( Options.noiseType == NOISE_TYPE::POISSON ) {
        if ( Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::ANSCOMBE ) {
          AnscombeInverse<RealType> anscombeInverse;
          scaleBackAndApplyVSTInverse ( *_estimate, anscombeInverse );
        }
      } else if ( Options.noiseType == NOISE_TYPE::CMOS ) {
        GeneralizedAnscombeInverse<RealType> gatInverse ( _alpha, _GATb );
        scaleBackAndApplyVSTInverse ( *_estimate, gatInverse );
      } else throw aol::Exception ( "Did not recognize noise model!", __FILE__, __LINE__ );
    } else throw aol::Exception ( "Estimate points to NULL", __FILE__, __LINE__ );
  }
  
  virtual void applyVSTForwardAndScale ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest,
                                         const VSTOp<RealType> &VSTForward,
                                         bool ResetParams = false ) {
    VSTForward.apply ( Arg, Dest );
    
    // Scale image to [0.15,0.85]
    if ( ResetParams ) {
      _minTransformed = VSTForward.getMinTransformed ( );
      _maxTransformed = Dest.getMaxValue ( );
    }
    Dest.addToAll ( -_minTransformed );
    Dest /= ( _maxTransformed - _minTransformed );
    Dest *= _scaleRange;
    Dest.addToAll ( _scaleShift );
    
    // Set noise standard deviation
    if ( ResetParams ) _stdDev = 1.0 / ( _maxTransformed - _minTransformed ) * _scaleRange;
  }
  
  virtual void scaleBackAndApplyVSTInverse ( aol::Vector<RealType> &ArgDest,
                                             const aol::Op<aol::Vector<RealType> > &VSTInverse ) const {
    ArgDest.addToAll ( - _scaleShift );
    ArgDest /= _scaleRange;
    ArgDest *= ( _maxTransformed - _minTransformed );
    ArgDest.addToAll ( _minTransformed );
    
    // Apply inverse variance-stabilizing transformation
    VSTInverse.apply ( ArgDest, ArgDest );
  }
  
  virtual ReferenceBlockIterator<RealType, Dim>* getReferenceIterator ( const OptionsType &Options ) const {
    return new GlobalReferenceBlockIterator<RealType, Dim> ( _X0, _XEnd, _numNodes, Options.progressBar, Options.quietMode );
  }
  
  void getReferenceCoordinates ( const OptionsType &Options, aol::RandomAccessContainer<qc::CoordType> &RefCoords ) const {
    RefCoords.clear ( );
    OptionsType options ( Options );
    options.quietMode = true;
    for ( ReferenceBlockIterator<RealType, Dim> *refIt = this->getReferenceIterator ( options ); refIt->notAtEnd ( ) ; ++(*refIt) )
      RefCoords.pushBack ( *(*refIt) );
  }
  
  virtual SimilaritySearchIterator<RealType, Dim>* getSimilaritySearchIterator ( const OptionsType &Options, const qc::CoordType &XRef ) const {
    if ( Options.similaritySearchMethod == SIMILARITYSEARCH_METHOD::GLOBAL )
      return new GlobalSimilaritySearchIterator<RealType, Dim> ( _X0, _XEnd );
    else if ( Options.similaritySearchMethod == SIMILARITYSEARCH_METHOD::LOCAL )
      return new LocalSimilaritySearchIterator<RealType, Dim> ( XRef, _X0, _XEnd, _searchWindowSize );
    else
      throw aol::Exception ( "Neighborhood filter: Did not recognize similarity search method!", __FILE__, __LINE__ );
  }
  
  virtual void setBlockDistanceFunction ( const OptionsType &Options ) {
    if ( Options.noiseType == NOISE_TYPE::GAUSSIAN
      || ( Options.noiseType == NOISE_TYPE::POISSON && Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::ANSCOMBE )
      || Options.noiseType == NOISE_TYPE::CMOS )
      _setBlockDistance = &BaseNeighborhoodFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::setBlockDistanceL2NormSqr;
    else if ( Options.noiseType == NOISE_TYPE::POISSON )
      _setBlockDistance = &BaseNeighborhoodFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::setBlockDistancePoissonLikelihoodRatio;
    else
      throw aol::Exception ( "Could not determine appropriate block distance measure!", __FILE__, __LINE__ );
    
    setBlockDistanceNormFactorFunction ( Options );
  }
  
  virtual void setBlockDistanceNormFactorFunction ( const OptionsType &Options ) {
    if ( Options.noiseType == NOISE_TYPE::GAUSSIAN
        || ( Options.noiseType == NOISE_TYPE::POISSON && Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::ANSCOMBE )
        || Options.noiseType == NOISE_TYPE::CMOS )
      _setBlockDistanceNormFactor = &BaseNeighborhoodFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::setBlockDistanceNormFactorL2NormSqr;
    else if ( Options.noiseType == NOISE_TYPE::POISSON ) {
      _setBlockDistanceNormFactor = &BaseNeighborhoodFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::setBlockDistanceNormFactorL2NormSqr;
      std::cerr << "NeighborhoodFilter: WARNING, L2NormSqr norm factor is used, despite NoiseType being POISSON!" << std::endl;
    } else
      throw aol::Exception ( "Could not determine appropriate block distance measure!", __FILE__, __LINE__ );
  }
  
  virtual void setBlockDistance ( const BlockCollectionType &Blocks, RealType &Dist, const qc::CoordType &XRef, const qc::CoordType &X ) const {
    (this->*this->_setBlockDistance) ( Blocks, Dist, XRef, X );
  }
  
  virtual void setBlockDistanceNormFactor ( const BlockCollectionType &Blocks, RealType &NormFactor, const qc::CoordType &XRef ) const {
    (this->*this->_setBlockDistanceNormFactor) ( Blocks, NormFactor, XRef );
  }
  
  virtual void setBlockDistanceL2NormSqr ( const BlockCollectionType &Blocks, RealType &Dist, const qc::CoordType &XRef, const qc::CoordType &X ) const {
    Dist = 0.0;
    for ( int k=0; k<_blockNumPixels; ++k )
      Dist += aol::Sqr<RealType> ( Blocks.get ( XRef, k ) - Blocks.get ( X, k ) );
  }
  
  virtual void setBlockDistanceNormFactorL2NormSqr ( const BlockCollectionType &Blocks, RealType &NormFactor, const qc::CoordType &XRef ) const {
    NormFactor = 0.0;
    for ( int k=0; k<_blockNumPixels; ++k )
      NormFactor += aol::Sqr<RealType> ( Blocks.get ( XRef, k ) );
  }
  
  virtual void setBlockDistancePoissonLikelihoodRatio ( const BlockCollectionType &Blocks, RealType &Dist, const qc::CoordType &XRef, const qc::CoordType &X ) const {
    Dist = 0.0;
    for ( int k=0; k<_blockNumPixels; ++k ) {
      // TODO: precompute log(k) for any k in [0, 2*max] if 2*max is not too large (less than the image size preferably, otherwise, one would assume Gaussian model anyways)
      const int k1 = Blocks.get ( XRef, k ), k2 = Blocks.get ( X, k );
      Dist += ( ( k1 > 0 ) ? k1 * log ( static_cast<RealType> ( k1 ) ) : 0 ) + ( ( k2 > 0 ) ? k2 * log ( static_cast<RealType> ( k2 ) ) : 0 ) - ( ( k1 + k2 > 0 ) ? ( k1 + k2 ) * log ( 0.5 * ( k1 + k2 ) ) : 0 );
    }
  }
  
  virtual void removeHotPixels ( OptionsType &Options ) {
    aol::Vector<int> positions;
    do {
      detectHotPixels ( Options, _preprocessedInput, positions );
      inpaint ( _preprocessedInput, positions );
    } while ( positions.size ( ) > 0 );
  }
  
  virtual void detectHotPixels ( const OptionsType &/*Options*/, const PictureType &/*Picture*/, aol::Vector<int> &/*Positions*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  virtual void inpaint ( PictureType &/*Input*/, const aol::Vector<int> &/*Positions*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void setCtrlCHandler () const {
    if (_catchCtrlC)
      _previousCtrlCHandler = signal ( InterruptSignal, aol::ctrlCHandler );
  }
  
  void unsetCtrlCHandler () const {
    if (_catchCtrlC)
      signal ( InterruptSignal, _previousCtrlCHandler );
  }
  
  bool wantsInterrupt() const {
    if (!_catchCtrlC || !aol::getCtrlCState())
      return false;
    else
      return true;
  }
  
  template <typename FilterType, typename FilterOptionsType>
  friend class BM3DFilterInterface;
};

/**
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, qc::Dimension Dim, typename _PictureType, typename _BlockType, typename _BlockSizeType, typename FilterTrait> class NeighborhoodFilter;

/**
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, typename FilterTrait>
class NeighborhoodFilter<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_2D>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D>, FilterTrait>
  : public BaseNeighborhoodFilter<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_2D>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D>, FilterTrait> {
public:
  typedef _RealType RealType;
  const static qc::Dimension Dim = qc::QC_2D;
  typedef qc::ScalarArray<_RealType, qc::QC_2D> PictureType;
  typedef qc::ScalarArray<_RealType, qc::QC_2D> BlockType;
  typedef qc::GridSize<qc::QC_2D> BlockSizeType;
  
  typedef typename FilterTrait::OptionsType OptionsType;
    
  typedef qc::MultiArray<RealType, Dim, 3> ColoredPictureType;
  NeighborhoodFilter ( ) : BaseNeighborhoodFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait> ( ) { }
protected:
  void initialize ( OptionsType &Options ) {
    this->_NX[0] = this->_input->getNumX ( );
    this->_NX[1] = this->_input->getNumY ( );
    this->_NX[2] = 1;
    
    BaseNeighborhoodFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::initialize ( Options );
  }
    
  void inpaint ( PictureType &Input, const aol::Vector<int> &Positions ) const {
    // Use local spiraling to find a set of healthy values (not needing to be inpainted)
    // Then define inpainted value as median of those values
    const int nVals = 8, nMax = 4 * aol::Max<int> ( Input.getNumX ( ), Input.getNumY ( ) );
    qc::FastILexMapper<qc::QC_2D> mapper ( this->_NX[0], this->_NX[1] );
    int x, y;
    for ( int i=0; i<Positions.size ( ) ; ++i ) {
      mapper.splitGlobalIndex ( Positions[i], x, y );
      aol::Vector<RealType> vals;
      int n = 0;
      for ( qc::LocalSpiralIterator<qc::QC_2D> it ( aol::Vec2<short> ( x, y ) ); vals.size ( ) < nVals ; ++it, ++n ) {
        if ( (*it)[0] >= 0 && (*it)[0] < this->_NX[0] && (*it)[1] >= 0 && (*it)[1] < this->_NX[1] && Positions.numOccurence ( mapper.getGlobalIndex ( (*it)[0], (*it)[1] ) ) == 0 )
          vals.pushBack ( Input.get ( *it ) );
        if ( n > nMax ) break;
      }
      if ( vals.size ( ) > 0 ) Input[Positions[i]] = vals.getMedianValue ( );
      else std::cerr << "NeighborhoodFilter: failed at inpainting pixel (" << x << "," << y << ")" << std::endl;
    }
  }
};

/**
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, typename FilterTrait>
class NeighborhoodFilter<_RealType, qc::QC_3D, qc::ScalarArray<_RealType, qc::QC_3D>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D>, FilterTrait>
: public BaseNeighborhoodFilter<_RealType, qc::QC_3D, qc::ScalarArray<_RealType, qc::QC_3D>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D>, FilterTrait> {
public:
  typedef _RealType RealType;
  const static qc::Dimension Dim = qc::QC_3D;
  typedef qc::ScalarArray<_RealType, qc::QC_3D> PictureType;
  typedef qc::ScalarArray<_RealType, qc::QC_2D> BlockType;
  typedef qc::GridSize<qc::QC_2D> BlockSizeType;
  
  typedef typename FilterTrait::OptionsType OptionsType;
  
  typedef qc::MultiArray<RealType, Dim, 3> ColoredPictureType;
  NeighborhoodFilter ( ) : BaseNeighborhoodFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait> ( ) { }
protected:
  void initialize ( OptionsType &Options ) {
    if ( this->_input != NULL ) {
      this->_NX[0] = this->_input->getNumX ( );
      this->_NX[1] = this->_input->getNumY ( );
      this->_NX[2] = this->_input->getNumZ ( );
      
      BaseNeighborhoodFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::initialize ( Options );
    } else throw aol::Exception ( "Input points to NULL", __FILE__, __LINE__ );
  }
  
  void inpaint ( PictureType &Input, const aol::Vector<int> &Positions ) const {
    // Use local spiraling to find a set of healthy values (not needing to be inpainted)
    // Then define inpainted value as median of those values
    const int nVals = 8, nMax = 4 * aol::Max<int> ( Input.getNumX ( ), Input.getNumY ( ) );
    qc::FastILexMapper<Dim> mapper ( qc::GridSize<Dim> ( this->_NX ) );
    int x, y, z;
    for ( int i=0; i<Positions.size ( ) ; ++i ) {
      mapper.splitGlobalIndex ( Positions[i], x, y, z );
      aol::Vector<RealType> vals;
      int n = 0;
      for ( qc::LocalSpiralIterator<qc::QC_2D> it ( aol::Vec2<short> ( x, y ) ); vals.size ( ) < nVals ; ++it, ++n ) {
        if ( (*it)[0] >= 0 && (*it)[0] < this->_NX[0] && (*it)[1] >= 0 && (*it)[1] < this->_NX[1] && Positions.numOccurence ( mapper.getGlobalIndex ( (*it)[0], (*it)[1], z ) ) == 0 )
          vals.pushBack ( Input.get ( (*it)[0], (*it)[1], z ) );
        if ( n > nMax ) break;
      }
      if ( vals.size ( ) > 0 ) Input[Positions[i]] = vals.getMedianValue ( );
      else std::cerr << "NeighborhoodFilter: failed at inpainting pixel (" << x << "," << y << ")" << std::endl;
    }
  }
};

/**
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, typename BlockType, typename BlockSizeType, typename FilterTrait>
class NeighborhoodFilter<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_3D>, BlockType, BlockSizeType, FilterTrait>
  : public BaseNeighborhoodFilter<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_3D>, BlockType, BlockSizeType, FilterTrait> {
public:
  typedef _RealType RealType;
  const static qc::Dimension Dim = qc::QC_2D;
  typedef qc::ScalarArray<_RealType, qc::QC_3D> PictureType;
  
  typedef typename FilterTrait::OptionsType OptionsType;
  NeighborhoodFilter ( ) : BaseNeighborhoodFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait> ( ) { }
protected:
  void initialize ( OptionsType &Options ) {
    this->_NX[0] = this->_input->getNumY ( );
    this->_NX[1] = this->_input->getNumZ ( );
    this->_NX[2] = 1;
    
    BaseNeighborhoodFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::initialize ( Options );
  }
};


/**
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template<typename _RealType, qc::Dimension Dim, typename _PictureType, typename _BlockType, typename _BlockSizeType, typename FilterTrait>
class CollaborativeNeighborhoodFilter : public NeighborhoodFilter<_RealType, Dim, _PictureType, _BlockType, _BlockSizeType, FilterTrait> {
  typedef _RealType RealType;
  typedef typename FilterTrait::OptionsType OptionsType;
protected:
  int _refStep;
  
  virtual ReferenceBlockIterator<RealType, Dim>* getReferenceIterator ( const OptionsType &Options ) const {
    return new GlobalReferenceBlockIterator<RealType, Dim> ( this->_X0, this->_XEnd, this->_numNodes, Options.progressBar, Options.quietMode, _refStep );
  }
};
  
} // namespace im

#endif