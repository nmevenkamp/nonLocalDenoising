#ifndef __NONLOCALMEANSFILTER_H
#define __NONLOCALMEANSFILTER_H


#include "neighborhoodFilter.h"
#include "blockStructures.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace im {

class NLM_SIMILARITYSEARCH_METHOD : public SIMILARITYSEARCH_METHOD {
public:
  static std::string getIdentifier ( const int SimilaritySearchMethod ) {
    return SIMILARITYSEARCH_METHOD::getIdentifier ( SimilaritySearchMethod );
  }
  
  static std::string toString ( const int SimilaritySearchMethod ) {
    return SIMILARITYSEARCH_METHOD::toString ( SimilaritySearchMethod );
  }
};


template <typename _RealType, typename _PictureType = qc::ScalarArray<_RealType, qc::QC_2D>, typename _BlockType = qc::ScalarArray<_RealType, qc::QC_2D>, typename _BlockSizeType = qc::GridSize<qc::QC_2D> >
struct NLMOptions : NHFilterOptions<_RealType, _PictureType, _BlockType, _BlockSizeType> {
  _RealType filterParameter;
  aol::MultiVector<short> pixels;
  
  NLMOptions ( const std::string OutputDir = "", const bool Quiet = true )
    : NHFilterOptions<_RealType, _PictureType, _BlockType, _BlockSizeType> ( OutputDir, Quiet ), filterParameter ( 0 ) { }
  
  NLMOptions ( const NLMOptions<_RealType, _PictureType, _BlockType, _BlockSizeType> &Options )
    : NHFilterOptions<_RealType, _PictureType, _BlockType, _BlockSizeType> ( Options ), filterParameter ( Options.filterParameter ) { }
};


template <typename _RealType, qc::Dimension Dim, typename _PictureType, typename _BlockType, typename _BlockSizeType>
class NonLocalMeansFilterTrait {
public:
  typedef NLMOptions<_RealType, _PictureType> OptionsType;
  typedef BlockCollection<_RealType, Dim, _PictureType, _BlockType, _BlockSizeType> BlockCollectionType;
  typedef NLM_SIMILARITYSEARCH_METHOD SimilaritySearchMethodType;
};


//  The following class NonLocalMeansFilter has been implemented by
//
//    Niklas Mevenkamp, email: mevenkamp@aices.rwth-aachen.de
//
//  based on the description of the Non-local Means filter as published in:
//
//      A. Buades, B. Coll, J.M. Morel "A review of image denoising methods, with a new one"
//      Multiscale Modeling and Simulation, Vol. 4 (2), pp: 490-530, 2006. DOI: 10.1137/040616024
//      available at http://epubs.siam.org/doi/abs/10.1137/040616024
//
//  Online applet and source code of patchwise implemenation: http://www.ipol.im/pub/art/2011/bcm_nlm/
//
//  AUTHORS:
//      Antoni Buades toni.buades@uib.es, CNRS-Paris Descartes
//      Bartomeu Coll tomeu.coll@uib.es, Universitat Illes Balears
//      Jean-Michel Morel morel@cmla.ens-cachan.fr, CMLA, ENS-Cachan
/**
 * \brief Abstract implementation of the Non-local means filter (needs specialization depending on dimension of input image, node space, blocks)
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, qc::Dimension Dim, typename _PictureType, typename _BlockType, typename _BlockSizeType, typename FilterTrait>
class BaseNonLocalMeansFilter : public NeighborhoodFilter<_RealType, Dim, _PictureType, _BlockType, _BlockSizeType, FilterTrait> {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
  typedef _BlockType BlockType;
  typedef _BlockSizeType BlockSizeType;
  
  typedef aol::FullMatrix<RealType> MatrixType;
  typedef qc::MultiArray<RealType, qc::QC_2D, 3> ColoredPictureType;
  
  typedef typename FilterTrait::OptionsType OptionsType;
  typedef typename FilterTrait::BlockCollectionType BlockCollectionType;
protected:
  RealType _filterParameter, _filterParameterSqr, _stdDevSqr;
public:  
  void apply ( const OptionsType &Options, const PictureType &Arg, PictureType &Dest ) {
    this->_input = &Arg;
    this->_estimate = &Dest;

    OptionsType options ( Options );
    this->initialize ( options );
    this->preprocess ( options );
    setParameters ( Options );
    initializeBlocks ( Options );
    
    denoise ( options );
    
    this->postprocess ( options );
    
    // (Optional) console output
    if ( !Options.quietMode && options.groundTruth != NULL && (*options.groundTruth).getNumX ( ) == this->_estimate->getNumX ( ) && (*options.groundTruth).getNumY ( ) == this->_estimate->getNumY ( ) ) {
      if ( options.noiseType == NOISE_TYPE::GAUSSIAN ) std::cerr << "PSNR = " << aol::PSNR<RealType> ( *options.groundTruth, *(this->_estimate), 255 ) << " dB" << std::endl;
      else if ( options.noiseType == NOISE_TYPE::POISSON ) std::cerr << "PSNR = " << aol::PSNR<RealType> ( *options.groundTruth, *(this->_estimate) ) << " dB" << std::endl;
    }
  }
  
  void setSimilaritySearchGrid ( const OptionsType &/*Options*/, const PictureType &/*Arg*/, const qc::CoordType &/*XRef*/,
                                 aol::RandomAccessContainer<qc::CoordType> &/*CoordsFinal*/, aol::RandomAccessContainer<qc::CoordType> &/*CoordsTmp*/ ) {
    throw aol::UnimplementedCodeException ( "Not implemented", __FILE__, __LINE__ );
  }
  
  virtual void setWeights ( const OptionsType &Options, const PictureType &Arg, const qc::CoordType &XRef, PictureType &Weights, const bool ReInitialize = true ) = 0;
  
  void cropBorders ( PictureType &Picture ) {
    Picture.crop ( aol::Vec2<int> ( this->_blockAnchor[0], this->_blockAnchor[1] ),
                   aol::Vec2<int> ( Picture.getNumX ( ) - 2 * this->_blockAnchor[0], Picture.getNumY ( ) - 2 * this->_blockAnchor[1] ) );
  }
  
  static std::string getMethodIdentifier ( const OptionsType &Options ) {
    return getMethodPrefix ( Options ) + "-" + name ( );
  }
  
  static std::string getMethodPrefix ( const OptionsType &Options ) {
    std::stringstream ss;
    ss << NeighborhoodFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::getMethodPrefix ( Options );
    ss << "-" << Options.filterParameter;
    return ss.str ( );
  }
  
  static std::string name ( ) { return "NLM"; }
protected:
  virtual void denoise ( const OptionsType &/*Options*/ ) = 0;
  
  void initializeBlocks ( const OptionsType &/*Options*/ ) {
    this->_blocksInitial.initialize ( this->_preprocessedInput, this->_blockSize, this->_blockAnchor );
  }
  
  virtual void setParameters ( const OptionsType &Options ) = 0;
  
  void setOperators ( const OptionsType &Options );
  
  inline RealType getWeight ( const RealType Dist ) {
    return exp ( -aol::Max<RealType> ( Dist / this->_blockNumPixels - 2 * _stdDevSqr, 0.0 ) / _filterParameterSqr  );
  }
};

  
/**
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, qc::Dimension Dim = qc::QC_2D, typename _PictureType = qc::ScalarArray<_RealType, qc::QC_2D>,
          typename _BlockType = qc::ScalarArray<_RealType, qc::QC_2D>, typename _BlockSizeType = qc::GridSize<qc::QC_2D>,
          typename FilterTrait = NonLocalMeansFilterTrait<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_2D>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D> > >
class NonLocalMeansFilter;

/**
 * \brief Implements the Non-local means filter for 2D scalar images
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, qc::Dimension Dim, typename FilterTrait>
class NonLocalMeansFilter<_RealType, Dim, qc::ScalarArray<_RealType, Dim>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D>, FilterTrait>
  : public BaseNonLocalMeansFilter<_RealType, Dim, qc::ScalarArray<_RealType, Dim>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D>, FilterTrait> {
  typedef _RealType RealType;
  typedef qc::ScalarArray<_RealType, Dim> PictureType;
  typedef qc::ScalarArray<_RealType, qc::QC_2D> BlockType;
  typedef qc::GridSize<qc::QC_2D> BlockSizeType;
    
  typedef typename FilterTrait::OptionsType OptionsType;
public:
  void setWeights ( const OptionsType &Options, const PictureType &Arg, const qc::CoordType &XRef, PictureType &Weights, const bool ReInitialize = true ) {
    if ( ReInitialize ) {
      this->_input = &Arg;

      OptionsType options ( Options );
      this->initialize ( options );
      this->preprocess ( options );
      setParameters ( Options );
      this->initializeBlocks ( Options );
    }
    
    if ( this->_input != NULL ) {
      if ( !Options.quietMode ) std::cerr << "NonLocalMeansFilter: Calculating weights for pixel " << XRef << ".." << std::endl;
      
      RealType dist;
      
      for ( SimilaritySearchIterator<RealType, Dim>* searchIt = this->getSimilaritySearchIterator ( Options, XRef ); searchIt->notAtEnd ( ) ; ++(*searchIt) ) {
        this->setBlockDistance ( this->_blocksInitial, dist, XRef, *(*searchIt) );
        Weights.set ( *(*searchIt), this->getWeight ( dist ) );
        searchIt->update ( dist );
      }
    } else throw aol::Exception ( "Input points to NULL", __FILE__, __LINE__ );
  }
protected:
  void denoise ( const OptionsType &Options ) {
    if ( this->_input != NULL && this->_estimate != NULL ) {
      // Assemble an iterable container of all reference coordinates
      aol::RandomAccessContainer<qc::CoordType> refCoords;
      this->getReferenceCoordinates ( Options, refCoords );
      const int numRefCoords = refCoords.size ( );
      if ( Options.progressBar != NULL ) {
        Options.progressBar->setText ( "NLM: Denoising" );
        Options.progressBar->start ( numRefCoords );
      }
  #ifdef _OPENMP
      #pragma omp parallel
      {
        const int numThreads = omp_get_num_threads ( );
        const int thrIdx = omp_get_thread_num ( );
  #endif
        RealType dist = 0, weight = 0, normFactor = 0;
  #ifdef _OPENMP
        #pragma omp for
  #endif
        for ( int k=0; k<numRefCoords ; ++k ) {
          this->_estimate->set ( refCoords[k], 0 );
          normFactor = 0;
        
          for ( SimilaritySearchIterator<RealType, Dim>* searchIt = this->getSimilaritySearchIterator ( Options, refCoords[k] ); searchIt->notAtEnd ( ) ; ++(*searchIt) ) {
            this->setBlockDistance ( this->_blocksInitial, dist, refCoords[k], *(*searchIt) );
            weight = this->getWeight ( dist );
            normFactor += weight;
            this->_estimate->set ( refCoords[k], this->_estimate->get ( refCoords[k] ) + weight * this->_preprocessedInput.get ( *(*searchIt) ) );
          }
        
          this->_estimate->set ( refCoords[k], this->_estimate->get ( refCoords[k] ) / normFactor );
  #ifdef _OPENMP
          if ( Options.progressBar != NULL && thrIdx == 0 ) for ( int thr=0; thr<numThreads ; ++thr) (*Options.progressBar)++;
  #else
          if ( Options.progressBar != NULL ) (*Options.progressBar)++;
  #endif
        }
  #ifdef _OPENMP
      }
  #endif
      if ( Options.progressBar != NULL ) Options.progressBar->finish ( );
    } else throw aol::Exception ( "Input and/or estimate point to NULL", __FILE__, __LINE__ );
  }
    
  void setParameters ( const OptionsType &Options ) {
    if ( Options.blockSize.getNumX ( ) && Options.blockSize.getNumY ( ) > 0 ) {
      if ( Options.blockSize.getNumX ( ) % 2 == 1 && Options.blockSize.getNumY ( ) % 2 == 1 ) this->_blockSize = Options.blockSize;
      else throw aol::Exception ( "Specified block size must be odd!", __FILE__, __LINE__ );
    } else if ( Options.noiseType == NOISE_TYPE::GAUSSIAN && this->_stdDev > 0 ) {
      if ( Options.sigma <= 15 ) this->_blockSize.setAll ( 3 );
      else if ( Options.sigma <= 30 ) this->_blockSize.setAll ( 5 );
      else if ( Options.sigma <= 45 ) this->_blockSize.setAll ( 7 );
      else if ( Options.sigma <= 75 ) this->_blockSize.setAll ( 9 );
      else this->_blockSize.setAll ( 11 );
    } else this->_blockSize.setAll ( 7 );
    this->_blockAnchor.set ( ( this->_blockSize.getNumX ( ) - 1 ) / 2, ( this->_blockSize.getNumY ( ) - 1 ) / 2 );
    this->_blockNumPixels = this->_blockSize.getNumX ( ) * this->_blockSize.getNumY ( );
    this->_X0[0] = this->_blockAnchor[0];
    this->_X0[1] = this->_blockAnchor[1];
    this->_XEnd[0] = this->_X0[0] + this->_NX[0] - this->_blockSize.getNumX ( ) + 1;
    this->_XEnd[1] = this->_X0[1] + this->_NX[1] - this->_blockSize.getNumY ( ) + 1;
    this->_numNodes = ( this->_NX[0] - this->_blockSize.getNumX ( ) + 1 ) * ( this->_NX[1] - this->_blockSize.getNumY ( ) + 1 ) * this->_input->getNumZ ( );
    
    if ( Options.similaritySearchMethod == NLM_SIMILARITYSEARCH_METHOD::LOCAL ) {
      if ( Options.searchWindowSize.getNumX ( ) > 0 && Options.searchWindowSize.getNumY ( ) > 0 ) {
        if ( Options.searchWindowSize.getNumX ( ) % 2 == 1 && Options.searchWindowSize.getNumY ( ) % 2 == 1 ) this->_searchWindowSize = Options.searchWindowSize;
        else throw aol::Exception ( "Specified search window size must be odd!", __FILE__, __LINE__ );
      } else if ( this->_stdDev > 0 && Options.noiseType == NOISE_TYPE::GAUSSIAN ) {
        if ( Options.sigma <= 30 ) this->_searchWindowSize.setAll ( 21 );
        else this->_searchWindowSize.setAll ( 35 );
      } else this->_searchWindowSize.setAll ( 35 );
    }
    
    if ( Options.filterParameter > 0 ) this->_filterParameter = Options.filterParameter;
    else if ( Options.sigma > 0 || this->_stdDev > 0 ) {
      this->_filterParameter = this->_stdDev;
      if ( Options.noiseType == NOISE_TYPE::GAUSSIAN ) {
        if ( Options.sigma <= 30 ) this->_filterParameter *= 0.4;
        else if ( Options.sigma <= 75 ) this->_filterParameter *= 0.35;
        else if ( Options.sigma <= 100 ) this->_filterParameter *= 0.3;
        else throw aol::Exception ( "Removing Gaussian noise with standard deviation exceeding 100 not feasible!", __FILE__, __LINE__ );
      } else if ( Options.noiseType == NOISE_TYPE::POISSON && Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::ANSCOMBE ) {
        this->_filterParameter *= 0.35;
      } else if ( Options.noiseType == NOISE_TYPE::CMOS ) {
        this->_filterParameter *= 0.35;
      } else throw aol::Exception ( "Automatic filter parameter settings are only provided for Gaussian noise (or Anscombe transformed Poisson noise)!", __FILE__, __LINE__ );
    } else throw aol::Exception ( "Neither filter parameter nor noise standard deviation was specified!", __FILE__, __LINE__ );
    this->_filterParameterSqr = aol::Sqr<RealType> ( this->_filterParameter );
    this->_stdDevSqr = aol::Sqr<RealType> ( this->_stdDev );
  }
    
//  RealType getPUniformityOptimalFilterParameter ( const OptionsType &Options ) {
//    if ( Options.noiseType == NOISE_TYPE::POISSON )
//      throw aol::Exception ( "P-values are not properly defined for Poisson noise statistics!", __FILE__, __LINE__ );
//    
//    if ( !this->_quietMode )
//      std::cerr << "NonLocalMeansFilter: searching for p-uniformity optimal filter parameter.." << std::endl;
//    
//    const int numBins = 100;
//    const RealType p = 0.5, alpha = 0.05;
//    RealType lBound = 0, rBound = aol::Sqr<RealType> ( this->_input.getMaxValue ( ) ), h0 = 0.5 * ( lBound + rBound ), h1 = h0, chiSqr0 = numBins * 2, chiSqr1 = chiSqr0;
//    int numPixels = 128, numItInner = 0, numItOuter = 1, pValuesQuantileThreshold;
//    aol::RandomGenerator randomGenerator;
//    randomGenerator.randomize ( );
//    
//    if ( !this->_quietMode )
//      std::cerr << "Iteration " << numItOuter << "." << numItInner << ": h0=" << h0 << " in [" << lBound << ", " << rBound << "]" << std::endl;
//    
//    do {
//      chiSqr0 = chiSqr1;
//      aol::MultiVector<int> pixels ( numPixels, 2 );
//      aol::Vector<int> min ( 2 ), max ( 2 );
//      min[0] = this->_blockAnchor[0]; min[1] = this->_blockAnchor[1];
//      max[0] = this->_NX - this->_blockAnchor[0]; max[1] = this->_NY - this->_blockAnchor[1];
//      randomGenerator.rIntMultiVecPairwiseDifferent ( pixels, min, max );
//      pValuesQuantileThreshold = DistributionTester<RealType>::getPValuesQuantileThreshold ( numPixels, p, alpha );
//      PictureType means ( this->_NX, this->_NY );
//      OptionsType options ( Options );
//      options.pixels.resize ( pixels.numComponents ( ), pixels[0].size ( ) );
//      for ( int c=0; c<pixels.numComponents ( ) ; ++c )
//        for ( int i=0; i<pixels[c].size ( ) ; ++i )
//          options.pixels[c][i] = pixels[c][i];
//      do {
//        h0 = h1;
//        options.filterParameter = h1;
//        this->apply ( options, means );
//        
//        GaussianNoiseImageAnalyzer<RealType, PictureType> noiseImageAnalyzer ( this->_input, means, numBins, this->_outputDir );
//        const std::pair<int, int> pValuesInSymmetricQuantiles ( noiseImageAnalyzer.getNumPValuesInSymmetricQuantiles ( p, options.pixels ) );
//        chiSqr1 = noiseImageAnalyzer.getChiSquareOfUniformDistributionOfPValues ( options.pixels );
//        
//        if ( pValuesInSymmetricQuantiles.first > pValuesQuantileThreshold ) {
//          h1 = 0.5 * ( lBound + h0 );
//          rBound = h0;
//        } else if ( pValuesInSymmetricQuantiles.second > pValuesQuantileThreshold ) {
//          h1 = 0.5 * ( h0 + rBound );
//          lBound = h0;
//        }
//        numItInner++;
//        
//        if ( !this->_quietMode )
//          std::cerr << "Iteration " << numItOuter << "." << numItInner << ": h1=" << h1 << " in [" << lBound << ", " << rBound << "]; (pv<, pv>)=(" << pValuesInSymmetricQuantiles.first << ", " << pValuesInSymmetricQuantiles.second << "); pv_thres=" << pValuesQuantileThreshold << "; chi^2=" << chiSqr1 << std::endl;
//      } while ( aol::Abs<RealType> ( h1 - h0 ) / h1 > 1e-3 );
//      numItInner = 0;
//      numItOuter++;
//      numPixels *= 2;
//    } while ( numPixels < this->_input.size ( ) / 4 );
//    
//    if ( !this->_quietMode )
//      std::cerr << "Finished binary search. h_opt=" << h1 << std::endl;
//    
//    return h1;
//  }
};
  
} // namespace im

#endif