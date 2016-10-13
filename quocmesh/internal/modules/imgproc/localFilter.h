#ifndef __LOCALFILTER_H
#define __LOCALFILTER_H

#include <neighborhoodFilter.h>
#include <denoising.h>


namespace im {
  
class LOCALFILTER_TYPE {
public:
  static const int MOVING_AVERAGES = 0;  // Remove additive Gaussian white noise (AGWN)
  static const int BILATERAL       = 1;  // Remove Poisson noise
  
  static const int NUM = 2;
  
  static std::string getIdentifier ( const int FilterType ) {
    if ( FilterType == MOVING_AVERAGES ) return "MA";
    else if ( FilterType == BILATERAL ) return "BL";
    else throw aol::Exception ( "Did not recognize local filter type!", __FILE__, __LINE__ );
  }
  
  static std::string toString ( const int FilterType ) {
    if ( FilterType == MOVING_AVERAGES ) return "Moving averages";
    else if ( FilterType == BILATERAL ) return "Poisson noise";
    else throw aol::Exception ( "Did not recognize local filter type!", __FILE__, __LINE__ );
  }
};

  

template <typename _RealType>
struct LocalFilterOptions : public NHFilterOptions<_RealType, qc::ScalarArray<_RealType, qc::QC_2D>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D> > {
  int filterType;
  _RealType sigmaPosition, sigmaIntensity;
  
  LocalFilterOptions ( const std::string OutputDir = "", const bool Quiet = true )
    : NHFilterOptions<_RealType, qc::ScalarArray<_RealType, qc::QC_2D>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D> > ( OutputDir, Quiet ),
      filterType ( 0 ),
      sigmaPosition ( 0 ), sigmaIntensity ( 0 ) { }
  
  
  LocalFilterOptions ( const LocalFilterOptions<_RealType> &Options )
    : NHFilterOptions<_RealType, qc::ScalarArray<_RealType, qc::QC_2D>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D> > ( Options ),
      filterType ( Options.filterType ),
      sigmaPosition ( Options.sigmaPosition ), sigmaIntensity ( Options.sigmaIntensity ) { }
};
  
  
template <typename _RealType>
class LocalFilterTrait : public NHFilterTrait<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_2D>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D> > {
public:
  typedef LocalFilterOptions<_RealType> OptionsType;
};
  
  

/**
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType>
class LocalFilter
  : public NeighborhoodFilter<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_2D>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D>, LocalFilterTrait<_RealType> > {
                                
  typedef _RealType RealType;
  typedef qc::ScalarArray<RealType, qc::QC_2D> PictureType;
                                
  const static qc::Dimension Dim = qc::QC_2D;
  typedef qc::ScalarArray<_RealType, qc::QC_2D> BlockType;
  typedef qc::GridSize<qc::QC_2D> BlockSizeType;
  
  typedef LocalFilterOptions<RealType> OptionsType;
                                
  typedef LocalFilterTrait<RealType> FilterTrait;
protected:
  RealType _sigmaIntensity, _sigmaPosition;
public:
  LocalFilter<RealType> ( ) : _sigmaIntensity ( 0 ), _sigmaPosition ( 0 ) { }
    
  virtual void apply ( const OptionsType &Options, const PictureType &Arg, PictureType &Dest ) {
    this->_input = &Arg;
    this->_estimate = &Dest;
    
    OptionsType options ( Options );
    this->initialize ( options );
    this->preprocess ( options );
    setParameters ( options );
    
    denoise ( options );
    
    this->postprocess ( options );
    
    // (Optional) console output
    if ( !Options.quietMode && options.groundTruth != NULL && (*options.groundTruth).getNumX ( ) == this->_estimate->getNumX ( ) && (*options.groundTruth).getNumY ( ) == this->_estimate->getNumY ( ) ) {
      if ( options.noiseType == NOISE_TYPE::GAUSSIAN ) std::cerr << "PSNR = " << aol::PSNR<RealType> ( *options.groundTruth, *(this->_estimate), 255 ) << " dB" << std::endl;
      else if ( options.noiseType == NOISE_TYPE::POISSON ) std::cerr << "PSNR = " << aol::PSNR<RealType> ( *options.groundTruth, *(this->_estimate) ) << " dB" << std::endl;
    }
  }
                                
  virtual void setSimilaritySearchGrid ( const OptionsType &/*Options*/, const PictureType &/*Arg*/, const qc::CoordType &/*XRef*/,
                                         aol::RandomAccessContainer<qc::CoordType> &/*CoordsFinal*/, aol::RandomAccessContainer<qc::CoordType> &/*CoordsTMP*/ ) {
    throw aol::Exception ( "Local filters do not support similarity search!", __FILE__, __LINE__ );
  }
                                
  static std::string getMethodIdentifier ( const OptionsType &Options ) {
    return getMethodPrefix ( Options ) + "-" + LOCALFILTER_TYPE::getIdentifier ( Options.filterType );
  }
  
  static std::string getMethodPrefix ( const OptionsType &Options ) {
    std::stringstream ss;
    ss << NeighborhoodFilter<RealType, Dim, PictureType, BlockType, BlockSizeType, FilterTrait>::getMethodPrefix ( Options );
    if ( Options.filterType == LOCALFILTER_TYPE::BILATERAL ) {
      if ( Options.sigmaPosition > 0 )
        ss << "-" << Options.sigmaPosition;
      if ( Options.sigmaIntensity > 0 )
        ss << "-" << Options.sigmaIntensity;
    }
    return ss.str ( );
  }
                                
  static std::string name ( ) { return "LF"; }
protected:
  void denoise ( const OptionsType &Options ) {
    if ( Options.filterType == LOCALFILTER_TYPE::MOVING_AVERAGES )
      im::getMeanFilteredImage ( this->_preprocessedInput, *(this->_estimate), false, this->_blockSize[0] );
    else if ( Options.filterType == LOCALFILTER_TYPE::BILATERAL )
      im::bilateralFilter ( this->_preprocessedInput, *(this->_estimate), this->_blockSize[0], _sigmaPosition, _sigmaIntensity );
    else throw aol::Exception ( "Did not recognize local filter type!", __FILE__, __LINE__ );
  }
                                
  void setParameters ( const OptionsType &Options ) {
    this->_blockSize.setAll ( Options.blockSize[0] > 0 ? Options.blockSize[0] : 7 );
    
    if ( Options.filterType == LOCALFILTER_TYPE::BILATERAL ) {
      _sigmaPosition = ( Options.sigmaPosition > 0 ) ? Options.sigmaPosition : this->_blockSize[0];
      _sigmaIntensity = ( Options.sigmaIntensity > 0 ) ? Options.sigmaIntensity : 2 * this->_stdDev;
    }
  }
};

  
} // namespace im

#endif