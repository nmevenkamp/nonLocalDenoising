#ifndef __DENOISING_H
#define __DENOISING_H

#include <quoc.h>
#include <firstOrderTVAlgos.h>

namespace im {

/**
 * Based on http://stackoverflow.com/questions/5695865/bilateral-filter
 *
 * \author Berkels
 */
template <typename RealType>
void bilateralFilter ( const qc::ScalarArray<RealType, qc::QC_2D> &Input,
                       qc::ScalarArray<RealType, qc::QC_2D> &Output,
                       const int KernelRadius,
                       const RealType Sigma,
                       const RealType SigmaIntensity ) {
  const int numX = Input.getNumX();
  const int numY = Input.getNumY();
  for ( int y = 0; y < numY; ++y ) {
    for ( int x = 0; x < numX; ++x ) {

      RealType sumWeight = 0;
      RealType sum = 0;

      const RealType ctrPix = Input.get( x, y );

      const int kernelStartX = aol::Max ( x - KernelRadius, 0 );
      const int kernelEndX   = aol::Min ( x + KernelRadius, numX-1 );
      const int kernelStartY = aol::Max ( y - KernelRadius, 0 );
      const int kernelEndY   = aol::Min ( y + KernelRadius, numY-1 );

      for ( int j = kernelStartY; j <= kernelEndY; ++j ) {
        for ( int i = kernelStartX; i <= kernelEndX; ++i ) {

          const RealType curPix = Input.get( i, j );
          const RealType imageDist = sqrt ( static_cast<RealType> ( aol::Sqr ( i - x ) + aol::Sqr ( j - y ) ) );
          const RealType colorDist = sqrt ( static_cast<RealType> ( aol::Sqr ( curPix - ctrPix ) ) );

          const RealType currWeight = exp ( - aol::Sqr ( imageDist / Sigma ) * 0.5 ) * exp ( - aol::Sqr ( colorDist / SigmaIntensity ) * 0.5 );
          sumWeight += currWeight;

          sum += currWeight * curPix;
        }
      }
      Output.set ( x, y, sum / sumWeight );
    }
  }
}

/**
 * 1D version of the 2D bilateral filter above
 *
 * \author Berkels
 */
template <typename RealType>
void bilateralFilter ( const aol::Vector<RealType> &Input,
                       aol::Vector<RealType> &Output,
                       const int KernelRadius,
                       const RealType Sigma,
                       const RealType SigmaIntensity ) {
  const int length = Input.size();
  for ( int x = 0; x < length; ++x ) {
    
    RealType sumWeight = 0;
    RealType sum = 0;
    
    const RealType ctrPix = Input[x];
    
    const int kernelStartX = aol::Max ( x - KernelRadius, 0 );
    const int kernelEndX   = aol::Min ( x + KernelRadius, length-1 );
  
    for ( int i = kernelStartX; i <= kernelEndX; ++i ) {
      
      const RealType curPix = Input[i];
      const RealType imageDist = sqrt ( static_cast<RealType> ( aol::Sqr ( i - x ) ) );
      const RealType colorDist = sqrt ( static_cast<RealType> ( aol::Sqr ( curPix - ctrPix ) ) );
      
      const RealType currWeight = exp ( - aol::Sqr ( imageDist / Sigma ) * 0.5 ) * exp ( - aol::Sqr ( colorDist / SigmaIntensity ) * 0.5 );
      sumWeight += currWeight;
      
      sum += currWeight * curPix;
    }
    Output[x] = sum / sumWeight;
  }
}
  
  
/**
 * 1D median filter
 *
 * \author Mevenkamp
 */
template <typename RealType>
void medianFilter ( const aol::Vector<RealType> &Input,
                    aol::Vector<RealType> &Output,
                    const int FilterSize = 3 ) {
  const int Off = ( FilterSize - 1 ) >> 1;
  const int length = Input.size();
  for ( int x = 0; x < length; ++x ) {
    const int XMin = aol::Max ( 0, x - Off );
    const int XMax = aol::Min ( length - 1, x + Off );
    aol::Vector<RealType> tmp ( XMax - XMin + 1 );
    for ( int xx=XMin; xx<=XMax ; ++xx )
      tmp[xx-XMin] = Input[xx];
    Output[x] = tmp.getMedianValue ( );
  }
}
  

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class FirstOrderPrimalDualROFMinimizer : public qc::FirstOrderChambollePockTVAlgorithmType2<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  
private:
  const ArrayType &_image;
  
  void applyResolventOfDataTermSingle ( const RealType TauOverGamma, ArrayType &ArgDest ) const {
    ArgDest.addMultiple ( _image, TauOverGamma );
    ArgDest /= ( 1 + TauOverGamma );
  }
  
public:
  FirstOrderPrimalDualROFMinimizer ( const typename ConfiguratorType::InitType &Initializer,
                                    const RealType Gamma,
                                    const ArrayType &Image,
                                    const int MaxIterations = 1000,
                                    const RealType StopEpsilon = 0 )
  : qc::FirstOrderChambollePockTVAlgorithmType2<ConfiguratorType> ( Initializer, Gamma, MaxIterations, StopEpsilon ),
  _image ( Image ) {}
};
  
} // namespace im

#endif
