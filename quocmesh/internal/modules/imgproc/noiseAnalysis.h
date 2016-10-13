#ifndef __NOISEANALYSIS_H
#define __NOISEANALYSIS_H

#include <aol.h>
#include <atomFinder.h>
#include <regression.h>
#include <randomGenerator.h>
#include <derivativeFreeOptimization.h>
#include <gradientDescent.h>
#include <anscombe.h>
#include <probDistributionFunction.h>
#ifdef USE_MODULES_QT
#include <customPlotHandler.h>
#endif

namespace im {

template <typename RealType>
class RegularizedL1Dirichlet1DEnergy : public aol::Op<aol::Vector<RealType>, aol::Scalar<RealType> > {
protected:
  const aol::Vector<RealType> &_vals;
  const RealType _lambda;
  const RealType _epsSqr;
public:
  RegularizedL1Dirichlet1DEnergy ( const aol::Vector<RealType> &Vals, const RealType Lambda = 1, const RealType Eps = 1e-4 )
  : _vals ( Vals ),
  _lambda ( Lambda ),
  _epsSqr ( aol::Sqr<RealType> ( Eps ) ) { }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    if ( Arg.size ( ) != _vals.size ( ) ) throw aol::Exception ( "Argument dimension does not match signal length!", __FILE__, __LINE__ );
    
    Dest = 0;
    for ( int i=0; i<_vals.size ( ) ; ++i ) Dest += regularizedAbs ( Arg[i] - _vals[i] ) + 0.5 * _lambda * aol::Sqr<RealType> ( ( i < _vals.size ( ) - 1 ) ? Arg[i+1] - Arg[i] : Arg[0]-Arg[i] );
  }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::Scalar<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException ( "Not implemented", __FILE__, __LINE__ );
  }
protected:
  RealType regularizedAbs ( const RealType &Arg ) const {
    return sqrt ( aol::Sqr<RealType> ( Arg ) + _epsSqr );
  }
};

template <typename RealType>
class RegularizedL1Dirichlet1DDerivative : public aol::Op<aol::Vector<RealType> > {
protected:
  const aol::Vector<RealType> &_vals;
  const RealType _lambda;
  const RealType _epsSqr;
public:
  RegularizedL1Dirichlet1DDerivative ( const aol::Vector<RealType> &Vals, const RealType Lambda = 1, const RealType Eps = 1e-4 )
  : _vals ( Vals ),
  _lambda ( Lambda ),
  _epsSqr ( aol::Sqr<RealType> ( Eps ) ) { }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( Arg.size ( ) != _vals.size ( ) ) throw aol::Exception ( "Argument dimension does not match signal length!", __FILE__, __LINE__ );
    if ( Dest.size ( ) != Arg.size ( ) ) throw aol::Exception ( "Destination and argument dimensions do not match!", __FILE__, __LINE__ );
    
    Dest.setZero ( );
    for ( int j=0; j<Arg.size ( ) ; ++j ) {
      const RealType Argm1 = ( j > 0 ) ? Arg[j-1] : Arg[Arg.size ( ) -1];
      const RealType Argp1 = ( j < Arg.size ( ) - 1 ) ? Arg[j+1] : Arg[0];
      Dest[j] += ( Arg[j] - _vals[j] ) / regularizedAbs ( Arg[j] - _vals[j] ) + _lambda * ( 2 * Arg[j] - Argp1 - Argm1 );
    }
  }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::Vector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException ( "Not implemented", __FILE__, __LINE__ );
  }
protected:
  RealType regularizedAbs ( const RealType &Arg ) const {
    return sqrt ( aol::Sqr<RealType> ( Arg ) + _epsSqr );
  }
};
  
template <typename RealType>
static void regularizeSignalL1Dirichlet ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest,
                                          const RealType Lambda = 1 ) {
  RegularizedL1Dirichlet1DEnergy<RealType> E ( Arg, Lambda );
  RegularizedL1Dirichlet1DDerivative<RealType> DE ( Arg, Lambda );
  aol::GridlessGradientDescent<RealType, aol::Vector<RealType> > gradientDescentOp ( E, DE );
  gradientDescentOp.apply ( Arg, Dest );
}

template <typename RealType, typename PictureType>
struct doRegularizeScanLines {
  static void apply ( const PictureType &Arg, PictureType &Dest, const RealType Lambda = 1, const std::string &OutputDir = "" );
};

template <typename RealType>
struct doRegularizeScanLines<RealType, qc::ScalarArray<RealType, qc::QC_2D> > {
  static void apply ( const qc::ScalarArray<RealType, qc::QC_2D> &Arg, qc::ScalarArray<RealType, qc::QC_2D> &Dest, const RealType Lambda = 1, const std::string &OutputDir = "" ) {
#ifdef USE_MODULES_QT
    CustomPlotHandler<RealType> qcpHandler ( "Scan lines", true );
#endif
    int k = 0;
    for ( int y=0; y<Arg.getNumY ( ) ; ++y ) {
      aol::Vector<RealType> signal  ( Arg.getNumX ( ) ), regularizedSignal ( signal, aol::STRUCT_COPY );
      for ( int x=0; x<Arg.getNumX ( ) ; ++x ) signal[x] = Arg.get ( x, y );
      regularizeSignalL1Dirichlet<RealType> ( signal, regularizedSignal, Lambda );
      for ( int x=0; x<Arg.getNumX ( ) ; ++x ) Dest.set ( x, y, regularizedSignal[x] );
      ++k;
      
#ifdef USE_MODULES_QT
      if ( k < 10 && OutputDir != "" ) {
        qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( signal, aol::strprintf ( "line %d (raw)", k ).c_str ( ) );
        qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( regularizedSignal, aol::strprintf ( "line %d (regularized)", k ).c_str ( ) );
      }
#endif
    }
#ifdef USE_MODULES_QT
    if ( OutputDir != "" )
      qcpHandler.saveToFile ( aol::strprintf ( "%s/scanLines.q1cp", OutputDir.c_str ( ) ).c_str ( ) );
#endif
  }
};

template <typename RealType>
struct doRegularizeScanLines<RealType, qc::ScalarArray<RealType, qc::QC_3D> > {
  static void apply ( const qc::ScalarArray<RealType, qc::QC_3D> &Arg, qc::ScalarArray<RealType, qc::QC_3D> &Dest, const RealType Lambda = 1, const std::string &OutputDir = "" ) {
#ifdef USE_MODULES_QT
    CustomPlotHandler<RealType> qcpHandler ( "Scan lines", true );
#endif
    int k = 0;
    for ( int z=0; z<Arg.getNumZ ( ) ; ++z ) {
      for ( int y=0; y<Arg.getNumY ( ) ; ++y ) {
        aol::Vector<RealType> signal  ( Arg.getNumX ( ) ), regularizedSignal ( signal, aol::STRUCT_COPY );
        for ( int x=0; x<Arg.getNumX ( ) ; ++x ) signal[x] = Arg.get ( x, y, z );
        regularizeSignalL1Dirichlet<RealType> ( signal, regularizedSignal, Lambda );
        for ( int x=0; x<Arg.getNumX ( ) ; ++x ) Dest.set ( x, y, z, regularizedSignal[x] );
        ++k;
      
#ifdef USE_MODULES_QT
        if ( k < 10 && OutputDir != "" ) {
          qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( signal, aol::strprintf ( "line %d (raw)", k ).c_str ( ) );
          qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( regularizedSignal, aol::strprintf ( "line %d (regularized)", k ).c_str ( ) );
        }
#endif
      }
    }
#ifdef USE_MODULES_QT
    if ( OutputDir != "" )
      qcpHandler.saveToFile ( aol::strprintf ( "%s/scanLines.q1cp", OutputDir.c_str ( ) ).c_str ( ) );
#endif
  }
};
  
  
  
template <typename RealType>
class VSTVarianceCostFunction : public aol::Op<aol::Vector<RealType>, aol::Scalar<RealType> > {
protected:
  const aol::VectorContainer<aol::Vector<RealType> > &_samples;
  const RealType _eps;
  const aol::Vector<RealType> &_deltas;
public:
  VSTVarianceCostFunction ( const aol::VectorContainer<aol::Vector<RealType> > &Samples, const aol::Vector<RealType> &Deltas ) : _samples ( Samples ), _eps ( 0.125 ), _deltas ( Deltas ) { }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    if ( Arg.size ( ) != 2 ) throw aol::Exception ( "Expecting argument of size 2!", __FILE__, __LINE__ );
    
    aol::VectorContainer<aol::Vector<RealType> > samplesCorrected ( _samples );
    if ( _deltas.size ( ) == _samples.size ( ) ) {
      for ( int i=0; i<_samples.size ( ) ; ++i ) {
        const RealType variance = Arg[0] * _samples[i].getMeanValue ( ) + Arg[1];
        samplesCorrected[i] *= sqrt ( variance / ( variance + 1.0 / 12.0 * aol::Sqr<RealType> ( _deltas[i] ) ) );
        samplesCorrected[i].addToAll ( _samples[i].getMeanValue ( ) - samplesCorrected[i].getMeanValue ( ) );
      }
    }
    aol::VectorContainer<aol::Vector<RealType> > VSTSamples ( samplesCorrected );
    GeneralizedAnscombeForward<RealType> GATForward ( Arg[0], Arg[1] );
    for ( int i=0; i<samplesCorrected.size ( ) ; ++i ) GATForward.apply ( samplesCorrected[i], VSTSamples[i] );
    aol::Vector<RealType> VSTSampleResiduals ( samplesCorrected.size ( ) );
    for ( int i=0; i<samplesCorrected.size ( ) ; ++i ) VSTSampleResiduals[i] = aol::Abs ( getNonLinearResponse ( VSTSamples[i].getVariance ( ) - 1.0 ) );
    Dest = VSTSampleResiduals.getMeanValue ( );
  }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::Scalar<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException ( "Not implemented", __FILE__, __LINE__ );
  }
protected:
  RealType getNonLinearResponse ( const RealType Arg ) const {
    return ( aol::Abs<RealType> ( Arg ) < _eps ) ? 0.5 / _eps * aol::Sqr<RealType> ( Arg ) : Arg - 0.5 * _eps;
  }
};
  
template <typename RealType>
class NoiseVarianceNegativeLogLikelihoodFunction : public aol::Op<aol::Vector<RealType>, aol::Scalar<RealType> > {
protected:
  const aol::VectorContainer<aol::Vector<RealType> > &_samples;
  const std::vector<std::pair<RealType, RealType> > &_priorProbabilities;
  const RealType _epsilon;
  int _N;
  aol::Vector<RealType> _means, _stdDevs, _c, _d, _coeff;
public:
  NoiseVarianceNegativeLogLikelihoodFunction ( const aol::VectorContainer<aol::Vector<RealType> > &Samples,
                                               const std::vector<std::pair<RealType, RealType> > &PriorProbabilities,
                                               const RealType Epsilon = 1e-3 )
    : _samples ( Samples ),
      _priorProbabilities ( PriorProbabilities ),
      _epsilon ( Epsilon ),
      _means ( Samples.size ( ) ), _stdDevs ( Samples.size ( ) ),
      _c ( Samples.size ( ) ), _d ( Samples.size ( ) ), _coeff ( Samples.size ( ) ) {
    for ( int i=0; i<_samples.size ( ) ; ++i ) {
      const RealType ni = static_cast<RealType> ( _samples[i].size ( ) );
      _means[i] = _samples[i].getMeanValue ( );
      _stdDevs[i] = _samples[i].getStdDevUnbiased ( );
      _c[i] = 1.0 / ni;
      _d[i] = 0.5 / ni + 0.625 / aol::Sqr<RealType> ( ni );
      _coeff[i] = 1.0 / ( 2.0 * aol::NumberTrait<RealType>::pi * sqrt ( _c[i] * _d[i] ) );
    }
  }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    if ( Arg.size ( ) != 2 ) throw aol::Exception ( "Expecting argument of size 2!", __FILE__, __LINE__ );
    
    Dest = 0.0;
    for ( int i=0; i<_samples.size ( ) ; ++i ) {
      RealType integral = 0.0;
      for ( unsigned int j=0; j<_priorProbabilities.size ( ) ; ++j ) {
        const RealType y = _priorProbabilities[j].first;
        const RealType sigmaSqrReg = aol::Max<RealType> ( _epsilon, Arg[0] * y + Arg[1] ), sigmaReg = sqrt ( sigmaSqrReg );
        const RealType pGivenY = _coeff[i] / sigmaSqrReg * exp ( -0.5 / sigmaSqrReg * ( aol::Sqr<RealType> ( _means[i] - y ) / _c[i] + aol::Sqr<RealType> ( _stdDevs[i] - sigmaReg ) / _d[i] ) );
        integral += pGivenY * _priorProbabilities[i].second;
      }
      integral *= ( _priorProbabilities[_priorProbabilities.size ( )-1].first - _priorProbabilities[0].first ) / static_cast<RealType> ( _priorProbabilities.size ( ) );
      Dest -= log ( integral );
    }
  }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::Scalar<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException ( "Not implemented", __FILE__, __LINE__ );
  }
};

template <typename RealType>
RealType F ( const RealType Lambda, const int K,
             const RealType Alpha, const RealType Sigma, const RealType Mu,
             const RealType R ) {
  return -aol::NormalDistribution<RealType>::PDF ( R - Alpha * K, Mu, aol::Sqr<RealType> ( Sigma ) ) * aol::PoissonDistribution<RealType>::CDF ( K, aol::Max<RealType> ( Lambda, 1e-3 ) );
}

template <typename RealType>
RealType DeltaLikelihood ( const aol::Vector<RealType> &Realizations, const RealType Delta,
                           const RealType Alpha, const RealType Sigma, const RealType Mu ) {
  const RealType lambda = ( Realizations.getMeanValue ( ) - Mu ) / Alpha;
  
  RealType logLikelihood = 0;
  for ( int i=0; i<Realizations.size ( ) ; ++i ) {
    RealType pdf = 0, summand = 1;
    int k = 0;
    while ( k < 50 ) {
      summand = F ( lambda + 0.5 * Delta, k, Alpha, Sigma, Mu, Realizations[i] ) - F ( lambda - 0.5 * Delta, k, Alpha, Sigma, Mu, Realizations[i] );
      pdf += summand;
      ++k;
    }
    pdf /= Delta;
    logLikelihood += log ( pdf );
  }
  return -logLikelihood;
}
  
template <typename RealType>
RealType MuLikelihood ( const aol::Vector<RealType> &Realizations, const RealType Mu,
                        const RealType a, const RealType b,
                        const RealType Lambda = 0 ) {
  const RealType lambda = ( Lambda > 0 ) ? Lambda : ( Realizations.getMeanValue ( ) - Mu ) / a;
  if ( lambda <= 0 || b + a * Mu <= 0 ) return aol::NumberTrait<RealType>::Inf;
  
  const RealType sigma = sqrt ( b + a * Mu );
  
  RealType logLikelihood = 0;
  for ( int i=0; i<Realizations.size ( ) ; ++i )
    logLikelihood += log ( aol::MixedPoissonGaussianDistribution<RealType>::PDF ( Realizations[i], lambda, a, Mu, sigma ) );
  return -logLikelihood;
}
  
template <typename RealType>
RealType MuLikelihood ( const aol::VectorContainer<aol::Vector<RealType> > &Samples, const RealType Mu,
                        const RealType a, const RealType b,
                        const aol::Vector<RealType> &Lambdas = aol::Vector<RealType> ( ) ) {
  aol::Vector<RealType> likelihoods ( Samples.size ( ) );
  for ( int i=0; i<Samples.size ( ) ; ++i ) likelihoods[i] = MuLikelihood ( Samples[i], Mu, a, b, Lambdas[i] );
  return likelihoods.getMeanValue ( );
}
  
template <typename RealType>
class MuLikelihoodCostFunction : public aol::Op<aol::Vector<RealType>, aol::Scalar<RealType> > {
protected:
  const aol::VectorContainer<aol::Vector<RealType> > &_samples;
  const RealType _A, _B;
  RealType _muMin, _muMax;
public:
  MuLikelihoodCostFunction ( const aol::VectorContainer<aol::Vector<RealType> > &Samples, const RealType A, const RealType B )
    : _samples ( Samples ), _A ( A ), _B ( B ) {
    const int m = Samples.size ( );
    aol::Vector<RealType> meanVals ( m );
    for ( int i=0; i<m ; ++i ) meanVals[i] = Samples[i].getMeanValue ( );
    _muMin = -B / A;
    _muMax = meanVals.getMinValue ( );
  }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    if ( Arg.size ( ) != 1 ) throw aol::Exception ( "Expecting argument of size 1!", __FILE__, __LINE__ );
    
    if ( Arg[0] < _muMin || Arg[0] > _muMax )
      Dest = aol::NumberTrait<RealType>::Inf;
    else
      Dest = MuLikelihood ( _samples, Arg[0], _A, _B );
  }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::Scalar<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException ( "Not implemented", __FILE__, __LINE__ );
  }
};
  
template <typename RealType>
RealType NoiseLikelihood ( const aol::VectorContainer<aol::Vector<RealType> > &Samples,
                           const RealType Alpha, const RealType Sigma, const RealType Mu,
                           const aol::Vector<RealType> &Deltas ) {
  aol::Vector<RealType> negLogLikelihoods ( Samples.size ( ) );
  for ( int i=0; i<Samples.size ( ) ; ++i ) {
    const RealType lambda = ( Samples[i].getMeanValue ( ) - Mu ) / Alpha;
    if ( lambda <= 0 || Sigma <= 0 ) return aol::NumberTrait<RealType>::Inf;
    
    RealType logLikelihood = 0;
    for ( int j=0; j<Samples[i].size ( ) ; ++j ) {
      RealType pdf = 0, summand = 1;
      int k = ceil ( lambda );
//      int k = 0;
      do {
        summand = F ( lambda + 0.5 * Deltas[i], k, Alpha, Sigma, Mu, Samples[i][j] ) - F ( lambda - 0.5 * Deltas[i], k, Alpha, Sigma, Mu, Samples[i][j] );
        pdf += summand;
        ++k;
      } while ( summand  > 1e-300 );
      k = ceil ( lambda ) - 1;
      do {
        summand = F ( lambda + 0.5 * Deltas[i], k, Alpha, Sigma, Mu, Samples[i][j] ) - F ( lambda - 0.5 * Deltas[i], k, Alpha, Sigma, Mu, Samples[i][j] );
        pdf += summand;
        --k;
      } while ( k >= 0 && summand  > 1e-300 );
      logLikelihood += log ( pdf );
    }
    negLogLikelihoods[i] = -logLikelihood;
  }
  
  return negLogLikelihoods.getMeanValue ( );
}
  
template <typename RealType>
class NoiseCostFunction : public aol::Op<aol::Vector<RealType>, aol::Scalar<RealType> > {
protected:
  const aol::VectorContainer<aol::Vector<RealType> > &_samples;
  const aol::Op<aol::Vector<RealType> > *_deltaMap;
public:
  NoiseCostFunction ( const aol::VectorContainer<aol::Vector<RealType> > &Samples, const aol::Op<aol::Vector<RealType> > *DeltaMap = NULL )
    : _samples ( Samples ),
      _deltaMap ( DeltaMap ) { }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    aol::Vector<RealType> deltas ( _samples.size ( ) );
    if ( _deltaMap != NULL ) {
      aol::Vector<RealType> deltaMapArg ( Arg.size ( ) - 3 );
      for ( int i=3; i<Arg.size ( ) ; ++i ) deltaMapArg[i-3] = Arg[i];
      _deltaMap->apply ( deltaMapArg, deltas );
    } else {
      for ( int i=0; i<deltas.size ( ) ; ++i ) deltas[i] = Arg[i+3];
    }
    Dest = NoiseLikelihood ( _samples, Arg[0], Arg[1], Arg[2], deltas );
  }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::Scalar<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException ( "Not implemented", __FILE__, __LINE__ );
  }
};
  
template <typename RealType>
class DeltaPolynomialMap : public aol::Op<aol::Vector<RealType> > {
protected:
  const aol::Vector<RealType> &_lambdas;
  const int _polyDegree;
public:
  DeltaPolynomialMap ( const aol::Vector<RealType> &Lambdas, const int PolyDegree = 2 )
    : _lambdas ( Lambdas ),
      _polyDegree ( PolyDegree ) { }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( Arg.size ( ) != _polyDegree + 1 ) throw aol::Exception ( "Argument size not equal to dof of specified polynomial!", __FILE__, __LINE__ );
    if ( Dest.size ( ) != _lambdas.size ( ) ) throw aol::Exception ( "Destiation size not equal to number of mean values to be mapped!", __FILE__, __LINE__ );
    Dest.setZero ( );
    for ( int i=0; i<_lambdas.size ( ) ; ++i )
      for ( int j=0; j<=_polyDegree ; ++j ) Dest[i] += Arg[j] * pow ( _lambdas[i], j );
  }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::Vector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException ( "Not implemented", __FILE__, __LINE__ );
  }
};

  
/**
 * \brief Noise analysis based on level sets (noisy samples belonging roughly to the same mean)
 *
 *        Currently supports CMOS noise model: \f$z = \alpha p + \eta, p \sim \text{Pois}(\lambda), \eta \sim \mathcal{N}(\mu, \sigma^2)\f$
 *        Note: the CMOS noise model implies a linear variance: \f$\text{Var}(z) = a \mathbb{E}[u] + b, a = \alpha, b = \sigma^2 - \alpha \mu \f$
 * \author mevenkamp
 * \ingroup imgproc
 */
template <typename RealType>
class NoiseAnalyzer {
public:
  static void getGaussianVarianceParams ( const aol::VectorContainer<aol::Vector<RealType> > &Samples, RealType &Sigma,
                                          const bool QuietMode = true, const std::string &OutputDir = "",
                                          const int NumSamplesMin = 0 ) {
    const RealType trueVariance = aol::Sqr<RealType> ( Sigma );
    
    aol::VectorContainer<aol::Vector<RealType> > samples;
    dropSamples ( Samples, samples, NumSamplesMin );
    const int m = samples.size ( );
    
    // Solve explicitly for sigma in the least-squares sense => mean value
    aol::Vector<RealType> varianceEstimates ( m );
    for ( int i=0; i<m ; ++i )
      varianceEstimates[i] = samples[i].getVariance ( );
    const RealType varianceLS = varianceEstimates.getMeanValue ( );
    Sigma = sqrt ( varianceLS );
    if ( !QuietMode ) std::cerr << "NoiseAnalyzer: least-squares solution: sigma = " << Sigma << std::endl;
    
    if ( OutputDir != "" ) {
      std::ofstream txtFile ( aol::strprintf ( "%s/noiseAnalysis.txt", OutputDir.c_str ( ) ).c_str ( ) );
      txtFile << "Least-squares solution" << std::endl;
      txtFile << "sigma = " << Sigma << std::endl;
      txtFile.close ( );
      
      // Prepare different data structures of the estimates for convenience
      aol::Vector<RealType> meanEstimates ( m );
      std::vector<std::pair<RealType, RealType> > fittedVariances, fittedVariancesTrue;
      for ( int i=0; i<m ; ++i ) {
        meanEstimates[i] = samples[i].getMeanValue ( );
        fittedVariances.push_back ( std::pair<RealType, RealType> ( meanEstimates[i], varianceLS ) );
        if ( trueVariance != 0 )
          fittedVariancesTrue.push_back ( std::pair<RealType, RealType> ( meanEstimates[i], trueVariance ) );
      }
      
#ifdef USE_MODULES_QT
      // Plot estimated variances versus the respective mean values of the signal
      {
        CustomPlotHandler<RealType> qcpHandler ( "Noise analysis based on raw samples", true );
        qcpHandler.setAxesLabels ( "mean", "variance" );
        qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( meanEstimates, varianceEstimates, "samples", QuocQCPGraphStyle<RealType> ( "", -1, -1, 0, 8, 3 ) );
        qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( fittedVariances, "least squares fit" );
        if ( trueVariance != 0 ) qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( fittedVariancesTrue, "ground truth" );
        qcpHandler.saveToFile ( aol::strprintf ( "%s/noiseAnalysis.q1cp", OutputDir.c_str ( ) ).c_str ( ) );
      }
#endif
    }
  }
  
  static void getCMOSVarianceParams ( const aol::VectorContainer<aol::Vector<RealType> > &Samples, RealType &Alpha, RealType &SigmaSqr,
                                      const bool QuietMode = true, const std::string &OutputDir = "",
                                      const int NumSamplesMin = 0,
                                      const aol::Vector<RealType> &Deltas = aol::Vector<RealType> ( ),
                                      const std::vector<std::pair<RealType, RealType> > &/*PriorProbabilities*/ = std::vector<std::pair<RealType, RealType> > ( ) ) {
    const RealType trueAlpha = Alpha;
    const RealType trueSigmaSqr = SigmaSqr;
    
    aol::VectorContainer<aol::Vector<RealType> > samples;
    dropSamples ( Samples, samples, NumSamplesMin );
    const int m = samples.size ( );
    
    // Fit model var(mu) = a mu + b to pairs (mu_i,sigma_i^2) using least squares
    aol::FullMatrix<RealType> A ( m, 2 );
    aol::Vector<RealType> b ( m ), xLS ( 2 );
    for ( int i=0; i<m ; ++i ) {
      A.set ( i, 0, samples[i].getMeanValue ( ) );
      A.set ( i, 1, 1.0 );
      b[i] = samples[i].getVariance ( ) - ( ( Deltas.size ( ) == Samples.size ( ) ) ?  1.0 / 12.0 * aol::Sqr<RealType> ( Deltas[i] ) : 0 );
    }
    aol::LinearRegressionNormalEquations<RealType> normEquations ( A );
    normEquations.apply ( b, xLS );
    if ( !QuietMode ) std::cerr << "NoiseAnalyzer: least-squares solution: a = " << xLS[0] << "; b = " << xLS[1] << std::endl;
    
    aol::Vector<RealType> abVSTVarMin ( 2 );
    VSTVarianceCostFunction<RealType> f ( samples, Deltas );
    aol::NelderMeadDownhillSimplexAlgorithm<RealType> optimizationAlg ( f, 1e-16, 100, QuietMode );
    abVSTVarMin = xLS;
    optimizationAlg.solve ( abVSTVarMin );
    if ( !QuietMode ) std::cerr << "NoiseAnalyzer: vst-variance equalization solution: a = " << abVSTVarMin[0] << "; b = " << abVSTVarMin[1] << std::endl;
    
    Alpha = abVSTVarMin[0];
    SigmaSqr = abVSTVarMin[1];
    
    if ( !QuietMode ) std::cerr << "Alpha = " << abVSTVarMin[0] << std::endl;
    
    if ( OutputDir != "" ) {
      std::ofstream txtFile ( aol::strprintf ( "%s/noiseAnalysis.txt", OutputDir.c_str ( ) ).c_str ( ) );
      txtFile << "Least-squares solution" << std::endl;
      txtFile << "a = " << xLS[0] << std::endl << "b = " << xLS[1] << std::endl << std::endl;
      txtFile << "VST variance equalization solution" << std::endl;
      txtFile << "a = " << Alpha << std::endl << "b = " << SigmaSqr << std::endl;
      txtFile.close ( );
      
      // Prepare different data structures of the estimates for convenience
      aol::Vector<RealType> meanEstimates ( samples.size ( ) ), varianceEstimates ( samples.size ( ) ), varianceEstimatesCorrected ( samples.size ( ) );
      for ( int i=0; i<samples.size ( ) ; ++i ) {
        meanEstimates[i] = samples[i].getMeanValue ( );
        varianceEstimates[i] = samples[i].getVariance ( );
        if ( Deltas.size ( ) == samples.size ( ) )
          varianceEstimatesCorrected[i] = varianceEstimates[i] - 1.0 / 12.0 * aol::Sqr<RealType> ( Deltas[i] );
      }
      
      // Sample least-squares, maximum-likelihood, and VST-variance minimization variance models
      const int numSamplePoints = 100;
      const RealType minMu = meanEstimates.getMinValue ( ), maxMu = meanEstimates.getMaxValue ( );
      std::vector<std::pair<RealType, RealType> > fittedVariancesLS, fittedVariancesVSTVarMin, fittedVariancesTrue;
      for ( int i=0; i<=numSamplePoints ; ++i ) {
        const RealType mui = minMu + i / static_cast<RealType> ( numSamplePoints ) * ( maxMu - minMu );
        fittedVariancesLS.push_back ( std::pair<RealType, RealType> ( mui, xLS[0] * mui + xLS[1] ) );
        fittedVariancesVSTVarMin.push_back ( std::pair<RealType, RealType> ( mui, abVSTVarMin[0] * mui + abVSTVarMin[1] ) );
        if ( trueAlpha != 0 ) fittedVariancesTrue.push_back ( std::pair<RealType, RealType> ( mui, trueAlpha * mui + trueSigmaSqr ) );
      }
#ifdef USE_MODULES_QT
      // Plot estimated variances versus the respective mean values of the signal
      {
        CustomPlotHandler<RealType> qcpHandler ( "Noise analysis based on raw samples", true );
        qcpHandler.setAxesLabels ( "mean", "variance" );
        qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( meanEstimates, varianceEstimates, "samples", QuocQCPGraphStyle<RealType> ( "", -1, -1, 0, 8, 3 ) );
        if ( Deltas.size ( ) == samples.size ( ) )
          qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( meanEstimates, varianceEstimatesCorrected, "samples (corr.)", QuocQCPGraphStyle<RealType> ( "", -1, -1, 0, 8, 3 ) );
        qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( fittedVariancesLS, "least squares fit" );
        qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( fittedVariancesVSTVarMin, "VST optimization fit" );
        if ( trueAlpha != 0 ) qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( fittedVariancesTrue, "ground truth" );
        qcpHandler.saveToFile ( aol::strprintf ( "%s/noiseAnalysis.q1cp", OutputDir.c_str ( ) ).c_str ( ) );
      }
      
      // Plot variances of the generalized Anscombe transformed samples
      {
        aol::VectorContainer<aol::Vector<RealType> > samplesVST ( samples );
        GeneralizedAnscombeForward<RealType> GATForward ( Alpha, SigmaSqr );
        for ( int i=0; i<samples.size ( ) ; ++i ) GATForward.apply ( samples[i], samplesVST[i] );
        std::vector<std::pair<RealType, RealType> > vstSamples;
        for ( int i=0; i<samples.size ( ) ; ++i ) vstSamples.push_back ( std::pair<RealType, RealType> ( samplesVST[i].getMeanValue ( ), samplesVST[i].getVariance ( ) ) );
        CustomPlotHandler<RealType> qcpHandler ( "Variance-stabilized samples" );
        qcpHandler.setAxesLabels ( "mean", "variance" );
        qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( vstSamples, "", QuocQCPGraphStyle<RealType> ( "", -1, -1, 0, 8, 3 ) );
        qcpHandler.saveToFile ( aol::strprintf ( "%s/samplesVST.q1cp", OutputDir.c_str ( ) ).c_str ( ) );
      }
      
      if ( Deltas.size ( ) == Samples.size ( ) ) {
        CustomPlotHandler<RealType> qcpHandler ( "Deltas", true );
        qcpHandler.setAxesLabels ( "mean", "delta" );
        qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( meanEstimates, Deltas, "", QuocQCPGraphStyle<RealType> ( "", -1, -1, 0, 8, 3 ) );
        qcpHandler.saveToFile ( aol::strprintf ( "%s/deltas.q1cp", OutputDir.c_str ( ) ).c_str ( ) );
      }
#endif
    }
  }
  
  static void getCMOSPedestalParam ( const aol::VectorContainer<aol::Vector<RealType> > &Samples, RealType &Mu,
                                     const RealType A, const RealType /*B*/,
                                     const bool QuietMode = true, const std::string &OutputDir = "",
                                     const int NumSamplesMin = 0 ) {
    aol::VectorContainer<aol::Vector<RealType> > samples;
    dropSamples ( Samples, samples, NumSamplesMin );
    const int m = samples.size ( );

    const RealType ASqr = aol::Sqr<RealType> ( A );
    
    Mu = 0;
    for ( int i=0; i<m ; ++i ) {
      const int n = samples[i].size ( );
      const RealType m1 = samples[i].getMeanValue ( );
      RealType m3 = 0;
      for ( int j=0; j<n ; ++j )
        m3 += aol::Sqr<RealType> ( samples[i][j] - m1 ) * ( samples[i][j] - m1 );
      m3 /= static_cast<RealType> ( n );
      Mu += m1 - m3 / ASqr;
    }
    Mu /= static_cast<RealType> ( m );
    
    if ( !QuietMode ) std::cerr << "Mu = " << Mu << std::endl;
    
    if ( OutputDir != "" ) {
      std::ofstream txtFile;
      txtFile.open ( aol::strprintf ( "%s/noiseAnalysis.txt", OutputDir.c_str ( ) ).c_str ( ), std::ofstream::app );
      txtFile << std::endl << "Method of moments solution for mu" << std::endl;
      txtFile << "mu = " << Mu << std::endl;
      txtFile.close ( );
    }
  }
  
  static void getCMOSGaussianStdDevParam ( RealType &Sigma, const RealType A, const RealType B, const RealType Mu,
                                           const bool QuietMode = true, const std::string &OutputDir = "" ) {
    Sigma = sqrt ( aol::Max<RealType> ( 0, B + Mu * A ) );
    
    if ( !QuietMode ) std::cerr << "Sigma = " << Sigma << std::endl;
    
    if ( OutputDir != "" ) {
      std::ofstream txtFile;
      txtFile.open ( aol::strprintf ( "%s/noiseAnalysis.txt", OutputDir.c_str ( ) ).c_str ( ), std::ofstream::app );
      txtFile << std::endl << "Direct solution for sigma" << std::endl;
      txtFile << "sigma = " << Sigma << std::endl;
      txtFile.close ( );
    }
  }
  
  template <typename PictureType>
  static void estimateDelta ( RealType &Delta,
                              const aol::Vector<int> &SampleCoords,
                              const PictureType &SmoothInput ) {
    aol::Vector<RealType> sampleEstimate ( SampleCoords.size ( ) );
    for ( int i=0; i<SampleCoords.size ( ) ; ++i ) sampleEstimate[i] = SmoothInput[SampleCoords[i]];
    aol::Vec2<RealType> minMax = sampleEstimate.getSaturatedMinMaxValue ( 10 );
    Delta = minMax[1] - minMax[0];
  }
  
  template <typename PictureType>
  static void estimateDeltas ( aol::Vector<RealType> &Deltas,
                               const aol::VectorContainer<aol::Vector<int> > &SamplesCoords,
                               const PictureType &RAWInput,
                               const std::string &OutputDir = "" ) {
    PictureType normalizedInput ( RAWInput );
    normalizedInput.scaleValuesTo01 ( );
    PictureType regularizedInput ( normalizedInput, aol::STRUCT_COPY );
    doRegularizeScanLines<RealType, PictureType>::apply ( normalizedInput, regularizedInput, 100, OutputDir );
    Deltas.resize ( SamplesCoords.size ( ) );
    aol::ProgressBar<> progressBar ( "Estimating deltas", std::cerr );
    progressBar.start ( SamplesCoords.size ( ) );
    for ( int i=0; i<SamplesCoords.size ( ) ; ++i, progressBar++ ) estimateDelta ( Deltas[i], SamplesCoords[i], RAWInput );
    progressBar.finish ( );
    Deltas.addToAll ( -Deltas.getMinValue ( ) );
  }

  static void dropSamples ( const aol::VectorContainer<aol::Vector<RealType> > &Samples, aol::VectorContainer<aol::Vector<RealType> > &SamplesReduced,
                            const int NumSamplesMin ) {
    SamplesReduced.clear ( );
    for ( int i=0; i<Samples.size ( ) ; ++i ) {
      if ( Samples[i].size ( ) >= NumSamplesMin )
        SamplesReduced.pushBack ( Samples[i] );
    }
    
    if ( SamplesReduced.size ( ) == 0 ) {
      aol::Vector<int> sizes;
      Samples.getSizes ( sizes );
      throw aol::Exception ( aol::strprintf ( "NoiseAnalyzer: no statistically significant samples found (max size = %d )", ( sizes.size ( ) > 0 ) ? sizes.getMaxValue ( ) : 0 ).c_str ( ), __FILE__, __LINE__ );
    }
  }
};
  
} // namespace im

#endif
