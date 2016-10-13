#ifndef __BLOCKOPERATORS_H
#define __BLOCKOPERATORS_H

#include <array.h>
#include <convolution.h>
#ifdef USE_BOOST
#ifndef Q_MOC_RUN // See: https://bugreports.qt-project.org/browse/QTBUG-22829
#include <boost/math/special_functions/bessel.hpp>
#endif
#endif
#include <blockStructures.h>
#include <matrixInverse.h>

/*
 * 2D windows
 */
// beta = pi * alpha
template <typename _RealType>
void setKaiserWindow ( qc::ScalarArray<_RealType, qc::QC_2D> &Window, const _RealType Beta ) {
  aol::Vec2<short> size ( Window.getNumX ( ), Window.getNumY ( ) );
  aol::Vector<_RealType> kaiserX ( size[0] ), kaiserY ( size[1] );
#ifdef USE_BOOST
  if ( size[0] > 1 )
    for ( int i=0; i<size[0] ; ++i )
      kaiserX[i] = static_cast<_RealType> ( boost::math::cyl_bessel_i<int, _RealType> ( 0, Beta * sqrt ( 1 - aol::Sqr<_RealType> ( 2.0 * i / ( size[0] - 1 ) - 1 ) ) )
                                          / boost::math::cyl_bessel_i<int, _RealType> ( 0, Beta ) );
  else kaiserX[0] = 1.0;
  if ( size[1] > 1 )
    for ( int i=0; i<size[1] ; ++i )
      kaiserY[i] = static_cast<_RealType> ( boost::math::cyl_bessel_i<int, _RealType> ( 0, Beta * sqrt ( 1 - aol::Sqr<_RealType> ( 2.0 * i / ( size[1] - 1 ) - 1 ) ) )
                                          / boost::math::cyl_bessel_i<int, _RealType> ( 0, Beta ) );
  else kaiserY[0] = 1.0;
#else
    throw aol::Exception ( "Boost required! Compile with -DUSE_BOOST=1", __FILE__, __LINE__ );
#endif
  
  for ( int x=0; x<size[0] ; ++x )
    for ( int y=0; y<size[1] ; ++y )
      Window.set ( x, y, kaiserX[x] * kaiserY[y] );
}


template <typename _RealType>
void setGaussianWindow ( qc::ScalarArray<_RealType, qc::QC_2D> &Window ) {
  if ( Window.getNumXYZ ( ) % 2 == 0 )
    throw aol::Exception ( "Window size must be odd!", __FILE__, __LINE__ );

  const short offset = ( Window.getNumXYZ ( ) - 1 ) / 2;
  const _RealType stdDev = 0.5 * Window.getNumXYZ ( ), c1 = 1.0 / ( stdDev * sqrt ( 2 * aol::NumberTrait<_RealType>::pi ) ), c2 = 2 * aol::Sqr<_RealType> ( stdDev );
  for ( int y=-offset; y<=offset ; ++y )
    for ( int x=-offset; x<=offset ; ++x )
      Window.set ( x + offset, y + offset, c1 * exp ( -( aol::Pow ( x, 2 ) + aol::Pow ( y, 2 ) ) / c2 ) );
}


/*
 * Block transform operators
 */
template <typename _RealType>
class MatrixFormTransformation {
  typedef _RealType RealType;
protected:
  aol::FullMatrix<RealType, false> _mTForward, _mTInverse;
public:
  MatrixFormTransformation ( const int SignalLength, const std::string &Transform ) {
    if ( SignalLength >= 2 ) {
      // Set matrix coefficients of forward transform
      _mTForward.reallocate ( SignalLength, SignalLength );
      if ( Transform == "dct" ) {
        for ( short row=0; row<SignalLength; ++row )
          for ( short col=0; col<SignalLength ; ++col )
            _mTForward.set ( row, col, cos ( aol::NumberTrait<RealType>::pi * ( 2 * col + 1 ) * row / ( 2 * SignalLength ) ) );
      } else if ( qc::DWT::isDWT ( Transform ) ) {
#ifdef USE_LIB_WAVELET
        for ( short col=0; col<SignalLength; ++col ) {
          std::vector<double> signal ( SignalLength ), dwtOutput, flags;
          signal[col] = 1.0;
          dwt ( signal, log2 ( SignalLength ), Transform, dwtOutput, flags );
          for ( short row=0; row<SignalLength ; ++row )
            _mTForward.set ( row, col, dwtOutput[row] );
        }
#else
        throw aol::Exception ( "External library wavelet1d required!", __FILE__, __LINE__ );
#endif
      } else throw aol::Exception ( "Did not recognize transformation!", __FILE__, __LINE__ );
      
      // Normalize basis elements (normalize 2-Norm of the rows of the forward transform)
      RealType rowNormInv;
      for ( short row=0; row<SignalLength ; ++row ) {
        rowNormInv = 0;
        for ( short col=0; col<SignalLength ; ++col ) rowNormInv += aol::Sqr<RealType> ( _mTForward.get ( row, col ) );
        rowNormInv = 1 / sqrt ( rowNormInv );
        for ( short col=0; col<SignalLength ; ++col ) _mTForward.set ( row, col, _mTForward.get ( row, col ) * rowNormInv );
      }
      
      // Compute inverse transform matrix
      _mTInverse.reallocate ( SignalLength, SignalLength );
      aol::QRInverse<RealType> qrInverse ( _mTForward );
      _mTInverse = qrInverse.getFM ( );
    } else {
      _mTForward.reallocate ( 1, 1 );
      _mTForward.set ( 0, 0, 1.0 );
      _mTInverse.reallocate ( 1, 1 );
      _mTInverse.set ( 0, 0, 1.0 );
    }
  }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest, enum qc::FourierTransformDirection Direction ) const {
    if ( Direction == qc::FTForward ) _mTForward.apply ( Arg, Dest );
    else _mTInverse.apply ( Arg, Dest );
  }
  
  void apply ( const qc::ScalarArray<RealType, qc::QC_2D> &Arg, qc::ScalarArray<RealType, qc::QC_2D> &Dest, enum qc::FourierTransformDirection Direction ) const {
    const int nxy = Arg.getNumXYZ ( );
    aol::Vector<RealType> vecArg ( nxy ), vecDest ( nxy );
    for ( int y=0; y<nxy ; ++y ) {
      for ( int x=0; x<nxy ; ++x ) vecArg[x] = Arg.get ( x, y );
      apply ( vecArg, vecDest, Direction );
      for ( int x=0; x<nxy ; ++x ) Dest.set ( x, y, vecDest[x] );
    }
    for ( int x=0; x<nxy ; ++x ) {
      for ( int y=0; y<nxy ; ++y ) vecArg[y] = Dest.get ( x, y );
      apply ( vecArg, vecDest, Direction );
      for ( int y=0; y<nxy ; ++y ) Dest.set ( x, y, vecDest[y] );
    }
  }
};


template <typename _RealType, typename _BlockType> struct DimensionTrait { static const qc::Dimension Dim = _BlockType::Dim; };
template <typename _RealType> struct DimensionTrait<_RealType, aol::Vector<_RealType> > { static const qc::Dimension Dim = qc::QC_1D; };

template <typename _RealType, typename _BlockSizeType> struct BlockSizeTrait;
template <typename _RealType> struct BlockSizeTrait<_RealType, int> {
  static bool isSquare ( const int /*BlockSize*/ ) { return true; }
  static int blockLength ( const int BlockSize ) { return BlockSize; }
};
template <typename _RealType> struct BlockSizeTrait<_RealType, qc::GridSize<qc::QC_2D> > {
  static bool isSquare ( const qc::GridSize<qc::QC_2D> &BlockSize ) { return ( BlockSize.getNumX ( ) == BlockSize.getNumY ( ) ); }
  static int blockLength ( const qc::GridSize<qc::QC_2D> &BlockSize ) {
    if ( !isSquare ( BlockSize ) ) throw aol::Exception ( "Block must be square!", __FILE__, __LINE__ );
    else return BlockSize.getNumX ( );
  }
};
template <typename _RealType> struct BlockSizeTrait<_RealType, qc::GridSize<qc::QC_3D> > {
  static bool isSquare ( const qc::GridSize<qc::QC_3D> &BlockSize ) { return ( BlockSize.getNumY ( ) == BlockSize.getNumZ ( ) ); }
  static int blockLength ( const qc::GridSize<qc::QC_3D> &BlockSize ) {
    if ( !isSquare ( BlockSize ) ) throw aol::Exception ( "Block must be square!", __FILE__, __LINE__ );
    else return BlockSize.getNumY ( );
  }
};


template <typename _RealType, typename _BlockType, typename _BlockSizeType>
class BlockTransformOp {
  typedef _RealType RealType;
  typedef _BlockType BlockType;
  typedef _BlockSizeType BlockSizeType;
  typedef void ( BlockTransformOp<RealType, BlockType, BlockSizeType>::*BlockTransformType )( const BlockType&, BlockType&, enum qc::FourierTransformDirection );
  static const qc::Dimension Dim = DimensionTrait<_RealType, BlockType>::Dim;
protected:
  std::string _waveletName;
  BlockTransformType _blockTransform;
  
  aol::DeleteFlagPointer<qc::FastCosineTransform<Dim, RealType> > _fct, _fctInv;
  aol::DeleteFlagPointer<MatrixFormTransformation<RealType> > _mFct;
public:
  BlockTransformOp ( ) { setTransform ( BlockSizeType ( static_cast<short> ( 1 ) ), "dct" ); }
  
  BlockTransformOp ( const BlockSizeType &BlockSize, const std::string &Transform ) { setTransform ( BlockSize, Transform ); }
  
  void apply ( const BlockType &Arg, BlockType &Dest ) {
   (this->*this->_blockTransform) ( Arg, Dest, qc::FTForward );
  }
  
  void applyInverse ( const BlockType &Arg, BlockType &Dest ) {
    (this->*this->_blockTransform) ( Arg, Dest, qc::FTBackward );
  }
  
  void setTransform ( const BlockSizeType &BlockSize, const std::string &Transform, bool EnforceDirectTransform = false ) {
    if ( !EnforceDirectTransform && BlockSizeTrait<RealType, BlockSizeType>::isSquare ( BlockSize ) && BlockSizeTrait<RealType, BlockSizeType>::blockLength ( BlockSize ) <= 32 ) {
      _mFct.reset ( new MatrixFormTransformation<RealType> ( BlockSizeTrait<RealType, BlockSizeType>::blockLength ( BlockSize ), Transform ), true );
      _blockTransform = &BlockTransformOp::matrixFormTransform;
    } else {
      if ( Transform == "dct" ) {
        _fct.reset ( new qc::FastCosineTransform<Dim, RealType> ( BlockSize, qc::FTForward ), true );
        _fctInv.reset ( new qc::FastCosineTransform<Dim, RealType> ( BlockSize, qc::FTBackward ), true );
        _blockTransform = &BlockTransformOp::cosineTransform;
      } else if ( qc::DWT::isDWT ( Transform ) ) {
        _waveletName = Transform;
        _blockTransform = &BlockTransformOp::waveletTransform;
      } else throw aol::Exception ( "Specified transformation not avalable!", __FILE__, __LINE__ );
    }
  }
protected:
  void cosineTransform ( const BlockType &Arg, BlockType &Dest, enum qc::FourierTransformDirection Direction ) {
    if ( Direction == qc::FTForward ) _fct->apply ( Arg, Dest, true );
    else _fctInv->apply ( Arg, Dest, true );
  }
  
  void waveletTransform ( const BlockType &Arg, BlockType &Dest, enum qc::FourierTransformDirection Direction ) {
    qc::WaveletTransform<RealType> ( Arg, Dest, _waveletName, Direction );
  }
  
  void matrixFormTransform ( const BlockType &Arg, BlockType &Dest, enum qc::FourierTransformDirection Direction ) {
    _mFct->apply ( Arg, Dest, Direction );
  }
};

template <typename _RealType>
class BlockTransformOp<_RealType, qc::ScalarArray<_RealType, qc::QC_3D>, qc::GridSize<qc::QC_3D> > {
  typedef _RealType RealType;
  typedef qc::ScalarArray<_RealType, qc::QC_3D> BlockType;
  typedef qc::GridSize<qc::QC_3D> BlockSizeType;
  typedef void ( BlockTransformOp<RealType, BlockType, BlockSizeType>::*BlockTransformType )( const aol::Vector<RealType>&, aol::Vector<RealType>&, enum qc::FourierTransformDirection );
protected:
  qc::GridSize<qc::QC_3D> _blockSize;
  
  std::string _waveletName;
  BlockTransformType _blockTransform;
  
  BlockTransformOp<_RealType, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D> > _blockTransformOp2D;
  
  aol::DeleteFlagPointer<qc::FastCosineTransform<qc::QC_1D, RealType> > _fct, _fctInv;
public:
  BlockTransformOp ( ) { setTransform ( BlockSizeType ( static_cast<short> ( 1 ) ), "dct" ); }
  
  BlockTransformOp ( const BlockSizeType &BlockSize, const std::string &Transform )
    : _blockSize ( BlockSize ), _blockTransformOp2D ( qc::GridSize<qc::QC_2D> ( BlockSize[1], BlockSize[2] ), Transform ) { setTransform ( BlockSize, Transform ); }
  
  void apply ( const BlockType &Arg, BlockType &Dest ) {
    qc::ScalarArray<RealType, qc::QC_2D> arg2D ( _blockSize[1], _blockSize[2] ), dest2D ( arg2D );
    for ( short x=0; x<_blockSize[0] ; ++x ) {
      Arg.getSlice ( qc::QC_X, x, arg2D );
      _blockTransformOp2D.apply ( arg2D, dest2D );
      Dest.putSlice ( qc::QC_X, x, dest2D );
    }
    
    aol::Vector<RealType> arg1D ( _blockSize[0] ), dest1D ( arg1D );
    for ( short yBlock=0; yBlock<_blockSize[2] ; ++yBlock ) {
      for ( short xBlock=0; xBlock<_blockSize[1] ; ++xBlock ) {
        for ( short x=0; x<_blockSize[0] ; ++x ) arg1D[x] = Dest.get ( x, xBlock, yBlock );
        (this->*this->_blockTransform) ( arg1D, dest1D, qc::FTForward );
        for ( short x=0; x<_blockSize[0] ; ++x ) Dest.set ( x, xBlock, yBlock, dest1D[x] );
      }
    }
  }
  
  void applyInverse ( const BlockType &Arg, BlockType &Dest ) {
    qc::ScalarArray<RealType, qc::QC_2D> arg2D ( _blockSize[1], _blockSize[2] ), dest2D ( arg2D );
    for ( short x=0; x<_blockSize[0] ; ++x ) {
      Arg.getSlice ( qc::QC_X, x, arg2D );
      _blockTransformOp2D.applyInverse ( arg2D, dest2D );
      Dest.putSlice ( qc::QC_X, x, dest2D );
    }
    
    aol::Vector<RealType> arg1D ( _blockSize[0] ), dest1D ( arg1D );
    for ( short yBlock=0; yBlock<_blockSize[2] ; ++yBlock ) {
      for ( short xBlock=0; xBlock<_blockSize[1] ; ++xBlock ) {
        for ( short x=0; x<_blockSize[0] ; ++x ) arg1D[x] = Dest.get ( x, xBlock, yBlock );
        (this->*this->_blockTransform) ( arg1D, dest1D, qc::FTBackward );
        for ( short x=0; x<_blockSize[0] ; ++x ) Dest.set ( x, xBlock, yBlock, dest1D[x] );
      }
    }
  }
  
  void setTransform ( const BlockSizeType &BlockSize, const std::string &Transform, bool /*EnforceDirectTransform*/ = false ) {
    _blockTransformOp2D.setTransform ( qc::GridSize<qc::QC_2D> ( BlockSize[1], BlockSize[2] ), Transform );
    _blockSize = BlockSize;
    
    if ( Transform == "dct" ) {
      _fct.reset ( new qc::FastCosineTransform<qc::QC_1D, RealType> ( BlockSize[0], qc::FTForward ), true );
      _fctInv.reset ( new qc::FastCosineTransform<qc::QC_1D, RealType> ( BlockSize[0], qc::FTBackward ), true );
      _blockTransform = &BlockTransformOp::cosineTransform;
    } else if ( qc::DWT::isDWT ( Transform ) ) {
      _waveletName = Transform;
      _blockTransform = &BlockTransformOp::waveletTransform;
    } else throw aol::Exception ( "Specified transformation not avalable!", __FILE__, __LINE__ );
  }
protected:
  void cosineTransform ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest, enum qc::FourierTransformDirection Direction ) {
    if ( Direction == qc::FTForward ) _fct->apply ( Arg, Dest, true );
    else _fctInv->apply ( Arg, Dest, true );
  }
  
  void waveletTransform ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest, enum qc::FourierTransformDirection Direction ) {
    qc::WaveletTransform<RealType> ( Arg, Dest, _waveletName, Direction );
  }
};



/*
 * Block stack transform operators
 */
template <typename _RealType, typename _BlockType, typename _BockSizeType>
class BlockStackTransformOp {
  typedef _RealType RealType;
  typedef _BlockType BlockType;
  typedef _BockSizeType BlockSizeType;
  typedef BlockStack<RealType, BlockType, BlockSizeType> BlockStackType;
protected:
  BlockTransformOp<RealType, BlockType, BlockSizeType> _blockTransformOp;
  aol::RandomAccessContainer<aol::DeleteFlagPointer<BlockTransformOp<RealType, aol::Vector<RealType>, int> > > _stackTransformOps;
public:
  BlockStackTransformOp ( ) : _blockTransformOp ( BlockSizeType ( static_cast<short> ( 1 ) ), "dct" ), _stackTransformOps ( ) { setStackTransforms ( 1, "haar" ); }
  
  BlockStackTransformOp ( const BlockSizeType &BlockSize, const std::string &BlockTransform, const int MaxNumBlocks, const std::string &StackTransform )
    : _blockTransformOp ( BlockSize, BlockTransform ), _stackTransformOps ( ) { setStackTransforms ( MaxNumBlocks, StackTransform ); }
  
  void apply ( const short NumBlocks, const BlockStackType &Arg, BlockStackType &Dest ) {
    for ( short z=0; z<NumBlocks ; ++z ) _blockTransformOp.apply ( Arg[z], Dest[z] );
    
    aol::Vector<RealType> tmpArg ( NumBlocks ), tmpDest ( NumBlocks );
    const int stackTransformOpIdx = floor ( log2 ( NumBlocks ) );
    for ( int k=0; k<Arg[0].size ( ) ; ++k ) {
      for ( short z=0; z<NumBlocks ; ++z ) tmpArg[z] = Dest.get ( z, k );
      _stackTransformOps[stackTransformOpIdx]->apply ( tmpArg, tmpDest );
      for ( short z=0; z<NumBlocks ; ++z ) Dest.set ( z, k, tmpDest[z] );
    }
  }
  
  void applyInverse ( const short NumBlocks, const BlockStackType &Arg, BlockStackType &Dest ) {
    for ( short z=0; z<NumBlocks ; ++z ) _blockTransformOp.applyInverse ( Arg[z], Dest[z] );
    
    aol::Vector<RealType> tmpArg ( NumBlocks ), tmpDest ( NumBlocks );
    const int stackTransformOpIdx = floor ( log2 ( NumBlocks ) );
    for ( int k=0; k<Arg[0].size ( ) ; ++k ) {
      for ( short z=0; z<NumBlocks ; ++z ) tmpArg[z] = Dest.get ( z, k );
      _stackTransformOps[stackTransformOpIdx]->applyInverse ( tmpArg, tmpDest );
      for ( short z=0; z<NumBlocks ; ++z ) Dest.set ( z, k, tmpDest[z] );
    }
  }
  
  void applyStackTransform ( const short NumBlocks, const BlockStackType &Arg, BlockStackType &Dest ) {
    aol::Vector<RealType> tmpArg ( NumBlocks ), tmpDest ( NumBlocks );
    const int stackTransformOpIdx = floor ( log2 ( NumBlocks ) );
    for ( int k=0; k<Arg[0].size ( ) ; ++k ) {
      for ( short z=0; z<NumBlocks ; ++z ) tmpArg[z] = Arg.get ( z, k );
      _stackTransformOps[stackTransformOpIdx]->apply ( tmpArg, tmpDest );
      for ( short z=0; z<NumBlocks ; ++z ) Dest.set ( z, k, tmpDest[z] );
    }
  }
  
  void setBlockTransform ( const BlockSizeType &BlockSize, const std::string &Transform ) {
    _blockTransformOp.setTransform ( BlockSize, Transform );
  }
  
  void setStackTransforms ( const int MaxNumBlocks, const std::string &Transform ) {
    _stackTransformOps.reallocate ( floor ( log2 ( MaxNumBlocks ) ) + 1 );
    for ( int z=0; z<=log2 ( MaxNumBlocks ) ; ++z )
      _stackTransformOps[z].reset ( new BlockTransformOp<RealType, aol::Vector<RealType>, int> ( aol::Pow ( 2, z ), Transform ), true );
  }
  
  void setTransforms ( const BlockSizeType &BlockSize, const std::string &BlockTransform, const int MaxNumBlocks, const std::string &StackTransform ) {
    setBlockTransform ( BlockSize, BlockTransform );
    setStackTransforms ( MaxNumBlocks, StackTransform );
  }
};



/*
 * Block stack filtering operators
 */
template <typename _RealType, typename _BlockType, typename _BlockSizeType>
class HTBlockStackFilteringOp {
  typedef _RealType RealType;
  typedef BlockStack<_RealType, _BlockType, _BlockSizeType> BlockStackType;
protected:
  RealType _threshold;
public:
  HTBlockStackFilteringOp ( ) : _threshold ( 0.0 ) { }
  
  HTBlockStackFilteringOp ( const RealType Threshold ) : _threshold ( Threshold ) { }
  
  void apply ( const short NumBlocks, const BlockStackType &Arg, BlockStackType &Dest, RealType &Weight ) const {
    Weight = 0;
    for ( int k=0; k<NumBlocks ; ++k ) {
      for ( int l=0; l<Arg[k].size ( ) ; ++l ) {
        if ( aol::Abs<RealType> ( Arg[k][l] ) > _threshold ) {
          Dest[k][l] = Arg[k][l];
          Weight += 1;
        } else {
          Dest[k][l] = 0.0;
        }
      }
    }
    Weight = ( Weight >= 1 ) ? 1.0 / Weight : 1.0;
  }
  
  void setThreshold ( const RealType Threshold ) {
    _threshold = Threshold;
  }
};


template <typename _RealType, typename _BlockType, typename _BlockSizeType>
class WienerBlockStackFilteringOp {
  typedef _RealType RealType;
  typedef BlockStack<_RealType, _BlockType, _BlockSizeType> BlockStackType;
protected:
  RealType _stdDev, _stdDevSqr;
public:
  WienerBlockStackFilteringOp ( ) : _stdDev ( 0.0 ), _stdDevSqr ( 0.0 ) { }
  
  WienerBlockStackFilteringOp ( const RealType NoiseStdDev ) : _stdDev ( NoiseStdDev ), _stdDevSqr ( aol::Sqr<RealType> ( NoiseStdDev ) ) { }
  
  void apply ( const short NumBlocks, const BlockStackType &Initial, const BlockStackType &Estimate, BlockStackType &Dest, RealType &Weight ) const {
    Weight = 0;
    RealType weight, estimateBlockEntrySqr;
    for ( int k=0; k<NumBlocks ; ++k ) {
      for ( int l=0; l<Initial[k].size ( ) ; ++l ) {
        estimateBlockEntrySqr = aol::Sqr<RealType> ( Estimate[k][l] );
        weight = estimateBlockEntrySqr / ( estimateBlockEntrySqr + _stdDevSqr );
        Dest[k][l] = weight * Initial[k][l];
        Weight += aol::Sqr<RealType> ( weight );
      }
    }
    Weight = 1.0 / Weight;
  }
  
  void setNoiseStdDev ( const RealType NoiseStdDev ) {
    _stdDev = NoiseStdDev;
    _stdDevSqr = aol::Sqr<RealType> ( NoiseStdDev );
  }
};



/*
 * Block stack denoising operators
 */
template <typename _RealType, typename _BlockType, typename _BlockSizeType>
class HTBlockStackDenoisingOp {
  typedef _RealType RealType;
  typedef _BlockType BlockType;
  typedef _BlockSizeType BlockSizeType;
  typedef BlockStack<RealType, BlockType, BlockSizeType> BlockStackType;
protected:
  BlockStackTransformOp<RealType, BlockType, BlockSizeType> _blockStackTransformOp;
  HTBlockStackFilteringOp<RealType, BlockType, BlockSizeType> _blockStackHTOp;
  BlockStackType _tmpBlockStack;
public:
  HTBlockStackDenoisingOp ( ) : _blockStackTransformOp ( ), _blockStackHTOp ( ) { }
  
  HTBlockStackDenoisingOp ( const BlockSizeType &BlockSize, const std::string &BlockTransform, const short MaxNumBlocks, const std::string &StackTransform,
                            const _RealType Threshold )
    : _blockStackTransformOp ( BlockSize, BlockTransform, MaxNumBlocks, StackTransform ), _blockStackHTOp ( Threshold ) {
    reallocate ( BlockSize, MaxNumBlocks );
  }
  
  void apply ( const short NumBlocks, const BlockStackType &Arg, BlockStackType &Dest, RealType &Weight ) {
    _blockStackTransformOp.apply ( NumBlocks, Arg, Dest );
    _blockStackHTOp.apply ( NumBlocks, Dest, _tmpBlockStack, Weight );
    _blockStackTransformOp.applyInverse ( NumBlocks, _tmpBlockStack, Dest );
  }
  
  void applyUsingOnlyStackTransform ( const short NumBlocks, const BlockStackType &Arg, BlockStackType &Dest, RealType &Weight ) {
    _blockStackTransformOp.applyStackTransform ( NumBlocks, Arg, Dest );
    _blockStackHTOp.apply ( NumBlocks, Dest, _tmpBlockStack, Weight );
    _blockStackTransformOp.applyInverse ( NumBlocks, _tmpBlockStack, Dest );
  }
  
  void setTransforms ( const BlockSizeType &BlockSize, const std::string &BlockTransform, const short MaxNumBlocks, const std::string &StackTransform ) {
    reallocate ( BlockSize, MaxNumBlocks );
    _blockStackTransformOp.setTransforms ( BlockSize, BlockTransform, MaxNumBlocks, StackTransform );
  }
  
  void setThreshold ( const RealType Threshold ) {
    _blockStackHTOp.setThreshold ( Threshold );
  }
protected:
  void reallocate ( const BlockSizeType &BlockSize, const short MaxNumBlocks ) {
    _tmpBlockStack.reallocate ( MaxNumBlocks, BlockSize );
  }
};


template <typename _RealType, typename _BlockType,  typename _BlockSizeType>
class WienerBlockStackDenoisingOp {
  typedef _RealType RealType;
  typedef _BlockType BlockType;
  typedef _BlockSizeType BlockSizeType;
  typedef BlockStack<RealType, BlockType, BlockSizeType> BlockStackType;
protected:
  BlockStackTransformOp<RealType, BlockType, BlockSizeType> _blockStackTransformOp;
  WienerBlockStackFilteringOp<RealType, BlockType, BlockSizeType> _blockStackWienerOp;
  BlockStackType _filteredBlockStack, _initialTransformedBlockStack, _estimateTransformedBlockStack;
public:
  WienerBlockStackDenoisingOp ( ) : _blockStackTransformOp ( ), _blockStackWienerOp ( ) { }
  
  WienerBlockStackDenoisingOp ( const BlockSizeType& BlockSize, const std::string &BlockTransform, const short MaxNumBlocks, const std::string &StackTransform,
                                const _RealType NoiseStdDev )
  : _blockStackTransformOp ( BlockSize, BlockTransform, MaxNumBlocks, StackTransform ), _blockStackWienerOp ( NoiseStdDev ) {
    reallocate ( BlockSize, MaxNumBlocks );
  }
  
  void apply ( const short NumBlocks, const BlockStackType &Initial, const BlockStackType &Estimate, BlockStackType &Dest, RealType &Weight ) {
    _blockStackTransformOp.apply ( NumBlocks, Initial, _initialTransformedBlockStack );
    _blockStackTransformOp.apply ( NumBlocks, Estimate, _estimateTransformedBlockStack );
    _blockStackWienerOp.apply ( NumBlocks, _initialTransformedBlockStack, _estimateTransformedBlockStack, _filteredBlockStack, Weight );
    _blockStackTransformOp.applyInverse ( NumBlocks, _filteredBlockStack, Dest );
  }
  
  void applyUsingOnlyStackTransform ( const short NumBlocks, const BlockStackType &Initial, const BlockStackType &Estimate, BlockStackType &Dest, RealType &Weight ) {
    _blockStackTransformOp.applyStackTransform ( NumBlocks, Initial, _initialTransformedBlockStack );
    _blockStackTransformOp.applyStackTransform ( NumBlocks, Estimate, _estimateTransformedBlockStack );
    _blockStackWienerOp.apply ( NumBlocks, _initialTransformedBlockStack, _estimateTransformedBlockStack, _filteredBlockStack, Weight );
    _blockStackTransformOp.applyInverse ( NumBlocks, _filteredBlockStack, Dest );
  }
  
  void setTransforms ( const BlockSizeType &BlockSize, const std::string &BlockTransform, const short MaxNumBlocks, const std::string &StackTransform ) {
    reallocate ( BlockSize, MaxNumBlocks );
    this->_blockStackTransformOp.setTransforms ( BlockSize, BlockTransform, MaxNumBlocks, StackTransform );
  }
  
  void setNoiseStdDev ( const RealType NoiseStdDev ) {
    _blockStackWienerOp.setNoiseStdDev ( NoiseStdDev );
  }
protected:
  void reallocate ( const BlockSizeType &BlockSize, const short MaxNumBlocks ) {
    _initialTransformedBlockStack.reallocate ( MaxNumBlocks, BlockSize );
    _estimateTransformedBlockStack.reallocate ( MaxNumBlocks, BlockSize );
    _filteredBlockStack.reallocate ( MaxNumBlocks, BlockSize );
  }
};

#endif