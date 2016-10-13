#ifndef __BLOCKSTRUCTURES_H
#define __BLOCKSTRUCTURES_H

#include <scalarArray.h>

template <typename _RealType, qc::Dimension Dim, typename _PictureType, typename _BlockType, typename _BlockSizeType>
class BaseConstBlockCollection {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
  typedef _BlockType BlockType;
  typedef _BlockSizeType BlockSizeType;
protected:
  const PictureType *_data;
  BlockSizeType _blockSize;
  int _blockNumPixels;
  aol::Vec2<short> _blockAnchor;
  BlockType _block;
public:
  BaseConstBlockCollection ( ) : _data ( NULL ), _blockSize ( static_cast<short> ( 0 ) ) { }
  
  virtual ~BaseConstBlockCollection ( ) { }
  
  void initialize ( const PictureType &Data, const BlockSizeType &BlockSize, const aol::Vec2<short> &BlockAnchor ) {
    if ( !blockAnchorInsideBlock ( BlockSize, BlockAnchor ) )
      throw aol::Exception ( "Block anchor outside block!", __FILE__, __LINE__ );
    
    _data = &Data;
    _blockSize = BlockSize;
    _blockNumPixels = getBlockNumPixels ( BlockSize );
    _blockAnchor = BlockAnchor;
    _block.reallocate ( BlockSize );
  }
  
  void reallocate ( const PictureType &/*Data*/, const BlockSizeType &/*BlockSize*/ ) {
    throw aol::Exception( "Reallocating is not supported by this class! Use NonConstBlockCollection instead!", __FILE__, __LINE__ );
  }
  
  RealType get ( qc::CoordType &/*XRef*/, const int /*K*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  virtual void getBlock ( BlockType &/*Block*/, const qc::CoordType &/*XRef*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  const BlockType& getBlockConstRef ( const qc::CoordType &XRef ) {
    getBlock ( this->_block, XRef );
    return this->_block;
  }
  
  BlockType& getBlockRef ( const qc::CoordType &/*XRef*/ ) {
    throw aol::Exception( "Non-const block references not available in this class! Use NonConstBlockCollection instead!", __FILE__, __LINE__ );
  }
  
  void setBlock ( const qc::CoordType &/*XRef*/, const BlockType &/*Block*/ ) {
    throw aol::Exception( "Setting blocks is not supported by this class! Use NonConstBlockCollection instead!", __FILE__, __LINE__ );
  }
  
  const PictureType& getData ( ) const {
    return (*_data);
  }
  
  short getNumX ( ) const {
    return _data->getNumX ( );
  }
  
  short getNumY ( ) const {
    return _data->getNumY ( );
  }
  
  short getNumZ ( ) const {
    return _data->getNumZ ( );
  }
  
  virtual short getNumXEff ( ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  virtual short getNumYEff ( ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  short getX0 ( ) const {
    return _blockAnchor[0];
  }
  
  short getY0 ( ) const {
    return _blockAnchor[1];
  }

  short getXEnd ( ) const {
    return getX0 ( ) + getNumXEff ( );
  }
  
  short getYEnd ( ) const {
    return getY0 ( ) + getNumYEff ( );
  }
  
  int numNodes ( ) const {
    return getNumX ( ) * getNumY ( ) * getNumZ ( );
  }
  
  int numInnerNodes ( ) const {
    return getNumXEff ( ) * getNumYEff ( ) * getNumZ ( );
  }
  
  const BlockSizeType& getBlockSize ( ) const {
    return _blockSize;
  }
  
  const aol::Vec2<short>& getBlockAnchor ( ) const {
    return _blockAnchor;
  }
  
  short getBlockAnchorX ( ) const {
	  return _blockAnchor[0];
  }
  
  short getBlockAnchorY ( ) const {
    return _blockAnchor[1];
  }
protected:
  virtual bool blockAnchorInsideBlock ( const BlockSizeType &/*BlockSize*/, const aol::Vec2<short> &/*BlockAnchor*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  virtual int getBlockNumPixels ( const BlockSizeType &/*BlockSize*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};


template <typename _RealType, qc::Dimension Dim, typename _PictureType, typename _BlockType, typename _BlockSizeType> class ConstBlockCollection;

template <typename _RealType>
class ConstBlockCollection<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_2D>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D> >
  : public BaseConstBlockCollection<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_2D>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D> > {
public:
  const static qc::Dimension Dim = qc::QC_2D;
  typedef qc::ScalarArray<_RealType, qc::QC_2D> BlockType;
  typedef qc::GridSize<qc::QC_2D> BlockSizeType;
protected:
  aol::Vector<short> _XLookup, _YLookup;
public:
  void initialize ( const qc::ScalarArray<_RealType, qc::QC_2D> &Data, const qc::GridSize<qc::QC_2D> &BlockSize, const aol::Vec2<short> &BlockAnchor ) {
    BaseConstBlockCollection<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_2D>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D> >::initialize ( Data, BlockSize, BlockAnchor );
    _XLookup.resize ( this->_blockNumPixels );
    _YLookup.resize ( this->_blockNumPixels );
    short k = 0;
    for ( short y=0; y<this->_blockSize.getNumY ( ) ; ++y ) {
      for ( short x=0; x<this->_blockSize.getNumX ( ) ; ++x ) {
        _XLookup[k] = x;
        _YLookup[k] = y;
        ++k;
      }
    }
  }
  
  _RealType get ( const qc::CoordType &XRef, const int K ) const {
    return this->_data->get ( XRef[0] + _XLookup[K] - this->_blockAnchor[0], XRef[1] + _YLookup[K] - this->_blockAnchor[1] );
  }
    
  void getBlock ( BlockType &Block, const qc::CoordType &XRef ) const {
    for ( short yBlock=0; yBlock<this->_blockSize.getNumY ( ) ; ++yBlock )
      for ( short xBlock=0; xBlock<this->_blockSize.getNumX ( ) ; ++xBlock )
        Block.set ( xBlock, yBlock, this->_data->get ( XRef[0] + xBlock - this->_blockAnchor[0], XRef[1] + yBlock - this->_blockAnchor[1] ) );
  }
    
  short getNumXEff ( ) const {
    return this->getNumX ( ) - this->_blockSize.getNumX ( );
  }
  
  short getNumYEff ( ) const {
    return this->getNumY ( ) - this->_blockSize.getNumY ( );
  }
protected:
  bool blockAnchorInsideBlock ( const BlockSizeType &BlockSize, const aol::Vec2<short> &BlockAnchor ) const {
    return ( BlockAnchor[0] >= 0 && BlockAnchor[0] < BlockSize.getNumX ( ) && BlockAnchor[1] >= 0 && BlockAnchor[1] < BlockSize.getNumY ( ) );
  }
  
  int getBlockNumPixels ( const BlockSizeType &BlockSize ) const {
    return BlockSize.getNumberOfNodes ( );
  }
};

template <typename _RealType>
class ConstBlockCollection<_RealType, qc::QC_3D, qc::ScalarArray<_RealType, qc::QC_3D>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D> >
  : public BaseConstBlockCollection<_RealType, qc::QC_3D, qc::ScalarArray<_RealType, qc::QC_3D>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D> > {
public:
  const static qc::Dimension Dim = qc::QC_3D;
    typedef qc::ScalarArray<_RealType, qc::QC_2D> BlockType;
    typedef qc::GridSize<qc::QC_2D> BlockSizeType;
protected:
  aol::Vector<short> _XLookup, _YLookup;
public:
  void initialize ( const qc::ScalarArray<_RealType, qc::QC_3D> &Data, const qc::GridSize<qc::QC_2D> &BlockSize, const aol::Vec2<short> &BlockAnchor ) {
    BaseConstBlockCollection<_RealType, qc::QC_3D, qc::ScalarArray<_RealType, qc::QC_3D>, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D> >::initialize ( Data, BlockSize, BlockAnchor );
    _XLookup.resize ( this->_blockNumPixels );
    _YLookup.resize ( this->_blockNumPixels );
    short k = 0;
    for ( short y=0; y<this->_blockSize.getNumY ( ) ; ++y ) {
      for ( short x=0; x<this->_blockSize.getNumX ( ) ; ++x ) {
        _XLookup[k] = x;
        _YLookup[k] = y;
        ++k;
      }
    }
  }
  
  _RealType get ( const qc::CoordType &XRef, const int K ) const {
    return this->_data->get ( XRef[0] + _XLookup[K] - this->_blockAnchor[0], XRef[1] + _YLookup[K] - this->_blockAnchor[1], XRef[2] );
  }
    
  void getBlock ( BlockType &Block, const qc::CoordType &XRef ) const {
    for ( short yBlock=0; yBlock<this->_blockSize.getNumY ( ) ; ++yBlock )
      for ( short xBlock=0; xBlock<this->_blockSize.getNumX ( ) ; ++xBlock )
        Block.set ( xBlock, yBlock, this->_data->get ( XRef[0] + xBlock - this->_blockAnchor[0], XRef[1] + yBlock - this->_blockAnchor[1], XRef[2] ) );
  }
  
  short getNumXEff ( ) const {
    return this->getNumX ( ) - this->_blockSize.getNumX ( );
  }
  
  short getNumYEff ( ) const {
    return this->getNumY ( ) - this->_blockSize.getNumY ( );
  }
protected:
  bool blockAnchorInsideBlock ( const BlockSizeType &BlockSize, const aol::Vec2<short> &BlockAnchor ) const {
    return ( BlockAnchor[0] >= 0 && BlockAnchor[0] < BlockSize.getNumX ( ) && BlockAnchor[1] >= 0 && BlockAnchor[1] < BlockSize.getNumY ( ) );
  }
  
  int getBlockNumPixels ( const BlockSizeType &BlockSize ) const {
    return BlockSize.getNumberOfNodes ( );
  }
};

template <typename _RealType>
class ConstBlockCollection<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_3D>, aol::Vector<_RealType>, int>
  : public BaseConstBlockCollection<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_3D>, aol::Vector<_RealType>, int> {
public:
  const static qc::Dimension Dim = qc::QC_2D;
  typedef aol::Vector<_RealType> BlockType;
  typedef int BlockSizeType;
public:
  _RealType get ( const qc::CoordType &XRef, const int K ) const {
    return this->_data->get ( K, XRef[0], XRef[1] );
  }
    
  void getBlock ( BlockType &Block, const qc::CoordType &XRef ) const {
    for ( int k=0; k<this->_blockSize ; ++k )
      Block[k] = this->_data->get ( k, XRef[0], XRef[1] );
  }
    
  short getNumXEff ( ) const {
    return this->getNumX ( );
  }
  
  short getNumYEff ( ) const {
    return this->getNumY ( );
  }
protected:
  bool blockAnchorInsideBlock ( const BlockSizeType &/*BlockSize*/, const aol::Vec2<short> &BlockAnchor ) const {
    return ( BlockAnchor[0] == 0 && BlockAnchor[1] == 0 );
  }
  
  int getBlockNumPixels ( const BlockSizeType &BlockSize ) const {
    return BlockSize;
  }
};

template <typename _RealType>
class ConstBlockCollection<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_3D>, qc::ScalarArray<_RealType, qc::QC_3D>, qc::GridSize<qc::QC_3D> >
: public BaseConstBlockCollection<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_3D>, qc::ScalarArray<_RealType, qc::QC_3D>, qc::GridSize<qc::QC_3D> > {
public:
  const static qc::Dimension Dim = qc::QC_2D;
  typedef qc::ScalarArray<_RealType, qc::QC_3D> BlockType;
  typedef qc::GridSize<qc::QC_3D> BlockSizeType;
protected:
  aol::Vector<short> _XLookup, _YLookup, _ZLookup;
public:
  void initialize ( const qc::ScalarArray<_RealType, qc::QC_3D> &Data, const qc::GridSize<qc::QC_3D> &BlockSize, const aol::Vec2<short> &BlockAnchor ) {
    BaseConstBlockCollection<_RealType, qc::QC_2D, qc::ScalarArray<_RealType, qc::QC_3D>, qc::ScalarArray<_RealType, qc::QC_3D>, qc::GridSize<qc::QC_3D> >::initialize ( Data, BlockSize, BlockAnchor );
    _XLookup.resize ( this->_blockNumPixels );
    _YLookup.resize ( this->_blockNumPixels );
    _ZLookup.resize ( this->_blockNumPixels );
    int k = 0;
    for ( short z=0; z<this->_blockSize.getNumZ ( ) ; ++z ) {
      for ( short y=0; y<this->_blockSize.getNumY ( ) ; ++y ) {
        for ( short x=0; x<this->_blockSize.getNumX ( ) ; ++x ) {
          _XLookup[k] = x;
          _YLookup[k] = y;
          _ZLookup[k] = z;
          ++k;
        }
      }
    }
  }
  
  _RealType get ( const qc::CoordType &XRef, const int K ) const {
    return this->_data->get ( _XLookup[K], XRef[0] + _YLookup[K] - this->_blockAnchor[0], XRef[1] + _ZLookup[K] - this->_blockAnchor[1] );
  }
  
  void getBlock ( BlockType &Block, const qc::CoordType &XRef ) const {
    for ( int yBlock=0; yBlock<this->_blockSize.getNumZ ( ) ; ++yBlock )
      for ( int xBlock=0; xBlock<this->_blockSize.getNumY ( ) ; ++xBlock )
        for ( int x=0; x<this->_blockSize.getNumX ( ) ; ++x )
          Block.set ( x, xBlock, yBlock, this->_data->get ( x, XRef[0] + xBlock, XRef[1] + yBlock ) );
  }
  
  short getNumXEff ( ) const {
    return this->getNumX ( ) - this->_blockSize.getNumY ( );
  }
  
  short getNumYEff ( ) const {
    return this->getNumY ( ) - this->_blockSize.getNumZ ( );
  }
protected:
  bool blockAnchorInsideBlock ( const BlockSizeType &BlockSize, const aol::Vec2<short> &BlockAnchor ) const {
    return ( BlockAnchor[0] >= 0 && BlockAnchor[0] < BlockSize.getNumY ( ) && BlockAnchor[1] >= 0 && BlockAnchor[1] < BlockSize.getNumZ ( ) );
  }
  
  int getBlockNumPixels ( const BlockSizeType &BlockSize ) const {
    return BlockSize.getNumberOfNodes ( );
  }
};



template <typename _RealType, qc::Dimension Dim, typename _PictureType, typename _BlockType, typename _BlockSizeType>
class BaseNonConstBlockCollection : public ConstBlockCollection<_RealType, Dim, _PictureType, _BlockType, _BlockSizeType> {
public:
  typedef _RealType RealType;
  typedef _PictureType PictureType;
  typedef _BlockType BlockType;
  typedef _BlockSizeType BlockSizeType;
protected:
  qc::AArray<BlockType, Dim> _blocks;
public:
  virtual void initialize ( const PictureType &Data, const BlockSizeType &BlockSize, const aol::Vec2<short> &BlockAnchor ) = 0;
  
  virtual void reallocate ( const PictureType &Data, const BlockSizeType &BlockSize ) = 0;
  
  RealType get ( const qc::CoordType &XRef, const int K ) const {
    return _blocks.getRef ( XRef )[K];
  }
  
  void getBlock ( BlockType &Block, const qc::CoordType &XRef ) const {
    Block = _blocks.getRef ( XRef );
  }
  
  const BlockType& getBlockConstRef ( const qc::CoordType &XRef ) {
    return _blocks.getRef ( XRef );
  }
  
  BlockType& getBlockRef ( const qc::CoordType &XRef ) {
    return _blocks.getRef ( XRef );
  }
  
  void setBlock ( const qc::CoordType &XRef, const BlockType& Block ) {
    _blocks.set ( XRef, Block );
  }
};

template <typename RealType, qc::Dimension Dim, typename PictureType, typename BlockType, typename BlockSizeType>
class NonConstBlockCollection;

template <typename RealType, typename PictureType, typename BlockType, typename BlockSizeType>
class NonConstBlockCollection<RealType, qc::QC_2D, PictureType, BlockType, BlockSizeType>
  : public BaseNonConstBlockCollection<RealType, qc::QC_2D, PictureType, BlockType, BlockSizeType> {
  const static qc::Dimension Dim = qc::QC_2D;
public:
  void initialize ( const PictureType &Data, const BlockSizeType &BlockSize, const aol::Vec2<short> &BlockAnchor ) {
    ConstBlockCollection<RealType, Dim, PictureType, BlockType, BlockSizeType>::initialize ( Data, BlockSize, BlockAnchor );
    this->_blocks.reallocate ( this->getNumX ( ), this->getNumY ( ) );
    for ( short y=this->getY0 ( ); y<this->getYEnd ( ) ; ++y ) {
      for ( short x=this->getX0 ( ); x<this->getXEnd ( ) ; ++x ) {
        this->_blocks.getRef ( x, y ) = BlockType ( this->_blockSize );
        for ( short k=0; k<this->_blockNumPixels ; ++k )
          this->_blocks.getRef ( x, y )[k] = ConstBlockCollection<RealType, Dim, PictureType, BlockType, BlockSizeType>::get ( qc::CoordType ( x, y, 0 ), k );
      }
    }
  }
    
  void reallocate ( const PictureType &Data, const BlockSizeType &BlockSize ) {
    this->_blocks.reallocate ( Data.getNumX ( ), Data.getNumY ( ) );
    for ( short y=this->getY0 ( ); y<this->getYEnd ( ) ; ++y )
      for ( short x=this->getX0 ( ); x<this->getXEnd ( ) ; ++x )
        this->_blocks.getRef ( x, y ).reallocate ( BlockSize );
  }
};

template <typename RealType, typename PictureType, typename BlockType, typename BlockSizeType>
class NonConstBlockCollection<RealType, qc::QC_3D, PictureType, BlockType, BlockSizeType>
  : public BaseNonConstBlockCollection<RealType, qc::QC_3D, PictureType, BlockType, BlockSizeType> {
  const static qc::Dimension Dim = qc::QC_3D;
public:
  void initialize ( const PictureType &Data, const BlockSizeType &BlockSize, const aol::Vec2<short> &BlockAnchor ) {
    ConstBlockCollection<RealType, Dim, PictureType, BlockType, BlockSizeType>::initialize ( Data, BlockSize, BlockAnchor );
    this->_blocks.reallocate ( this->getNumX ( ), this->getNumY ( ), this->getNumZ ( ) );
    for ( short z=0; z<this->getNumZ ( ) ; ++z ) {
      for ( short y=this->getY0 ( ); y<this->getYEnd ( ) ; ++y ) {
        for ( short x=this->getX0 ( ); x<this->getXEnd ( ) ; ++x ) {
          this->_blocks.getRef ( x, y, z ) = BlockType ( this->_blockSize );
          for ( short k=0; k<this->_blockNumPixels ; ++k )
            this->_blocks.getRef ( x, y, z )[k] = ConstBlockCollection<RealType, Dim, PictureType, BlockType, BlockSizeType>::get ( qc::CoordType ( x, y, z ), k );
        }
      }
    }
  }
  
  void reallocate ( const PictureType &Data, const BlockSizeType &BlockSize ) {
    this->_blocks.reallocate ( Data.getNumX ( ), Data.getNumY ( ), Data.getNumZ ( ) );
    for ( short z=0; z<this->getNumZ ( ) ; ++z )
      for ( short y=this->getY0 ( ); y<this->getYEnd ( ) ; ++y )
        for ( short x=this->getX0 ( ); x<this->getXEnd ( ) ; ++x )
          this->_blocks.getRef ( x, y, z ).reallocate ( BlockSize );
  }
};


template <typename _RealType, qc::Dimension Dim, typename _PictureType, typename _BlockType, typename _BlockSizeType>
class BlockCollection {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
  typedef _BlockType BlockType;
  typedef _BlockSizeType BlockSizeType;
protected:
  ConstBlockCollection<RealType, Dim, PictureType, BlockType, BlockSizeType> _constBlockCollection;
  NonConstBlockCollection<RealType, Dim, PictureType, BlockType, BlockSizeType> _nonConstBlockCollection;
  ConstBlockCollection<RealType, Dim, PictureType, BlockType, BlockSizeType> *_blockCollection;
  bool _isConst;
public:
  BlockCollection ( ) : _isConst ( true ) { }
  
  void initialize ( const PictureType &Data, const BlockSizeType &BlockSize, const aol::Vec2<short> &BlockAnchor,
                    bool ForceNonConst = false ) {
    if ( ForceNonConst || memorySufficientForNonConst ( Data, BlockSize ) ) {
      _nonConstBlockCollection.initialize ( Data, BlockSize, BlockAnchor );
      _blockCollection = &_nonConstBlockCollection;
      _isConst = false;
    } else {
      _constBlockCollection.initialize ( Data, BlockSize, BlockAnchor );
      _blockCollection = &_constBlockCollection;
      _isConst = true;
    }
  }
  
  bool reallocate ( const BlockCollection &Other, bool ForceNonConst = false ) {
    if ( ForceNonConst || memorySufficientForNonConst ( Other.getData ( ), Other.getBlockSize ( ) ) ) {
      _nonConstBlockCollection.reallocate ( Other.getData ( ), Other.getBlockSize ( ) );
      _blockCollection = &_nonConstBlockCollection;
      _isConst = false;
      return true;
    } else {
      _isConst = true;
      return false;
    }
  }
  
  RealType get ( const qc::CoordType &XRef, const int K ) const {
    return _blockCollection->get ( XRef, K );
  }
  
  void getBlock ( BlockType &Block, const qc::CoordType &XRef ) const {
    _blockCollection->getBlock ( Block, XRef );
  }
  
  const BlockType& getBlockConstRef ( const qc::CoordType &XRef ) {
    return _blockCollection->getBlockConstRef ( XRef );
  }
  
  BlockType& getBlockRef ( const qc::CoordType &XRef ) {
    return _blockCollection->getBlockRef ( XRef );
  }
  
  void setBlock ( qc::CoordType &XRef, const BlockType& Block ) {
    _blockCollection->setBlock ( XRef, Block );
  }
  
  const PictureType& getData ( ) const {
    return _blockCollection->getData ( );
  }
  
  short getNumX ( ) const {
    return _blockCollection->getNumX ( );
  }
  
  short getNumY ( ) const {
    return _blockCollection->getNumY ( );
  }
  
  short getNumZ ( ) const {
    return _blockCollection->getNumZ ( );
  }
  
  short getNumXEff ( ) const {
    return _blockCollection->getNumXEff ( );
  }
  
  short getNumYEff ( ) const {
    return _blockCollection->getNumYEff ( );
  }
  
  short getX0 ( ) const {
    return _blockCollection->getX0 ( );
  }
  
  short getY0 ( ) const {
    return _blockCollection->getY0 ( );
  }
  
  short getXEnd ( ) const {
    return _blockCollection->getYEnd ( );
  }
  
  short getYEnd ( ) const {
    return _blockCollection->getYEnd ( );
  }
  
  int numNodes ( ) const {
    return _blockCollection->numNodes ( );
  }
  
  int numInnerNodes ( ) const {
    return _blockCollection->numInnerNodes ( );
  }
  
  const BlockSizeType& getBlockSize ( ) const {
    return _blockCollection->getBlockSize ( );
  }
  
  const aol::Vec2<short>& getBlockAnchor ( ) const {
    return _blockCollection->getBlockAnchor ( );
  }
  
  short getBlockAnchorX ( ) const {
	  return _blockCollection->getBlockAnchorX ( );
  }
  
  short getBlockAnchorY ( ) const {
    return _blockCollection->getBlockAnchorY ( );
  }
  
  bool isConst ( ) const {
    return _isConst;
  }
protected:
  bool memorySufficientForNonConst ( const PictureType &/*Data*/, const BlockSizeType &/*BlockSize*/ ) const {
    return false; // TODO
  }
};





template <typename _RealType, typename _BlockType, typename _BlockSizeType>
class BaseBlockStack : public aol::VectorContainer<_BlockType> {
  typedef _RealType RealType;
  typedef _BlockType BlockType;
  typedef _BlockSizeType BlockSizeType;
protected:
  BlockSizeType _blockSize;
public:
  BaseBlockStack ( ) : _blockSize ( static_cast<short> ( 0 ) ) { }
  
  void reallocate ( const short NumBlocks, const BlockSizeType &BlockSize ) {
    _blockSize = BlockSize;
    this->clear();
    this->_data.reserve ( NumBlocks );
    for ( short i = 0; i < NumBlocks; ++i )
      this->_data.push_back ( new BlockType ( BlockSize ) );
  }
  
  void reallocate ( const BaseBlockStack &Other ) {
    reallocate ( Other.size ( ), Other.getBlockSize ( ) );
  }
  
  RealType get ( const short Z, const int K ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void set ( const short Z, const int K ) {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void setBlock ( const short K, const BlockType &Block ) {
    for ( short i=0; i<Block.size ( ) ; ++i )
      (*this)[K][i] = Block[i];
  }
  
  short getNumBlocks ( ) const {
    return this->size ( );
  }
  
  const BlockSizeType& getBlockSize ( ) const {
    return _blockSize;
  }
};

template <typename _RealType, typename _BlockType, typename _BlockSizeType> class BlockStack;

template <typename _RealType>
class BlockStack<_RealType, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D> >
  : public BaseBlockStack<_RealType, qc::ScalarArray<_RealType, qc::QC_2D>, qc::GridSize<qc::QC_2D> > {
  typedef _RealType RealType;
public:
  BlockStack ( ) { }
  
  BlockStack ( const short NumBlocks, const qc::GridSize<qc::QC_2D> &BlockSize ) {
    this->reallocate ( NumBlocks, BlockSize );
  }
  
  BlockStack ( const short NumBlocks, const short NumX, const short NumY ) {
    this->reallocate ( NumBlocks, qc::GridSize<qc::QC_2D> ( NumX, NumY ) );
  }
  
  RealType get ( const short Z, const int K ) const {
    return (*this)[Z][K];
  }
  
  RealType get ( const short Z, const short X, const short Y ) const {
    return (*this)[Z].get ( X, Y );
  }
  
  void set ( const short Z, const int K, const RealType Val ) {
    (*this)[Z][K] = Val;
  }
  
  void set ( const short Z, const short X, const short Y, const RealType Val ) {
    (*this)[Z].set ( X, Y, Val );
  }
};

template <typename _RealType>
class BlockStack<_RealType, aol::Vector<_RealType>, int>
  : public BaseBlockStack<_RealType, aol::Vector<_RealType>, int> {
  typedef _RealType RealType;
public:
  BlockStack ( ) { }
  
  BlockStack ( const short NumBlocks, const int BlockSize ) {
    this->reallocate ( NumBlocks, BlockSize );
  }
  
  RealType get ( const short Z, const int K ) const {
    return (*this)[Z][K];
  }
  
  void set ( const short Z, const int K, const RealType Val ) {
    (*this)[Z][K] = Val;
  }
};

template <typename _RealType>
class BlockStack<_RealType, qc::ScalarArray<_RealType, qc::QC_3D>, qc::GridSize<qc::QC_3D> >
  : public BaseBlockStack<_RealType, qc::ScalarArray<_RealType, qc::QC_3D>, qc::GridSize<qc::QC_3D> > {
  typedef _RealType RealType;
public:
  BlockStack ( ) { }
  
  BlockStack ( const short NumBlocks, const qc::GridSize<qc::QC_3D> BlockSize ) {
    this->reallocate ( NumBlocks, BlockSize );
  }
  
  RealType get ( const short ZStack, const int KBlock ) const {
    return (*this)[ZStack][KBlock];
  }
  
  RealType get ( const short ZStack, const short XBlock, const short YBlock, const short ZBlock ) const {
    return (*this)[ZStack].get ( XBlock, YBlock, ZBlock );
  }
  
  void set ( const short ZStack, const int KBlock, const RealType Val ) {
    (*this)[ZStack][KBlock] = Val;
  }
  
  void set ( const short ZStack, const short XBlock, const short YBlock, const short ZBlock, const RealType Val ) {
    (*this)[ZStack].set ( XBlock, YBlock, ZBlock, Val );
  }
};

#endif
