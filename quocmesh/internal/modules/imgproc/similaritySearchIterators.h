#ifndef __SIMILARITYSEARCHITERATORS_H
#define __SIMILARITYSEARCHITERATORS_H

#include <aol.h>
#include <vec.h>
#include <geom.h>

#include "patternAnalysis.h"
#include "featureSegmentation.h"


namespace im {

/**
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, qc::Dimension Dim>
class SimilaritySearchIterator {
protected:
  qc::CoordType _x0, _xEnd, _xL, _xR;
  qc::CoordType _cur;
public:
  SimilaritySearchIterator ( const qc::CoordType &X0, const qc::CoordType &XEnd )
    : _x0 ( X0 ), _xEnd ( XEnd ) { }
  
  virtual ~SimilaritySearchIterator ( ) { };
  
  virtual void update ( const _RealType Distance ) = 0;
  
  virtual bool notAtEnd ( ) const = 0;
  
  virtual SimilaritySearchIterator& operator++ ( ) = 0;
  
  qc::CoordType& operator* ( ) {
    return _cur;
  }
  
  virtual bool isCurFinal ( ) const {
    return true;
  }
  
  virtual int getNumNonLocalWindows ( ) const {
    return 0;
  }
};

/**
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, qc::Dimension Dim>
class LocalSimilaritySearchIterator : public SimilaritySearchIterator<_RealType, Dim> {
protected:
  qc::GridSize<qc::QC_2D> _searchWindowSize;
public:
  LocalSimilaritySearchIterator ( const qc::CoordType &XRef,
                                  const qc::CoordType &X0, const qc::CoordType &XEnd,
                                  const qc::GridSize<qc::QC_2D> &SearchWindowSize )
    : SimilaritySearchIterator<_RealType, Dim> ( X0, XEnd ),
      _searchWindowSize ( SearchWindowSize ) {
    aol::Vec2<short> searchWindowOffset ( ( _searchWindowSize.getNumX ( ) - 1 ) / 2, ( _searchWindowSize.getNumY ( ) - 1 ) / 2 );
    for ( int d=0; d<2 ; ++d ) {
      this->_xL[d] = aol::Max<int> ( XRef[d] - searchWindowOffset[d], this->_x0[d] );
      this->_xR[d] = aol::Min<int> ( XRef[d] + searchWindowOffset[d], this->_xEnd[d] - 1 );
    }
    if ( Dim == qc::QC_3D ) {
      this->_xL[2] = this->_x0[2];
      this->_xR[2] = this->_xEnd[2]-1;
    }
    this->_cur = this->_xL;
  }
  
  void update ( const _RealType /*Distance*/ ) { }
  
  bool notAtEnd ( ) const {
    return this->_cur[Dim-1] <= this->_xR[Dim-1];
  }
  
  LocalSimilaritySearchIterator& operator++ ( ) {
    ++this->_cur[0];
    for ( int d=1; d<Dim ; ++d ) {
      if ( this->_cur[d-1] <= this->_xR[d-1] ) break;
      this->_cur[d-1] = this->_xL[d-1];
      ++this->_cur[d];
    }
    
    return *this;
  }
};

/**
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, qc::Dimension Dim>
class GlobalSimilaritySearchIterator : public LocalSimilaritySearchIterator<_RealType, Dim> {
public:
  GlobalSimilaritySearchIterator ( const qc::CoordType &X0, const qc::CoordType &XEnd )
   : LocalSimilaritySearchIterator<_RealType, Dim> ( X0, X0, XEnd, qc::GridSize<qc::QC_2D> ( static_cast<short> ( 2 * XEnd.getMaxValue ( ) + 1 ) ) ) { }
};


/**
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, qc::Dimension Dim>
class BasePiecewisePeriodicSimilaritySearchIterator : public SimilaritySearchIterator<_RealType, Dim> {
  typedef _RealType RealType;
protected:
  int _localNs;
  int _localNsOffset;
  bool _curFinal, _notAtEnd;
  
  std::list<qc::CoordType> _latticeNodes;
  std::list<std::pair<qc::CoordType, short> > _latticeLeafs;
  const int _leafMaxAge;
  qc::BitArray<Dim> _visited;
  qc::CoordType _localMinCoords;
  RealType _localMinDist;
  qc::CoordType _localXL, _localXR;
  
  int _numNonLocalWindows;
  const int _minNumNonLocalWindows, _maxNumNonLocalWindows;
  
  aol::RandomAccessContainer<aol::MultiVector<_RealType> > _latticeVectors;
  const qc::ScalarArray<int, Dim> &_hardSegmentation;
  int _curSegment;
public:
  BasePiecewisePeriodicSimilaritySearchIterator ( const qc::CoordType &XRef,
                                                  const qc::CoordType &X0, const qc::CoordType &XEnd, const int LocalNs,
                                                  const aol::RandomAccessContainer<aol::MultiVector<_RealType> > &LatticeVectors,
                                                  const qc::ScalarArray<int, Dim> &HardSegmentation,
                                                  const int MinNumNonLocalWindows = 0, const int MaxNumNonLocalWindows = 0 )
    : SimilaritySearchIterator<RealType, Dim> ( X0, XEnd ),
      _localNs ( LocalNs ), _localNsOffset ( ( _localNs - 1 ) / 2 ),
      _notAtEnd ( true ),
      _leafMaxAge ( 3 ),
      _visited ( HardSegmentation.getSize ( ) ),
      _minNumNonLocalWindows ( MinNumNonLocalWindows ), _maxNumNonLocalWindows ( MaxNumNonLocalWindows ),
      _latticeVectors ( LatticeVectors ), _hardSegmentation ( HardSegmentation ),
      _curSegment ( _hardSegmentation.get ( XRef ) ) {
    this->_xL = this->_x0;
    for ( int d=0; d<Dim ; ++d ) this->_xR[d] = this->_xEnd[d] - 1;
    _numNonLocalWindows = 0;
    setNodeVisited ( XRef );
    localReset ( XRef );
  }
  
  void update ( const RealType Distance ) {
    if ( Distance < _localMinDist ) {
      _localMinDist = Distance;
      _localMinCoords = this->_cur;
    }
  }
  
  bool notAtEnd ( ) const {
    return _notAtEnd;
  }
  
  BasePiecewisePeriodicSimilaritySearchIterator& operator++ ( ) {
    if ( !_curFinal ) {
      // Local similarity search
      if ( this->_cur[0] < _localXR[0] )
        ++this->_cur[0];
      else if ( this->_cur[1] < _localXR[1] ) {
        this->_cur[0] = _localXL[0];
        ++this->_cur[1];
      } else {
        this->_cur = _localMinCoords;
        _curFinal = true;
        ++_numNonLocalWindows;
        if ( _maxNumNonLocalWindows > 0 && _numNonLocalWindows >= _maxNumNonLocalWindows )
          _notAtEnd = false;
        addNeighborsToList ( this->_cur );
      }
    } else {
      // If there are no more nodes and the minimum number of nodes has not been reached
      // try to find additional nodes by advancing along leafs
      if ( _latticeNodes.size ( ) == 0 && _numNonLocalWindows < _minNumNonLocalWindows ) {
        while ( _latticeNodes.size ( ) == 0 && _latticeLeafs.size ( ) > 0 ) {
          std::pair<qc::CoordType, short> latticeLeaf = _latticeLeafs.front ( );
          _latticeLeafs.pop_front ( );
          addNeighborsToList ( latticeLeaf.first, latticeLeaf.second );
        }
      }
      // If there are any nodes at this point, step to the first node in the list
      if ( _latticeNodes.size ( ) > 0 ) {
        this->_cur = _latticeNodes.front ( );
        _latticeNodes.pop_front ( );
        localReset ( this->_cur );
      } else {
        _notAtEnd = false;
      }
    }
    
    return *this;
  }
  
  virtual bool isCurFinal ( ) const {
    return _curFinal;
  }
  
  int getNumNonLocalWindows ( ) const {
    return _numNonLocalWindows;
  }
protected:
  void addNodeToList ( const qc::CoordType &Node, const short ParentAge = 0 ) {
    // Check whether node was not visited yet and is inside image
    if ( aol::InsideQuad<Dim> ( this->_x0, this->_xEnd, Node ) && !_visited.get ( Node ) ) {
      // Check whether node belongs to current segment
      if ( _hardSegmentation.get ( Node ) == _curSegment ) {
        _latticeNodes.push_back ( Node );
        setNodeVisited ( Node );
      } else {
        // Node is a leaf
        if ( ParentAge < _leafMaxAge ) {
          _latticeLeafs.push_back ( std::pair<qc::CoordType,short> ( Node, ParentAge + 1 ) );
          setNodeVisited ( Node );
        }
      }
    }
  }
  
  virtual void addNeighborsToList ( const qc::CoordType &ParentNode, const short ParentAge = 0 ) {
    qc::CoordType latticeNode;
    for ( int i=0; i<2 ; ++i ) {
      for ( int sign=-1; sign<=1 ; sign+=2 ) {
        latticeNode = ParentNode;
        for ( int d=0; d<2 ; ++d )
          latticeNode[d] += round ( sign * _latticeVectors[_curSegment][i][d] );
        addNodeToList ( latticeNode, ParentAge );
      }
    }
  }
  
  void setNodeVisited ( const qc::CoordType &Node ) {
    qc::CoordType xLocal ( Node );
    qc::CoordType localXL, localXR;
    for ( int d=0; d<2 ; ++d ) {
      localXL[d] = aol::Max<int> ( Node[d] - _localNsOffset, this->_x0[d] );
      localXR[d] = aol::Min<int> ( Node[d] + _localNsOffset, this->_xEnd[d] - 1 );
    }
    for ( short y=localXL[1] ; y<=localXR[1] ; ++y ) {
      for ( short x=localXL[0] ; x<=localXR[0] ; ++x ) {
        xLocal[0] = x;
        xLocal[1] = y;
        _visited.set ( xLocal, true );
      }
    }
  }
  
  void localReset ( const qc::CoordType &Center ) {
    _localMinCoords.set ( Center );
    _localMinDist = aol::NumberTrait<RealType>::getInf ( );
    for ( int d=0; d<2 ; ++d ) {
      _localXL[d] = aol::Max<int> ( Center[d] - _localNsOffset, this->_x0[d] );
      _localXR[d] = aol::Min<int> ( Center[d] + _localNsOffset, this->_xEnd[d] - 1 );
    }
    if ( Dim == qc::QC_3D ) {
      _localXL[2] = Center[2];
      _localXR[2] = Center[2];
    }
    this->_cur = _localXL;
    _curFinal = false;
  }
};

/**
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, qc::Dimension Dim>
class PiecewisePeriodicSimilaritySearchIterator : public BasePiecewisePeriodicSimilaritySearchIterator<_RealType, Dim> {
public:
  PiecewisePeriodicSimilaritySearchIterator ( const qc::CoordType &XRef,
                                              const qc::CoordType &X0, const qc::CoordType &XEnd, const int LocalNs,
                                              const aol::RandomAccessContainer<aol::MultiVector<_RealType> > &LatticeVectors,
                                              const qc::ScalarArray<int, Dim> &HardSegmentation,
                                              const int MinNumNonLocalWindows = 0, const int MaxNumNonLocalWindows = 0 )
    : BasePiecewisePeriodicSimilaritySearchIterator<_RealType, Dim> ( XRef, X0, XEnd, LocalNs, LatticeVectors, HardSegmentation, MinNumNonLocalWindows, MaxNumNonLocalWindows ) { }
};
  
/**
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType>
class PiecewisePeriodicSimilaritySearchIterator<_RealType, qc::QC_3D> : public BasePiecewisePeriodicSimilaritySearchIterator<_RealType, qc::QC_3D> {
  const static qc::Dimension Dim = qc::QC_3D;
protected:
  const aol::RandomAccessContainer<qc::MultiArray<_RealType, qc::QC_2D> > &_displacements, &_invDisplacements;
  const _RealType _h;
public:
  PiecewisePeriodicSimilaritySearchIterator ( const qc::CoordType &XRef,
                                              const qc::CoordType &X0, const qc::CoordType &XEnd, const int LocalNs,
                                              const aol::RandomAccessContainer<aol::MultiVector<_RealType> > &LatticeVectors,
                                              const qc::ScalarArray<int, Dim> &HardSegmentation,
                                              const aol::RandomAccessContainer<qc::MultiArray<_RealType, qc::QC_2D> > &Displacements,
                                              const aol::RandomAccessContainer<qc::MultiArray<_RealType, qc::QC_2D> > &InvDisplacements,
                                              const int MinNumNonLocalWindows = 0, const int MaxNumNonLocalWindows = 0 )
    : BasePiecewisePeriodicSimilaritySearchIterator<_RealType, Dim> ( XRef, X0, XEnd, LocalNs, LatticeVectors, HardSegmentation, MinNumNonLocalWindows, MaxNumNonLocalWindows ),
      _displacements ( Displacements ), _invDisplacements ( InvDisplacements ),
      _h ( 1.0 / static_cast<_RealType> ( aol::Max<int> ( Displacements[0].getNumX ( ), Displacements[0].getNumY ( ) ) ) ) { }
  
  virtual void addNeighborsToList ( const qc::CoordType &ParentNode, const short ParentAge = 0 ) {
    // Add neighbors in z-direction
    qc::CoordType latticeNode ( ParentNode );
    transformNode ( latticeNode, _displacements );
    ++latticeNode[2];
    if ( aol::InsideQuad<Dim> ( this->_x0, this->_xEnd, latticeNode ) ) {
      transformNode ( latticeNode, _invDisplacements );
      this->addNodeToList ( latticeNode, ParentAge );
    }
    
    // Add neighbors in x-y plane
    for ( int i=0; i<2 ; ++i ) {
      for ( int sign=-1; sign<=1 ; sign+=2 ) {
        latticeNode = ParentNode;
        transformNode ( latticeNode, _displacements );
        for ( int d=0; d<2 ; ++d )
          latticeNode[d] += round ( sign * this->_latticeVectors[this->_curSegment][i][d] );
        if ( aol::InsideQuad<Dim> ( this->_x0, this->_xEnd, latticeNode ) ) {
          transformNode ( latticeNode, _invDisplacements );
          this->addNodeToList ( latticeNode, ParentAge );
        }
      }
    }

    latticeNode = ParentNode;
    transformNode ( latticeNode, _displacements );
    --latticeNode[2];
    if ( aol::InsideQuad<Dim> ( this->_x0, this->_xEnd, latticeNode ) ) {
      transformNode ( latticeNode, _invDisplacements );
      this->addNodeToList ( latticeNode, ParentAge );
    }
  }
protected:
  void transformNode ( qc::CoordType &Node, const aol::RandomAccessContainer<qc::MultiArray<_RealType, qc::QC_2D> > &Displacements ) {
    qc::CoordType node ( Node );
    for ( int d=0; d<2 ; ++d )
      Node[d] += aol::Rint ( Displacements[node[2]][d].get ( node[0], node[1] ) / _h );
  }
};
  
  
/**
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, qc::Dimension Dim>
class PixelSimilaritySearchIterator : public SimilaritySearchIterator<_RealType, Dim> {
  typedef _RealType RealType;
protected:
  aol::RandomAccessContainer<qc::CoordType> _pixels;
  int _curIdx;
public:
  PixelSimilaritySearchIterator ( const qc::CoordType &XRef,
                                  const qc::CoordType &X0, const qc::CoordType &XEnd,
                                  const qc::ScalarArray<int, Dim> &HardSegmentation )
    : SimilaritySearchIterator<_RealType, Dim> ( X0, XEnd ),
      _curIdx ( 0 ) {
    if ( HardSegmentation.get ( XRef ) < 0 ) {
      for ( qc::RectangularIterator<Dim> it ( X0, XEnd ); it.notAtEnd ( ) ; ++it ) {
        if ( HardSegmentation.get ( *it ) < 0 )
          _pixels.pushBack ( *it );
      }
      if ( notAtEnd ( ) ) this->_cur = _pixels[0];
    }
  }
  
  void update ( const RealType /*Distance*/ ) { }
  
  bool notAtEnd ( ) const {
    return ( _curIdx < _pixels.size ( ) );
  }
  
  PixelSimilaritySearchIterator& operator++ ( ) {
    ++_curIdx;
    if ( notAtEnd ( ) ) this->_cur = _pixels[_curIdx];
    
    return *this;
  }
  
  virtual bool isCurFinal ( ) const {
    return true;
  }
};
  
  
} // namespace im
  

#endif
