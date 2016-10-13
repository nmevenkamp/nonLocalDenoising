#ifndef __REFERENCEBLOCKITERATORS_H
#define __REFERENCEBLOCKITERATORS_H

#include <aol.h>
#include <vec.h>
#include <multiVector.h>
#include <progressBar.h>
#include <vectorExtensions.h>


namespace im {

/**
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, qc::Dimension Dim>
class ReferenceBlockIterator {
protected:
  aol::ProgressBar<> *_progressBar;
  const bool _quietMode;
public:
  ReferenceBlockIterator ( aol::ProgressBar<> *ProgressBar = NULL, bool Quiet = false ) : _progressBar ( ProgressBar ), _quietMode ( Quiet ) { }
  
  virtual ~ReferenceBlockIterator ( ) { }
  
  virtual bool notAtEnd ( ) const = 0;
  
  virtual ReferenceBlockIterator& operator++ ( ) {
    if ( !_quietMode && _progressBar != NULL ) (*_progressBar)++;
    return *this;
  };

  virtual qc::CoordType& operator* ( ) = 0;
};
  
/**
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, qc::Dimension Dim>
class GlobalReferenceBlockIterator : public ReferenceBlockIterator<_RealType, Dim> {
protected:
  qc::CoordType _xL, _xR;
  int _numNodes;
  qc::CoordType _refStep;
  qc::CoordType _cur;
public:
  GlobalReferenceBlockIterator ( const qc::CoordType &X0, const qc::CoordType &XEnd,
                                 const int NumNodes,
                                 aol::ProgressBar<> *ProgressBar = NULL, bool Quiet = false,
                                 const int RefStep = 1 )
    : ReferenceBlockIterator<_RealType, Dim> ( ProgressBar, Quiet ),
      _xL ( X0 ), _xR ( XEnd ),
      _numNodes ( NumNodes ), _refStep ( RefStep, RefStep, 1 ),
      _cur ( _xL ) {
    _xR.addToAll ( -1 );
    if ( !this->_quietMode && this->_progressBar != NULL ) this->_progressBar->start ( _numNodes / ( _refStep * _refStep ) );
  }
  
  bool notAtEnd ( ) const {
    if ( _cur[Dim-1] <= _xR[Dim-1] ) return true;
    else {
      if ( !this->_quietMode && this->_progressBar != NULL ) this->_progressBar->finish ( );
      return false;
    }
  }
  
  GlobalReferenceBlockIterator& operator++ ( ) {
    for ( int d=0; d<Dim ; ++d ) {
      _cur[d] += _refStep[d];
      if ( _cur[d] <= _xR[d] ) break;
      if ( _cur[d] - _refStep[d] < _xR[d] ) {
        _cur[d] = _xR[d];
        break;
      }
      if ( d < Dim - 1 ) _cur[d] = _xL[d];
    }
    ReferenceBlockIterator<_RealType, Dim>::operator++ ( );
    
    return *this;
  }
  
  qc::CoordType& operator* ( ) {
    return this->_cur;
  }
};

/**
 * \author mevenkamp
 * \ingroup NonLocalNeighborhoodFilters
 */
template <typename _RealType, qc::Dimension Dim>
class PixelReferenceBlockIterator : public ReferenceBlockIterator<_RealType, Dim> {
protected:
  const aol::RandomAccessContainer<qc::CoordType> &_pixels;
  int _curIdx;
  qc::CoordType _curPos;
public:
  PixelReferenceBlockIterator ( const aol::RandomAccessContainer<qc::CoordType > &Pixels,
                                aol::ProgressBar<> *ProgressBar = NULL, bool Quiet = false )
    : ReferenceBlockIterator<_RealType, Dim> ( ProgressBar, Quiet ),
      _pixels ( Pixels ), _curIdx ( 0 ), _curPos ( _pixels[_curIdx] ) {
    if ( !this->_quietMode && this->_progressBar != NULL ) this->_progressBar->start ( _pixels.size ( ) );
  }
  
  bool notAtEnd ( ) const {
    if ( _curIdx < _pixels.size ( ) ) return true;
    else {
      if ( !this->_quietMode && this->_progressBar != NULL ) this->_progressBar->finish ( );
      return false;
    }
  }
  
  PixelReferenceBlockIterator& operator++ ( ) {
    ++_curIdx;
    if ( notAtEnd ( ) ) _curPos = _pixels[_curIdx];
    ReferenceBlockIterator<_RealType, Dim>::operator++ ( );
    
    return *this;
  }
  
  qc::CoordType& operator* ( ) {
    return _curPos;
  }
};


} // namespace im


#endif
