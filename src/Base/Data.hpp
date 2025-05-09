// *****************************************************************************
/*!
  \file      src/Base/Data.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Generic data storage with different memory layouts
  \details   Generic data storage with different memory layouts. See also the
    rationale discussed in the [design](layout.html) document.
*/
// *****************************************************************************
#ifndef Data_h
#define Data_h

#include <array>
#include <string>
#include <cstdint>
#include <vector>
#include <set>
#include <algorithm>

#include "Types.hpp"
#include "Exception.hpp"

#include "NoWarning/pup_stl.hpp"

namespace tk {

//! Tags for selecting data layout policies
const uint8_t UnkEqComp = 0;
const uint8_t EqCompUnk = 1;

//! Zero-runtime-cost data-layout wrappers with type-based compile-time dispatch
template< uint8_t Layout >
class Data {

  private:
    using ncomp_t = uint64_t;

  public:
    //! Default constructor (required for Charm++ migration)
    explicit Data() : m_vec(), m_nunk(), m_nprop() {}

    //! Constructor
    //! \param[in] nu Number of unknowns to allocate memory for
    //! \param[in] np Total number of properties, i.e., scalar variables or
    //!   components, per unknown
    explicit Data( ncomp_t nu, ncomp_t np ) :
      m_vec( nu*np ),
      m_nunk( nu ),
      m_nprop( np ) {}

    //! Const data access dispatch
    //! \details Public interface to const-ref data access to a single real
    //!   value. Use it as Data(p,c), where p is the unknown index, and c is
    //!   the component index specifying the scalar equation within a system of
    //!   equations. Requirement: component < nprop, unknown < nunk, enforced
    //!   with an assert in DEBUG mode, see also the constructor.
    //! \param[in] unknown Unknown index
    //! \param[in] component Component index, i.e., position of a scalar within
    //!   a system
    //! \return Const reference to data of type tk::real
    const tk::real&
    operator()( ncomp_t unknown, ncomp_t component ) const
    { return access( unknown, component, int2type< Layout >() ); }

    //! Non-const data access dispatch
    //! \details Public interface to non-const-ref data access to a single real
    //!   value. Use it as Data(p,c,o), where p is the unknown index, and c is
    //!   the component index specifying the scalar equation within a system of
    //!   equations. Requirement: component < nprop, unknown < nunk, enforced
    //!   with an assert in DEBUG mode, see also the constructor.
    //! \param[in] unknown Unknown index
    //! \param[in] component Component index, i.e., position of a scalar within
    //!   a system
    //! \return Non-const reference to data of type tk::real
    //! \see "Avoid Duplication in const and Non-const Member Function," and
    //!   "Use const whenever possible," Scott Meyers, Effective C++, 3d ed.
    tk::real&
    operator()( ncomp_t unknown, ncomp_t component ) {
      return const_cast< tk::real& >(
               static_cast< const Data& >( *this ).
                 operator()( unknown, component ) );
    }

    //! Access to number of unknowns
    //! \return Number of unknowns
    ncomp_t nunk() const noexcept { return m_nunk; }

    //! Access to number of properties
    //! \details This is the total number of scalar components per unknown
    //! \return Number of propertes/unknown
    ncomp_t nprop() const noexcept { return m_nprop; }

    //! Extract vector of unknowns given component
    //! \details Requirement: component < nprop, enforced with an assert in
    //!   DEBUG mode, see also the constructor.
    //! \param[in] component Component index, i.e., position of a scalar within
    //!   a system
    //! \return A vector of unknowns given by component (length: nunk(), i.e.,
    //!   the first constructor argument)
    std::vector< tk::real >
    extract( ncomp_t component ) const {
      std::vector< tk::real > w( m_nunk );
      for (ncomp_t i=0; i<m_nunk; ++i) w[i] = operator()( i, component );
      return w;
    }

    //! Extract all components for unknown
    //! \details Requirement: unknown < nunk, enforced with an assert in DEBUG
    //!   mode, see also the constructor.
    //! \param[in] unknown Index of unknown
    //! \return A vector of components for a single unknown (length: nprop,
    //!   i.e., the second constructor argument)
    std::vector< tk::real >
    operator[]( ncomp_t unknown ) const {
      std::vector< tk::real > w( m_nprop );
      for (ncomp_t i=0; i<m_nprop; ++i) w[i] = operator()( unknown, i );
      return w;
    }

    //! Extract (a copy of) four values of unknowns
    //! \details Requirement: component < nprop, for all N[i] < nunk, enforced
    //!   with an assert in DEBUG mode, see also the constructor.
    //! \param[in] component Component index, i.e., position of a scalar within
    //!   a system
    //! \param[in] N Indices of the 4 unknowns
    //! \return Array of the four values of component
    std::array< tk::real, 4 >
    extract( ncomp_t component, const std::array< ncomp_t, 4 >& N ) const {
      auto p = cptr( component );
      return {{ var(p,N[0]), var(p,N[1]), var(p,N[2]), var(p,N[3]) }};
    }

    //! Const-ref accessor to underlying raw data as a std::vector
    //! \return Constant reference to underlying raw data
    const std::vector< tk::real >& vec() const { return m_vec; }

    //! Non-const-ref accessor to underlying raw data as a std::vector
    //! \return Non-constant reference to underlying raw data
    std::vector< tk::real >& vec() { return m_vec; }

    //! Compound operator-=
    //! \param[in] rhs Data object to subtract
    //! \return Reference to ourselves after subtraction
    Data< Layout >& operator-= ( const Data< Layout >& rhs ) {
      Assert( rhs.nunk() == m_nunk, "Incorrect number of unknowns" );
      Assert( rhs.nprop() == m_nprop, "Incorrect number of properties" );
      std::transform( rhs.vec().cbegin(), rhs.vec().cend(),
                      m_vec.cbegin(), m_vec.begin(),
                      []( tk::real s, tk::real d ){ return d-s; } );
      return *this;
    }
    //! Operator -
    //! \param[in] rhs Data object to subtract
    //! \return Copy of Data object after rhs has been subtracted
    //! \details Implemented in terms of compound operator-=
    Data< Layout > operator- ( const Data< Layout >& rhs )
    const { return Data< Layout >( *this ) -= rhs; }

    //! Compound operator+=
    //! \param[in] rhs Data object to add
    //! \return Reference to ourselves after addition
    Data< Layout >& operator+= ( const Data< Layout >& rhs ) {
      Assert( rhs.nunk() == m_nunk, "Incorrect number of unknowns" );
      Assert( rhs.nprop() == m_nprop, "Incorrect number of properties" );
      std::transform( rhs.vec().cbegin(), rhs.vec().cend(),
                      m_vec.cbegin(), m_vec.begin(),
                      []( tk::real s, tk::real d ){ return d+s; } );
      return *this;
    }
    //! Operator +
    //! \param[in] rhs Data object to add
    //! \return Copy of Data object after rhs has been multiplied with
    //! \details Implemented in terms of compound operator+=
    Data< Layout > operator+ ( const Data< Layout >& rhs )
    const { return Data< Layout >( *this ) += rhs; }

    //! Compound operator*= multiplying by another Data object item by item
    //! \param[in] rhs Data object to multiply with
    //! \return Reference to ourselves after multiplication
    Data< Layout >& operator*= ( const Data< Layout >& rhs ) {
      Assert( rhs.nunk() == m_nunk, "Incorrect number of unknowns" );
      Assert( rhs.nprop() == m_nprop, "Incorrect number of properties" );
      std::transform( rhs.vec().cbegin(), rhs.vec().cend(),
                      m_vec.cbegin(), m_vec.begin(),
                      []( tk::real s, tk::real d ){ return d*s; } );
      return *this;
    }
    //! Operator * multiplying by another Data object item by item
    //! \param[in] rhs Data object to multiply with
    //! \return Copy of Data object after rhs has been multiplied with
    //! \details Implemented in terms of compound operator*=
    Data< Layout > operator* ( const Data< Layout >& rhs )
    const { return Data< Layout >( *this ) *= rhs; }

    //! Compound operator*= multiplying all items by a scalar
    //! \param[in] rhs Scalar to multiply with
    //! \return Reference to ourselves after multiplication
    Data< Layout >& operator*= ( tk::real rhs ) {
      // cppcheck-suppress useStlAlgorithm
      for (auto& v : m_vec) v *= rhs;
      return *this;
    }
    //! Operator * multiplying all items by a scalar
    //! \param[in] rhs Scalar to multiply with
    //! \return Copy of Data object after rhs has been multiplied with
    //! \details Implemented in terms of compound operator*=
    Data< Layout > operator* ( tk::real rhs )
    const { return Data< Layout >( *this ) *= rhs; }

    //! Compound operator/=
    //! \param[in] rhs Data object to divide by
    //! \return Reference to ourselves after division
    Data< Layout >& operator/= ( const Data< Layout >& rhs ) {
      Assert( rhs.nunk() == m_nunk, "Incorrect number of unknowns" );
      Assert( rhs.nprop() == m_nprop, "Incorrect number of properties" );
      std::transform( rhs.vec().cbegin(), rhs.vec().cend(),
                      m_vec.cbegin(), m_vec.begin(),
                      []( tk::real s, tk::real d ){ return d/s; } );
      return *this;
    }
    //! Operator /
    //! \param[in] rhs Data object to divide by
    //! \return Copy of Data object after rhs has been divided by
    //! \details Implemented in terms of compound operator/=
    Data< Layout > operator/ ( const Data< Layout >& rhs )
    const { return Data< Layout >( *this ) /= rhs; }

    //! Compound operator/= dividing all items by a scalar
    //! \param[in] rhs Scalar to divide with
    //! \return Reference to ourselves after division
    Data< Layout >& operator/= ( tk::real rhs ) {
      // cppcheck-suppress useStlAlgorithm
      for (auto& v : m_vec) v /= rhs;
      return *this;
    }
    //! Operator / dividing all items by a scalar
    //! \param[in] rhs Scalar to divide with
    //! \return Copy of Data object after rhs has been divided by
    //! \details Implemented in terms of compound operator/=
    Data< Layout > operator/ ( tk::real rhs )
    const { return Data< Layout >( *this ) /= rhs; }

    //! Add new unknown at the end of the container
    //! \param[in] prop Vector of properties to initialize the new unknown with
    void push_back( const std::vector< tk::real >& prop )
    { return push_back( prop, int2type< Layout >() ); }

    //! Resize data store to contain 'count' elements
    //! \param[in] count Resize store to contain count * nprop elements
    //! \param[in] value Value to initialize new data with (default: 0.0)
    //! \note This works for both shrinking and enlarging, as this simply
    //!   translates to std::vector::resize(). Note that count changes, nprop
    //!   does not, see the private overload resize().
    void resize( std::size_t count, tk::real value = 0.0 )
    { resize( count, value, int2type< Layout >() ); }

    //! Resize data given number of unknowns and number of properties
    //! \param[in] nu Number of unknowns to allocate memory for
    //! \param[in] np Total number of properties, i.e., scalar variables or
    //!   components, per unknown
    void resize( std::size_t nu, std::size_t np ) {
      m_vec.resize( nu * np );
      m_nunk = nu;
      m_nprop = np;
    }

    //! Remove a number of unknowns
    //! \param[in] unknown Set of indices of unknowns to remove
    void rm( const std::set< ncomp_t >& unknown ) {
      auto remove = [ &unknown ]( std::size_t i ) -> bool {
        if (unknown.find(i) != end(unknown)) return true;
        return false;
      };
      std::size_t last = 0;
      for(std::size_t i=0; i<m_nunk; ++i, ++last) {
        while( remove(i) ) ++i;
        if (i >= m_nunk) break;
        for (ncomp_t p = 0; p<m_nprop; ++p)
          m_vec[ last*m_nprop+p ] = m_vec[ i*m_nprop+p ];
      }
      m_vec.resize( last*m_nprop );
      m_nunk -= unknown.size();
    }

    //! Fill vector of unknowns with the same value
    //! \details Requirement: component < nprop, enforced with an assert in
    //!   DEBUG mode, see also the constructor.
    //! \param[in] component Component index, i.e., position of a scalar within
    //!   a system
    //! \param[in] value Value to fill vector of unknowns with
    inline void fill( ncomp_t component, tk::real value ) {
      auto p = cptr( component );
      for (ncomp_t i=0; i<m_nunk; ++i) var(p,i) = value;
    }

    //! Fill full data storage with value
    //! \param[in] value Value to fill data with
    void fill( tk::real value )
    { std::fill( begin(m_vec), end(m_vec), value ); }

    //! Check if vector of unknowns is empty
    bool empty() const noexcept { return m_vec.empty(); }

    //! Layout name dispatch
    //! \return The name of the data layout used
    static std::string layout() { return layout( int2type< Layout >() ); }

    /** @name Pack/Unpack: Serialize Data object for Charm++ */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    // cppcheck-suppress constParameter
    void pup( PUP::er &p ) {
      p | m_vec;
      p | m_nunk;
      p | m_nprop;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] d DataLyaout object reference
    friend void operator|( PUP::er& p, Data& d ) { d.pup(p); }
    //@}

  private:
    //! Transform a compile-time uint8_t into a type, used for dispatch
    //! \see A. Alexandrescu, Modern C++ Design: Generic Programming and Design
    //!   Patterns Applied, Addison-Wesley Professional, 2001.
    template< uint8_t m > struct int2type { enum { value = m }; };

    //! Const ptr to physical variable access dispatch
    //! \details Public interface to the first half of a physical variable
    //!   access. cptr() and var() are two member functions intended to be used
    //!   together in case when component would be expensive to compute for data
    //!   access via the function call operator, i.e., cptr(), can be used to
    //!   pre-compute part of the address, which returns a pointer and var() can
    //!   be used to finish the data access using the pointer returned by
    //!   cptr(). In other words, cptr() returns part of the address known based
    //!   on component and intended to be used in a setup phase. Then var()
    //!   takes this partial address and finishes the address calculation given
    //!   the unknown id. Thus the following two data accesses are equivalent
    //!   (modulo constness):
    //!   * real& value = operator()( unk, comp ); and
    //!   * const real* p = cptr( comp ); and
    //!     const real& value = var( p, unk ); or real& value = var( p, unk );
    //!   Requirement: component < nprop, enforced with an assert in DEBUG mode,
    //!   see also the constructor.
    //! \param[in] component Component index, i.e., position of a scalar within
    //!   a system
    //! \return Pointer to data of type tk::real for use with var()
    const tk::real*
    cptr( ncomp_t component ) const
    { return cptr( component, int2type< Layout >() ); }

    //! Const-ref data-access dispatch
    //! \details Public interface to the second half of a physical variable
    //!   access. cptr() and var() are two member functions intended to be used
    //!   together in case when component would be expensive to compute for data
    //!   access via the function call operator, i.e., cptr(), can be used to
    //!   pre-compute part of the address, which returns a pointer and var() can
    //!   be used to finish the data access using the pointer returned by
    //!   cptr(). In other words, cptr() returns part of the address known based
    //!   on component and intended to be used in a setup phase. Then var()
    //!   takes this partial address and finishes the address calculation given
    //!   the unknown id. Thus the following two data accesses are equivalent
    //!   (modulo constness):
    //!   * real& value = operator()( unk, comp ); and
    //!   * const real* p = cptr( comp ); and
    //!     const real& value = var( p, unk ); or real& value = var( p, unk );
    //!   Requirement: unknown < nunk, enforced with an assert in DEBUG mode,
    //!   see also the constructor.
    //! \param[in] pt Pointer to data of type tk::real as returned from cptr()
    //! \param[in] unknown Unknown index
    //! \return Const reference to data of type tk::real
    const tk::real&
    var( const tk::real* pt, ncomp_t unknown ) const
    { return var( pt, unknown, int2type< Layout >() ); }

    //! Non-const-ref data-access dispatch
    //! \details Public interface to the second half of a physical variable
    //!   access. cptr() and var() are two member functions intended to be used
    //!   together in case when component would be expensive to compute for data
    //!   access via the function call operator, i.e., cptr(), can be used to
    //!   pre-compute part of the address, which returns a pointer and var() can
    //!   be used to finish the data access using the pointer returned by
    //!   cptr(). In other words, cptr() returns part of the address known based
    //!   on component and intended to be used in a setup phase. Then var()
    //!   takes this partial address and finishes the address calculation given
    //!   the unknown id. Thus the following two data accesses are equivalent
    //!   (modulo constness):
    //!   * real& value = operator()( unk, comp ); and
    //!   * const real* p = cptr( comp ); and
    //!     const real& value = var( p, unk ); or real& value = var( p, unk );
    //!   Requirement: unknown < nunk, enforced with an assert in DEBUG mode,
    //!   see also the constructor.
    //! \param[in] pt Pointer to data of type tk::real as returned from cptr()
    //! \param[in] unknown Unknown index
    //! \return Non-const reference to data of type tk::real
    //! \see "Avoid Duplication in const and Non-const Member Function," and
    //!   "Use const whenever possible," Scott Meyers, Effective C++, 3d ed.
    tk::real&
    var( const tk::real* pt, ncomp_t unknown ) {
      return const_cast< tk::real& >(
               static_cast< const Data& >( *this ).var( pt, unknown ) );
    }

    //! Overloads for the various const data accesses
    //! \details Requirement: component < nprop, unknown < nunk, enforced with
    //!   an assert in DEBUG mode, see also the constructor.
    //! \param[in] unknown Unknown index
    //! \param[in] component Component index, i.e., position of a scalar within
    //!   a system
    //! \return Const reference to data of type tk::real
    //! \see A. Alexandrescu, Modern C++ Design: Generic Programming and Design
    //!   Patterns Applied, Addison-Wesley Professional, 2001.
    const tk::real&
    access( ncomp_t unknown, ncomp_t component, int2type< UnkEqComp > ) const {
      Assert( component < m_nprop,
              "Out-of-bounds access: component < number of properties" );
      Assert( unknown < m_nunk,
              "Out-of-bounds access: unknown < number of unknowns" );
      return m_vec[ unknown*m_nprop + component ];
    }
    const tk::real&
    access( ncomp_t unknown, ncomp_t component, int2type< EqCompUnk > ) const {
      Assert( component < m_nprop,
              "Out-of-bounds access: component < number of properties" );
      Assert( unknown < m_nunk,
              "Out-of-bounds access: unknown < number of unknowns" );
      return m_vec[ component*m_nunk + unknown ];
    }

    // Overloads for the various const ptr to physical variable accesses
    //! \details Requirement: component < nprop, unknown < nunk, enforced with
    //!   an assert in DEBUG mode, see also the constructor.
    //! \param[in] component Component index, i.e., position of a scalar within
    //!   a system
    //! \return Pointer to data of type tk::real for use with var()
    //! \see A. Alexandrescu, Modern C++ Design: Generic Programming and Design
    //!   Patterns Applied, Addison-Wesley Professional, 2001.
    const tk::real*
    cptr( ncomp_t component, int2type< UnkEqComp > ) const {
      Assert( component < m_nprop,
              "Out-of-bounds access: component < number of properties" );
      return m_vec.data() + component;
    }
    const tk::real*
    cptr( ncomp_t component, int2type< EqCompUnk > ) const {
      Assert( component < m_nprop,
              "Out-of-bounds access: component < number of properties" );
      return m_vec.data() + component*m_nunk;
    }

    // Overloads for the various const physical variable accesses
    //!   Requirement: unknown < nunk, enforced with an assert in DEBUG mode,
    //!   see also the constructor.
    //! \param[in] pt Pointer to data of type tk::real as returned from cptr()
    //! \param[in] unknown Unknown index
    //! \return Const reference to data of type tk::real
    //! \see A. Alexandrescu, Modern C++ Design: Generic Programming and Design
    //!   Patterns Applied, Addison-Wesley Professional, 2001.
    inline const tk::real&
    var( const tk::real* const pt, ncomp_t unknown, int2type< UnkEqComp > )
    const {
      Assert( unknown < m_nunk, "Out-of-bounds access: unknown < number of "
              "unknowns" );
      return *(pt + unknown*m_nprop);
    }
    inline const tk::real&
    var( const tk::real* const pt, ncomp_t unknown, int2type< EqCompUnk > )
    const {
      Assert( unknown < m_nunk, "Out-of-bounds access: unknown < number of "
              "unknowns" );
      return *(pt + unknown);
    }

    //! Add new unknown
    //! \param[in] prop Vector of properties to initialize the new unknown with
    //! \note Only the UnkEqComp overload is provided as this operation would be
    //!   too inefficient with the EqCompUnk data layout.
    void push_back( const std::vector< tk::real >& prop, int2type< UnkEqComp > )
    {
      Assert( prop.size() == m_nprop, "Incorrect number of properties" );
      m_vec.resize( (m_nunk+1) * m_nprop );
      ncomp_t u = m_nunk;
      ++m_nunk;
      for (ncomp_t i=0; i<m_nprop; ++i) operator()( u, i ) = prop[i];
    }

    void push_back( const std::vector< tk::real >&, int2type< EqCompUnk > )
    { Throw( "Not implented. It would be inefficient" ); }

    //! Resize data store to contain 'count' elements
    //! \param[in] count Resize store to contain 'count' elements
    //! \param[in] value Value to initialize new data with
    //! \note Only the UnkEqComp overload is provided as this operation would be
    //!   too inefficient with the EqCompUnk data layout.
    //! \note This works for both shrinking and enlarging, as this simply
    //!   translates to std::vector::resize().
    void resize( std::size_t count, tk::real value, int2type< UnkEqComp > ) {
      m_vec.resize( count * m_nprop, value );
      m_nunk = count;
    }

    void resize( std::size_t, tk::real, int2type< EqCompUnk > ) {
      Throw( "Not implemented. It would be inefficient" );
    }

    // Overloads for the name-queries of data lauouts
    //! \return The name of the data layout used
    //! \see A. Alexandrescu, Modern C++ Design: Generic Programming and Design
    //!   Patterns Applied, Addison-Wesley Professional, 2001.
    static std::string layout( int2type< UnkEqComp > )
    { return "unknown-major"; }
    static std::string layout( int2type< EqCompUnk > )
    { return "equation-major"; }

    std::vector< tk::real > m_vec;      //!< Data pointer
    ncomp_t m_nunk;                     //!< Number of unknowns
    ncomp_t m_nprop;                    //!< Number of properties/unknown
};

//! Operator * multiplying all items by a scalar from the left
//! \param[in] lhs Scalar to multiply with
//! \param[in] rhs Date object to multiply
//! \return New Data object with all items multipled with lhs
template< uint8_t Layout >
Data< Layout > operator* ( tk::real lhs, const Data< Layout >& rhs ) {
  return Data< Layout >( rhs ) *= lhs;
}

//! Operator min between two Data objects
//! \param[in] a 1st Data object
//! \param[in] b 2nd Data object
//! \return New Data object containing the minimum of all values for each
//!   value in _a_ and _b_
//! \note The Data objects _a_ and _b_ must have the same number of
//!   unknowns and properties.
//! \note As opposed to std::min, this function creates and returns a new object
//!   instead of returning a reference to one of the operands.
template< uint8_t Layout >
Data< Layout > min( const Data< Layout >& a, const Data< Layout >& b ) {
  Assert( a.nunk() == b.nunk(), "Number of unknowns unequal" );
  Assert( a.nprop() == b.nprop(), "Number of properties unequal" );
  Data< Layout > r( a.nunk(), a.nprop() );
  std::transform( a.vec().cbegin(), a.vec().cend(),
                  b.vec().cbegin(), r.vec().begin(),
                  []( tk::real s, tk::real d ){ return std::min(s,d); } );

  return r;
}

//! Operator max between two Data objects
//! \param[in] a 1st Data object
//! \param[in] b 2nd Data object
//! \return New Data object containing the maximum of all values for each
//!   value in _a_ and _b_
//! \note The Data objects _a_ and _b_ must have the same number of
//!   unknowns and properties.
//! \note As opposed to std::max, this function creates and returns a new object
//!   instead of returning a reference to one of the operands.
template< uint8_t Layout >
Data< Layout > max( const Data< Layout >& a, const Data< Layout >& b ) {
  Assert( a.nunk() == b.nunk(), "Number of unknowns unequal" );
  Assert( a.nprop() == b.nprop(), "Number of properties unequal" );
  Data< Layout > r( a.nunk(), a.nprop() );
  std::transform( a.vec().cbegin(), a.vec().cend(),
                  b.vec().cbegin(), r.vec().begin(),
                  []( tk::real s, tk::real d ){ return std::max(s,d); } );
  return r;
}

//! Operator == between two Data objects
//! \param[in] lhs Data object to compare
//! \param[in] rhs Data object to compare
//! \return True if all entries are equal up to epsilon
template< uint8_t Layout >
bool operator== ( const Data< Layout >& lhs, const Data< Layout >& rhs ) {
  Assert( rhs.nunk() == lhs.nunk(), "Incorrect number of unknowns" );
  Assert( rhs.nprop() == lhs.nprop(), "Incorrect number of properties" );
  auto l = lhs.vec().cbegin();
  auto r = rhs.vec().cbegin();
  while (l != lhs.vec().cend()) {
    if (std::abs(*l - *r) > std::numeric_limits< tk::real >::epsilon())
     return false;
    ++l; ++r;
  }
  return true;
}

//! Operator != between two Data objects
//! \param[in] lhs Data object to compare
//! \param[in] rhs Data object to compare
//! \return True if all entries are unequal up to epsilon
template< uint8_t Layout >
bool operator!= ( const Data< Layout >& lhs, const Data< Layout >& rhs )
{ return !(lhs == rhs); }

//! Compute the maximum difference between the elements of two Data objects
//! \param[in] lhs 1st Data object
//! \param[in] rhs 2nd Data object
//! \return The index, i.e., the raw position, of and the largest absolute value
//!   of the difference between all corresponding elements of _lhs_ and _rhs_.
//! \details The position returned is the position in the underlying raw data
//!   structure, independent of components. If lhs == rhs with precision
//!   std::numeric_limits< tk::real >::epsilon(), a pair of (0,0.0) is returned.
//! \note The Data objects _lhs_ and _rhs_ must have the same number of
//!   unknowns and properties.
template< uint8_t Layout >
std::pair< std::size_t, tk::real >
maxdiff( const Data< Layout >& lhs, const Data< Layout >& rhs ) {
  Assert( lhs.nunk() == rhs.nunk(), "Number of unknowns unequal" );
  Assert( lhs.nprop() == rhs.nprop(), "Number of properties unequal" );
  auto l = lhs.vec().cbegin();
  auto r = rhs.vec().cbegin();
  std::pair< std::size_t, tk::real > m( 0, std::abs(*l - *r) );
  ++l; ++r;
  while (l != lhs.vec().cend()) {
    const auto d = std::abs(*l - *r);
    if (d > m.second) m = { std::distance(lhs.vec().cbegin(),l), d };
    ++l; ++r;
  }
  return m;
}

} // tk::

#endif // Data_h
