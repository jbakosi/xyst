// *****************************************************************************
/*!
  \file      src/Control/StatCtr.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Types and associated functions to deal with moments and PDFs
  \details   Types and associated functions to deal with statistical moments and
    probability density functions.
*/
// *****************************************************************************
#ifndef StatControl_h
#define StatControl_h

#include "Types.hpp"
#include "Exception.hpp"
#include "PUPUtil.hpp"

namespace tk {
namespace ctr {

//! \brief Moment specifies which type of moment is computed for a quantity in
//!    a Term
enum class Moment : uint8_t { ORDINARY=0,      //!< Full variable
                              CENTRAL          //!< Fluctuation
};

//! \brief Term is a Moment of a quantity with a field ID to be ensemble
//!    averaged
//! \details Internally the numbering of field IDs starts from 0, but presented
//!    to the user, e.g., in screen-output, as starting from 1.
struct Term {
  char var;             //!< Variable name
  uint64_t field;       //!< Field ID
  Moment moment;        //!< Moment type: ordinary, central

  /** @name Pack/Unpack: Serialize Term object for Charm++ */
  ///@{
  //! Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er& p ) {
    p | var;
    p | field;
    PUP::pup( p, moment );
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] t Term object reference
  friend void operator|( PUP::er& p, Term& t ) { t.pup(p); } 
  ///@}

  //! \brief Constructor: initialize all state data
  //! \param[in] v Variable name
  //! \param[in] f Field ID
  //! \param[in] m Moment type enum: Moment::ORDINARY or Moment::CENTRAL
  explicit Term( char v = 0, uint64_t f = 0, Moment m = Moment::ORDINARY ) :
    var( v ), field( f ), moment( m ) {}

  //! \brief Equal operator for, e.g., finding unique elements, used by, e.g.,
  //!    std::unique().
  //! \details Test on field, moment, and var
  //! \param[in] term Term to compare
  //! \return Boolean indicating if term equals 'this'
  bool operator== ( const Term& term ) const {
    if (var == term.var && field == term.field && moment == term.moment)
      return true;
    else
      return false;
  }

  //! \brief Less-than operator for ordering, used by, e.g., std::sort().
  //! \details Test on var, field, and moment.
  //! \param[in] term Term to compare
  //! \return Boolean indicating if term is less than 'this'
  bool operator< ( const Term& term ) const {
    if (var < term.var)
      return true;
    else if (var == term.var && field < term.field)
      return true;
    else if (var == term.var && field == term.field && moment < term.moment)
      return true;
    else
      return false;
  }
};

//! \brief Pack/Unpack: Namespace-scope serialize Term object for Charm++
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] t Term object reference
inline void pup( PUP::er& p, Term& t ) { t.pup(p); }

//! \brief Products are arbitrary number of Terms to be multiplied and ensemble
//!   averaged.
//! \details An example is the scalar flux in x direction which needs two terms
//! for ensemble averaging: (Y-\<Y\>) and (U-\<U\>), then the central moment is
//! \<yu\> = <(Y-\<Y\>)(U-\<U\>)>, another example is the third mixed central
//! moment of three scalars which needs three terms for ensemble averaging:
//! (Y1-\<Y1\>), (Y2-\<Y2\>), and (Y3-\<Y3\>), then the central moment is
//! \<y1y2y3\> = \<(Y1-\<Y1\>)(Y2-\<Y2\>)(Y3-\<Y3\>)\>.
using Product = std::vector< Term >;

//! \brief Case-insensitive character comparison functor
struct CaseInsensitiveCharLess {
  //! Function call operator
  //! \param[in] lhs Left character of the comparitor operand
  //! \param[in] rhs Right character of the comparitor operand
  //! \return Boolean indicating the result of the comparison
  bool operator() ( char lhs, char rhs ) const {
    return std::tolower( lhs ) < std::tolower( rhs );
  }
};

//! \brief Find out if a vector of Terms only contains ordinary moment terms
//! \details If and only if all terms are ordinary, the vector of Terms is
//!    ordinary.
//! \param[in] vec Vector of Terms to check
//! \return Boolean indicating if all terms are ordinary
static inline bool
ordinary( const std::vector< ctr::Term >& vec ) {
  if (std::any_of( vec.cbegin(), vec.cend(),
        []( const ctr::Term& t ){ return t.moment == ctr::Moment::CENTRAL; } ))
    return false;
  else
    return true;
}

//! \brief Find out if a vector of Terms contains any central moment terms
//! \details If any term is central, the vector of Terms is central.
//! \param[in] vec Vector of Terms to check
//! \return Boolean indicating of any term is central
static inline bool
central( const std::vector< ctr::Term >& vec )
{ return !ordinary( vec ); }

//! \brief Probability density function (vector of sample space variables)
using Probability = std::vector< Term >;

//! \brief PDF information bundle
//! \note If the user did not specify extents for a PDF, the corresponding
//!   extents vector still exists but it is empty.
struct PDFInfo {
  const std::string& name;                  //!< PDF identifier, i.e., name
  const std::vector< tk::real >& exts;      //!< Sample space extents
  std::vector< std::string > vars;          //!< List of sample space variables
  std::uint64_t it;                         //!< Iteration count
  tk::real time;                            //!< Time stamp
};

//! \brief Find PDF information, see tk::ctr::PDFInfo
//! \note Size of binsizes, names, pdfs, and exts must all be equal
//! \note idx must be less than the length of binsizes, names, and pdfs
//! \param[in] binsizes Vector of sample space bin sizes for multiple PDFs
//! \param[in] names Vector of PDF names
//! \param[in] exts Vector of sample space extents. Note: if the user did not
//!   specify extents for a PDF, the corresponding extents vector still exists
//!   but it is empty.
//! \param[in] pdfs Vector of PDFs
//! \param[in] m ORDINARY or CENTRAL PDF we are looking for
//! \param[in] idx Index of the PDF within the list of matching (based on D and
//!   m) PDFs requested
//! \param[in] it Iteration count
//! \param[in] time Physical time
//! \return The PDF metadata requested
//! \details Find PDF information given the sample space dimension (template
//!   argument D), ordinary or central PDF (m), and the index within the list of
//!   matching (based on D and m) PDFs requested (idx). This function must find
//!   the PDF, if it does not, it throws an exception.
//! \see walker::Distributor
template< std::size_t D >
PDFInfo pdfInfo( const std::vector< std::vector< tk::real > >& binsizes,
                 const std::vector< std::string >& names,
                 const std::vector< std::vector< tk::real > >& exts,
                 const std::vector< Probability >& pdfs,
                 tk::ctr::Moment m,
                 std::size_t idx,
                 std::uint64_t it,
                 tk::real time )
{
  Assert( binsizes.size() == names.size(),
          "Length of binsizes vector and that of PDF names must equal" );
  Assert( binsizes.size() == pdfs.size(),
          "Length of binsizes vector and that of PDFs must equal" );
  Assert( binsizes.size() == exts.size(),
          "Length of binsizes vector and that of PDF extents must equal" );
  Assert( binsizes.size() > idx, "Indexing out of bounds" );

  std::size_t i = 0;  // will count all PDFs queried
  std::size_t n = 0;  // will count PDFs with sample space dimensions and type
                      // (orindary or central) we are looking for
  for (const auto& bs : binsizes) {
    if ( bs.size() == D &&
         ((m == Moment::ORDINARY && ordinary(pdfs[i])) ||
          (m == Moment::CENTRAL && central(pdfs[i]))) ) ++n;
    if (n == idx+1) {
      std::vector< std::string > vars;
      for (const auto& term : pdfs[i])
        // cppcheck-suppress useStlAlgorithm
        vars.push_back( term.var + std::to_string(term.field+1) );
      return { names[i], exts[i], std::move(vars), it, time };
    }
    ++i;
  }
  Throw( "Cannot find PDF." );
}

//! Construct mean
//! \param[in] var Variable
//! \param[in] c Component number
//! \return Constructed vector< Term > identifying the first ordinary moment
//!   (mean) of field (component) c of variable var
static inline Product
mean( char var, uint64_t c ) {
  tk::ctr::Term m( static_cast<char>(std::toupper(var)), c, Moment::ORDINARY );
  return tk::ctr::Product( { m } );
}

//! Construct variance
//! \param[in] var Variable
//! \param[in] c Component number
//! \return Constructed vector< Term > identifying the second central moment
//!   (variance) of field (component) c of variable var
static inline Product
variance( char var, uint64_t c ) {
  tk::ctr::Term f( static_cast<char>(std::tolower(var)), c, Moment::CENTRAL );
  return tk::ctr::Product( { f, f } );
}

} // ctr::
} // tk::

#endif // StatControl_h
