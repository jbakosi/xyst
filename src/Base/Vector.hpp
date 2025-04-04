// *****************************************************************************
/*!
  \file      src/Base/Vector.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Vector algebra
  \details   Vector algebra.
*/
// *****************************************************************************
#ifndef Vector_h
#define Vector_h

#include <array>
#include <cmath>
#include <vector>

#include "Types.hpp"
#include "Exception.hpp"

namespace tk {

//! Compute the cross-product of two vectors
//! \param[in] v1x x coordinate of the 1st vector
//! \param[in] v1y y coordinate of the 1st vector
//! \param[in] v1z z coordinate of the 1st vector
//! \param[in] v2x x coordinate of the 2nd vector
//! \param[in] v2y y coordinate of the 2nd vector
//! \param[in] v2z z coordinate of the 2nd vector
//! \param[out] rx x coordinate of the product vector
//! \param[out] ry y coordinate of the product vector
//! \param[out] rz z coordinate of the product vector
inline void
cross( real v1x, real v1y, real v1z,
       real v2x, real v2y, real v2z,
       real& rx, real& ry, real& rz )
{
  rx = v1y*v2z - v2y*v1z;
  ry = v1z*v2x - v2z*v1x;
  rz = v1x*v2y - v2x*v1y;
}

//! Compute the cross-product of two vectors
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \return Cross-product
inline std::array< real, 3 >
cross( const std::array< real, 3 >& v1, const std::array< real, 3 >& v2 ) {
  real rx, ry, rz;
  cross( v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], rx, ry, rz );
  return { std::move(rx), std::move(ry), std::move(rz) };
}

//! Compute the cross-product of two vectors divided by a scalar
//! \param[in] v1x x coordinate of the 1st vector
//! \param[in] v1y y coordinate of the 1st vector
//! \param[in] v1z z coordinate of the 1st vector
//! \param[in] v2x x coordinate of the 2nd vector
//! \param[in] v2y y coordinate of the 2nd vector
//! \param[in] v2z z coordinate of the 2nd vector
//! \param[in] j The scalar to divide the product with
//! \param[out] rx x coordinate of the product vector
//! \param[out] ry y coordinate of the product vector
//! \param[out] rz z coordinate of the product vector
inline void
crossdiv( real v1x, real v1y, real v1z,
          real v2x, real v2y, real v2z,
          real j,
          real& rx, real& ry, real& rz )
{
  cross( v1x, v1y, v1z, v2x, v2y, v2z, rx, ry, rz );
  rx /= j;
  ry /= j;
  rz /= j;
}

//! Compute the cross-product of two vectors divided by a scalar
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \param[in] j Scalar to divide each component by
//! \return Cross-product divided by scalar
inline std::array< real, 3 >
crossdiv( const std::array< real, 3 >& v1,
          const std::array< real, 3 >& v2,
          real j )
{
  return {{ (v1[1]*v2[2] - v2[1]*v1[2]) / j,
            (v1[2]*v2[0] - v2[2]*v1[0]) / j,
            (v1[0]*v2[1] - v2[0]*v1[1]) / j }};
}

//! Compute the dot-product of two vectors
//! \param[in] v1x x coordinate of 1st vector
//! \param[in] v1y y coordinate of 1st vector
//! \param[in] v1z z coordinate of 1st vector
//! \param[in] v2x x coordinate of 2nd vector
//! \param[in] v2y y coordinate of 2nd vector
//! \param[in] v2z z coordinate of 2ndt vector
//! \return Dot-product
inline real
dot( real v1x, real v1y, real v1z, real v2x, real v2y, real v2z ) {
  return v1x*v2x + v1y*v2y + v1z*v2z;
}

//! Compute the dot-product of two vectors
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \return Dot-product
inline real
dot( const std::array< real, 3 >& v1, const std::array< real, 3 >& v2 ) {
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

//! Compute length of a vector
//! \param[in] x X coordinate of vector
//! \param[in] y Y coordinate of vector
//! \param[in] z Z coordinate of vector
//! \return length
inline real
length( real x, real y, real z ) {
  return std::sqrt( x*x + y*y + z*z );
}

//! Compute length of a vector
//! \param[in] v vector
//! \return length
inline real
length( const std::array< real, 3 >& v ) {
  return std::sqrt( dot(v,v) );
}

//! Scale vector to unit length
//! \param[in,out] v Vector to normalize
inline void
unit( std::array< real, 3 >& v ) noexcept(ndebug) {
  auto len = length( v );
  // cppcheck-suppress throwInNoexceptFunction
  Assert( len > std::numeric_limits< tk::real >::epsilon(), "div by zero" );
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

//! Compute the triple-product of three vectors
//! \param[in] v1x x coordinate of the 1st vector
//! \param[in] v1y y coordinate of the 1st vector
//! \param[in] v1z z coordinate of the 1st vector
//! \param[in] v2x x coordinate of the 2nd vector
//! \param[in] v2y y coordinate of the 2nd vector
//! \param[in] v2z z coordinate of the 2nd vector
//! \param[in] v3x x coordinate of the 3rd vector
//! \param[in] v3y y coordinate of the 3rd vector
//! \param[in] v3z z coordinate of the 3rd vector
//! \return Scalar value of the triple product
inline tk::real
triple( real v1x, real v1y, real v1z,
        real v2x, real v2y, real v2z,
        real v3x, real v3y, real v3z )
{
  real cx, cy, cz;
  cross( v2x, v2y, v2z, v3x, v3y, v3z, cx, cy, cz );
  return v1x*cx + v1y*cy + v1z*cz;
}

//! Compute the triple-product of three vectors
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \param[in] v3 3rd vector
//! \return Triple-product
inline real
triple( const std::array< real, 3 >& v1,
        const std::array< real, 3 >& v2,
        const std::array< real, 3 >& v3 )
{
  return dot( v1, cross(v2,v3) );
}

//! Rotate vector about X axis
//! \param[in] v Vector to rotate
//! \param[in] angle Angle to use to rotate with
//! \return Rotated vector
inline std::array< real, 3 >
rotateX( const std::array< real, 3 >& v, real angle )
{
  using std::cos;  using std::sin;

  std::array< std::array< real, 3 >, 3 >
    R{{ {{ 1.0,         0.0,          0.0 }},
        {{ 0.0,   cos(angle), -sin(angle) }},
        {{ 0.0,   sin(angle),  cos(angle) }} }};

  return {{ dot(R[0],v), dot(R[1],v), dot(R[2],v) }};
}

//! Rotate vector about Y axis
//! \param[in] v Vector to rotate
//! \param[in] angle Angle to use to rotate with
//! \return Rotated vector
inline std::array< real, 3 >
rotateY( const std::array< real, 3 >& v, real angle )
{
  using std::cos;  using std::sin;

  std::array< std::array< real, 3 >, 3 >
    R{{ {{ cos(angle),  0.0, sin(angle) }},
        {{ 0.0,         1.0,        0.0 }},
        {{ -sin(angle), 0.0, cos(angle) }} }};

  return {{ dot(R[0],v), dot(R[1],v), dot(R[2],v) }};
}

//! Rotate vector about Z axis
//! \param[in] v Vector to rotate
//! \param[in] angle Angle to use to rotate with
//! \return Rotated vector
inline std::array< real, 3 >
rotateZ( const std::array< real, 3 >& v, real angle )
{
  using std::cos;  using std::sin;

  std::array< std::array< real, 3 >, 3 >
    R{{ {{ cos(angle), -sin(angle), 0.0 }},
        {{ sin(angle),  cos(angle), 0.0 }},
        {{ 0.0,         0.0,        1.0 }} }};

  return {{ dot(R[0],v), dot(R[1],v), dot(R[2],v) }};
}

//! Compute the unit normal vector of a triangle
//! \param[in] x1 x coordinate of the 1st vertex of the triangle
//! \param[in] x2 x coordinate of the 2nd vertex of the triangle
//! \param[in] x3 x coordinate of the 3rd vertex of the triangle
//! \param[in] y1 y coordinate of the 1st vertex of the triangle
//! \param[in] y2 y coordinate of the 2nd vertex of the triangle
//! \param[in] y3 y coordinate of the 3rd vertex of the triangle
//! \param[in] z1 z coordinate of the 1st vertex of the triangle
//! \param[in] z2 z coordinate of the 2nd vertex of the triangle
//! \param[in] z3 z coordinate of the 3rd vertex of the triangle
//! \param[out] nx x coordinate of the unit normal
//! \param[out] ny y coordinate of the unit normal
//! \param[out] nz z coordinate of the unit normal
//! \return Triangle area
inline real
normal( real x1, real x2, real x3,
        real y1, real y2, real y3,
        real z1, real z2, real z3,
        real& nx, real& ny, real& nz )
{
  real ax = x2 - x1;
  real ay = y2 - y1;
  real az = z2 - z1;

  real bx = x3 - x1;
  real by = y3 - y1;
  real bz = z3 - z1;

  real n1 =   ay*bz - az*by;
  real n2 = -(ax*bz - az*bx);
  real n3 =   ax*by - ay*bx;

  auto farea = std::sqrt( n1*n1 + n2*n2 + n3*n3 );

  nx = n1/farea;
  ny = n2/farea;
  nz = n3/farea;

  return farea;
}

} // tk::

#endif // Vector_h
