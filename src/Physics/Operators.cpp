// *****************************************************************************
/*!
  \file      src/Physics/Operators.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Compressible single-material flow using continuous Galerkin
  \details   Physics operators governing compressible single-material flow using
             continuous Galerkin finite elements.
*/
// *****************************************************************************

#include "Tags.hpp"
#include "Vector.hpp"
#include "Around.hpp"
#include "DerivedData.hpp"
#include "EOS.hpp"
#include "Operators.hpp"
#include "Problems.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

void
ICBoxNodes( const tk::UnsMesh::Coords& coord,
                     std::vector< std::unordered_set< std::size_t > >& inbox )
// *****************************************************************************
//! Determine nodes that lie inside the user-defined IC box
//! \param[in] coord Mesh node coordinates
//! \param[in,out] inbox List of nodes at which box user ICs are set for
//!    each IC box
// *****************************************************************************
{
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  const auto& icbox =
    g_inputdeck.get< tag::param, tag::compflow, tag::ic, tag::box >();

  if (icbox.size() > 0) {
    std::size_t bcnt = 0;
    for (const auto& b : icbox[0]) {
      inbox.emplace_back();
      std::vector< tk::real > box
        { b.get< tag::xmin >(), b.get< tag::xmax >(),
          b.get< tag::ymin >(), b.get< tag::ymax >(),
          b.get< tag::zmin >(), b.get< tag::zmax >() };

      const auto eps = std::numeric_limits< tk::real >::epsilon();
      // Determine which nodes lie in the IC box
      if ( std::any_of( begin(box), end(box), [=](auto p)
                        { return abs(p) > eps; } ) )
      {
        for (std::size_t i=0; i<x.size(); ++i) {
          if ( x[i]>box[0] && x[i]<box[1] && y[i]>box[2] && y[i]<box[3] &&
            z[i]>box[4] && z[i]<box[5] )
          {
            inbox[bcnt].insert( i );
          }
        }
      }
      ++bcnt;
    }
  }
}

#pragma omp declare simd
static bool
stagPoint( tk::real x, tk::real y, tk::real z )
// *****************************************************************************
//! Decide if point is a stagnation point
//! \param[in] x X mesh point coordinates to query
//! \param[in] y Y mesh point coordinates to query
//! \param[in] z Z mesh point coordinates to query
//! \return True if point is configured as a stagnation point by the user
// *****************************************************************************
{
  const auto& stagcnf = g_inputdeck.specialBC< tag::compflow, tag::stag >( 0 );
  const auto& pnt = std::get< 0 >( stagcnf );
  const auto& rad = std::get< 1 >( stagcnf );
  for (std::size_t i=0; i<rad.size(); ++i) {
    if (tk::length( x-pnt[i*3+0], y-pnt[i*3+1], z-pnt[i*3+2] ) < rad[i])
      return true;
  }
  return false;
}

#pragma omp declare simd
static bool
skipPoint( tk::real x, tk::real y, tk::real z )
// *****************************************************************************
//! Decide if point is a skip-BC point
//! \param[in] x X mesh point coordinates to query
//! \param[in] y Y mesh point coordinates to query
//! \param[in] z Z mesh point coordinates to query
//! \return True if point is configured as a skip-BC point by the user
// *****************************************************************************
{
  const auto& skipcnf = g_inputdeck.specialBC< tag::compflow, tag::skip >( 0 );
  const auto& pnt = std::get< 0 >( skipcnf );
  const auto& rad = std::get< 1 >( skipcnf );
  for (std::size_t i=0; i<rad.size(); ++i) {
    if (tk::length( x-pnt[i*3+0], y-pnt[i*3+1], z-pnt[i*3+2] ) < rad[i])
      return true;
  }
  return false;
}

std::function< std::vector< tk::real >
             ( tk::real, tk::real, tk::real, tk::real ) >
IC()
// *****************************************************************************
//  Query user config and assign function to set initial conditions
//! \return The function to call to set initial conditions
// *****************************************************************************
{
  const auto& problem =
    g_inputdeck.get< tag::param, tag::compflow, tag::problem >()[ 0 ];

  std::function< std::vector< tk::real >
               ( tk::real, tk::real, tk::real, tk::real ) > ic;

  if (problem == ctr::ProblemType::USER_DEFINED)
    ic = userdef::ic;
  else if (problem == ctr::ProblemType::NONLINEAR_ENERGY_GROWTH)
    ic = nonlinear_energy_growth::ic;
  else if (problem == ctr::ProblemType::RAYLEIGH_TAYLOR)
    ic = rayleigh_taylor::ic;
  else if (problem == ctr::ProblemType::SEDOV)
    ic = sedov::ic;
  else if (problem == ctr::ProblemType::SOD)
    ic = sod::ic;
  else if (problem == ctr::ProblemType::TAYLOR_GREEN)
    ic = taylor_green::ic;
  else if (problem == ctr::ProblemType::VORTICAL_FLOW)
    ic = vortical_flow::ic;
  else
    Throw( "problem type ic not hooked up" );

  return ic;
}

void
initialize( const std::array< std::vector< tk::real >, 3 >& coord,
            tk::Fields& U,
            tk::real t )
// *****************************************************************************
//! Initalize the compressible flow equations, prepare for time integration
//! \param[in] coord Mesh node coordinates
//! \param[in,out] U Array of unknowns
//! \param[in] t Physical time
// *****************************************************************************
{
  Assert( coord[0].size() == U.nunk(), "Size mismatch" );

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  auto ic = IC();

  // Set initial conditions dependeing on problem configured
  for (std::size_t i=0; i<x.size(); ++i) {
    auto s = ic( x[i], y[i], z[i], t );
    U(i,0,0) = s[0];    // rho
    if (!skipPoint(x[i],y[i],z[i]) && stagPoint(x[i],y[i],z[i])) {
      U(i,1,0) = U(i,2,0) = U(i,3,0) = 0.0;
    } else {
      U(i,1,0) = s[1]; // rho * u
      U(i,2,0) = s[2]; // rho * v
      U(i,3,0) = s[3]; // rho * w
    }
    U(i,4,0) = s[4]; // rho * e, e: total = kinetic + internal
  }
}

static
std::function< void( tk::real, tk::real, tk::real, tk::real,
                     tk::real&, tk::real&, tk::real&, tk::real&, tk::real& ) >
SRC()
// *****************************************************************************
//  Query user config and assign function to add a source term
//! \return The function to call to evaluate a problem-sepcific source term
// *****************************************************************************
{
  const auto& problem =
    g_inputdeck.get< tag::param, tag::compflow, tag::problem >()[ 0 ];

  std::function<
    void( tk::real, tk::real, tk::real, tk::real,
          tk::real&, tk::real&, tk::real&, tk::real&, tk::real& ) > src;

  if (problem == ctr::ProblemType::NONLINEAR_ENERGY_GROWTH)
    src = nonlinear_energy_growth::src;
  else if (problem == ctr::ProblemType::RAYLEIGH_TAYLOR)
    src = rayleigh_taylor::src;
  else if (problem == ctr::ProblemType::TAYLOR_GREEN)
    src = taylor_green::src;
  else if (problem == ctr::ProblemType::VORTICAL_FLOW)
    src = vortical_flow::src;

  return src;
}

void
bndgrad( const std::array< std::vector< tk::real >, 3 >& coord,
         const std::vector< std::size_t >& inpoel,
         const std::vector< std::size_t >& bndel,
         const std::vector< std::size_t >& gid,
         const std::unordered_map< std::size_t, std::size_t >& bid,
         const tk::Fields& U,
         tk::Fields& G )
// *****************************************************************************
//  Compute nodal gradients of primitive variables for along chare boundaries
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] bndel List of elements contributing to chare-boundary nodes
//! \param[in] gid Local->global node id map
//! \param[in] bid Local chare-boundary node ids (value) associated to
//!    global node ids (key)
//! \param[in] U Solution vector at recent time step
//! \param[in,out] G Nodal gradients of primitive variables
//! \details This function only computes local contributions to gradients
//!   at chare-boundary nodes. Internal node gradients are calculated as
//!   required, and do not need to be stored.
// *****************************************************************************
{
  Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
          "vector at recent time step incorrect" );

  G.fill( 0.0 );

  // access node cooordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (auto e : bndel) {  // elements contributing to chare boundary nodes
    // access node IDs
    std::size_t N[4] =
      { inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] };
    // compute element Jacobi determinant, J = 6V
    tk::real bax = x[N[1]]-x[N[0]];
    tk::real bay = y[N[1]]-y[N[0]];
    tk::real baz = z[N[1]]-z[N[0]];
    tk::real cax = x[N[2]]-x[N[0]];
    tk::real cay = y[N[2]]-y[N[0]];
    tk::real caz = z[N[2]]-z[N[0]];
    tk::real dax = x[N[3]]-x[N[0]];
    tk::real day = y[N[3]]-y[N[0]];
    tk::real daz = z[N[3]]-z[N[0]];
    auto J = tk::triple( bax, bay, baz, cax, cay, caz, dax, day, daz );
    ErrChk( J > 0, "Element Jacobian non-positive" );
    auto J24 = J/24.0;
    // shape function derivatives, nnode*ndim [4][3]
    tk::real g[4][3];
    tk::crossdiv( cax, cay, caz, dax, day, daz, J,
                  g[1][0], g[1][1], g[1][2] );
    tk::crossdiv( dax, day, daz, bax, bay, baz, J,
                  g[2][0], g[2][1], g[2][2] );
    tk::crossdiv( bax, bay, baz, cax, cay, caz, J,
                  g[3][0], g[3][1], g[3][2] );
    for (std::size_t i=0; i<3; ++i)
      g[0][i] = -g[1][i] - g[2][i] - g[3][i];
    // scatter-add gradient contributions to boundary nodes
    for (std::size_t a=0; a<4; ++a) {
      auto i = bid.find( gid[N[a]] );
      if (i != end(bid)) {
        tk::real u[5];
        for (std::size_t b=0; b<4; ++b) {
          u[0] = U(N[b],0,0);
          u[1] = U(N[b],1,0)/u[0];
          u[2] = U(N[b],2,0)/u[0];
          u[3] = U(N[b],3,0)/u[0];
          u[4] = U(N[b],4,0)/u[0]
                 - 0.5*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3]);
          if ( !skipPoint(x[N[b]],y[N[b]],z[N[b]]) &&
               stagPoint(x[N[b]],y[N[b]],z[N[b]]) )
          {
            u[1] = u[2] = u[3] = 0.0;
          }
          for (std::size_t c=0; c<5; ++c)
            for (std::size_t j=0; j<3; ++j)
              G( i->second, c*3+j, 0 ) += J24 * g[b][j] * u[c];
        }
      }
    }
  }
}

static tk::Fields
grad( const std::array< std::vector< tk::real >, 3 >& coord,
      const std::vector< std::size_t >& inpoel,
      const std::pair< std::vector< std::size_t >,
                       std::vector< std::size_t > >& esup,
      const std::unordered_map< std::size_t, std::size_t >& lid,
      const std::unordered_map< std::size_t, std::size_t >& bid,
      const std::vector< tk::real >& vol,
      const tk::Fields& U,
      const tk::Fields& bG )
// *****************************************************************************
//! \brief Compute/assemble nodal gradients of primitive variables in all points
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] esup Elements surrounding points
//! \param[in] lid Global->local node ids
//! \param[in] bid Local chare-boundary node ids (value) associated to
//!    global node ids (key)
//! \param[in] vol Nodal volumes
//! \param[in] U Solution vector at recent time step
//! \param[in] bG Nodal gradients in chare boundaries
//! \return Gradients of primitive variables in all mesh points
// *****************************************************************************
{
  // zero storage of nodal gradients of primitive variables
  tk::Fields G( U.nunk(), bG.nprop() );  
  G.fill( 0.0 );

  // access node cooordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // compute gradients of primitive variables in points
  auto npoin = U.nunk();
  #pragma omp simd
  for (std::size_t p=0; p<npoin; ++p)
    for (auto e : tk::Around(esup,p)) {
      // access node IDs
      std::size_t N[4] =
        { inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] };
      // compute element Jacobi determinant, J = 6V
      tk::real bax = x[N[1]]-x[N[0]];
      tk::real bay = y[N[1]]-y[N[0]];
      tk::real baz = z[N[1]]-z[N[0]];
      tk::real cax = x[N[2]]-x[N[0]];
      tk::real cay = y[N[2]]-y[N[0]];
      tk::real caz = z[N[2]]-z[N[0]];
      tk::real dax = x[N[3]]-x[N[0]];
      tk::real day = y[N[3]]-y[N[0]];
      tk::real daz = z[N[3]]-z[N[0]];
      auto J = tk::triple( bax, bay, baz, cax, cay, caz, dax, day, daz );
      auto J24 = J/24.0;
      // shape function derivatives, nnode*ndim [4][3]
      tk::real g[4][3];
      tk::crossdiv( cax, cay, caz, dax, day, daz, J,
                    g[1][0], g[1][1], g[1][2] );
      tk::crossdiv( dax, day, daz, bax, bay, baz, J,
                    g[2][0], g[2][1], g[2][2] );
      tk::crossdiv( bax, bay, baz, cax, cay, caz, J,
                    g[3][0], g[3][1], g[3][2] );
      for (std::size_t i=0; i<3; ++i)
        g[0][i] = -g[1][i] - g[2][i] - g[3][i];
      // scatter-add gradient contributions to boundary nodes
      tk::real u[5];
      for (std::size_t b=0; b<4; ++b) {
        u[0] = U(N[b],0,0);
        u[1] = U(N[b],1,0)/u[0];
        u[2] = U(N[b],2,0)/u[0];
        u[3] = U(N[b],3,0)/u[0];
        u[4] = U(N[b],4,0)/u[0] - 0.5*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3]);
        if ( !skipPoint(x[N[b]],y[N[b]],z[N[b]]) &&
             stagPoint(x[N[b]],y[N[b]],z[N[b]]) )
        {
          u[1] = u[2] = u[3] = 0.0;
        }
        for (std::size_t c=0; c<5; ++c)
          for (std::size_t i=0; i<3; ++i)
            G(p,c*3+i,0) += J24 * g[b][i] * u[c];
      }
    }

  // put in nodal gradients of chare-boundary points
  for (const auto& [g,b] : bid) {
    auto i = tk::cref_find( lid, g );
    for (std::size_t c=0; c<G.nprop(); ++c)
      G(i,c,0) = bG(b,c,0);
  }

  // divide weak result in gradients by nodal volume
  for (std::size_t p=0; p<npoin; ++p)
    for (std::size_t c=0; c<5*3; ++c)
      G(p,c,0) /= vol[p];

  return G;
}

static std::vector< tk::real >
spongePressures( const std::array< std::vector< tk::real >, 3 >& coord,
                 const std::unordered_set< std::size_t >& nodes )
// *****************************************************************************
//! Compute sponge pressure multiplers at sponge side sets
//! \param[in] coord Mesh node coordinates
//! \param[in] nodes Unique set of nodes for sponge conditions
//! \return Sponge ressure multiplers at nodes, one per sponge side set
//! \details This function computes a sponge-like multiplier that will be
//!   applied to nodes of side sets specified in the input file. This is
//!   used to reduce the pressure gradient normal to boundaries and thereby
//!   modeling the effect of a solid wall on the fluid via fluid-structure
//!   interaction.
//! \note If no sponge pressure coefficients are configured, an empty
//!   vector is returned.
// *****************************************************************************
{
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];
  std::vector< tk::real > spmult;
  std::size_t nset = 0;     // number of sponge side sets configured
  const auto& sponge = g_inputdeck.get< tag::param, tag::compflow, tag::sponge >();
  const auto& ss = sponge.get< tag::sideset >();
  if (ss.size() > 0) {  // if symbcs configured for this system
    const auto& sppre = sponge.get< tag::pressure >();
    nset = ss[0].size();  // number of sponge side sets configured
    spmult.resize( x.size() * nset, 0.0 );
    for (auto p : nodes) {
      if (not skipPoint(x[p],y[p],z[p]) && sppre.size() > 0) {
        Assert( nset == sppre[0].size(), "Size mismatch" );
        for (std::size_t s=0; s<nset; ++s)
          spmult[p*nset+s] = sppre[0][s];
      } else {
        for (std::size_t s=0; s<nset; ++s)
          spmult[p*nset+s] = 0.0;
      }
    }
  }
  Assert( ss.size() > 0 ?
          spmult.size() == x.size() * ss[0].size() :
          spmult.size() == 0, "Sponge pressure multipler wrong size" );
  return spmult;
}

static void
muscl( std::size_t p,
       std::size_t q,
       const tk::UnsMesh::Coords& coord,
       const tk::Fields& Grad,
       tk::real& rL, tk::real& uL, tk::real& vL, tk::real& wL, tk::real& eL,
       tk::real& rR, tk::real& uR, tk::real& vR, tk::real& wR, tk::real& eR )
// *****************************************************************************
//   Compute MUSCL reconstruction in edge-end points using a MUSCL procedure
//! \param[in] p Left node id of edge-end
//! \param[in] q Right node id of edge-end
//! \param[in] coord Array of nodal coordinates
//! \param[in] Grad Gradient of all unknowns in mesh points
//! \param[in,out] rL Left density
//! \param[in,out] uL Left X velocity
//! \param[in,out] vL Left Y velocity
//! \param[in,out] wL Left Z velocity
//! \param[in,out] eL Left internal energy
//! \param[in,out] rR Right density
//! \param[in,out] uR Right X velocity
//! \param[in,out] vR Right Y velocity
//! \param[in,out] wR Right Z velocity
//! \param[in,out] eR Right internal energy
// *****************************************************************************
{
  const tk::real muscl_eps = 1.0e-9;
  const tk::real muscl_const = 1.0/3.0;

  // access node coordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // edge vector
  std::array< tk::real, 3 > vw{ x[q]-x[p], y[q]-y[p], z[q]-z[p] };

  tk::real delta1[5], delta2[5], delta3[5];
  std::array< tk::real, 5 > ls{ rL, uL, vL, wL, eL };
  std::array< tk::real, 5 > rs{ rR, uR, vR, wR, eR };
  auto url = ls;
  auto urr = rs;

  // MUSCL reconstruction of edge-end-point primitive variables
  for (std::size_t c=0; c<5; ++c) {
    // gradients
    std::array< tk::real, 3 >
      g1{ Grad(p,c*3+0,0), Grad(p,c*3+1,0), Grad(p,c*3+2,0) },
      g2{ Grad(q,c*3+0,0), Grad(q,c*3+1,0), Grad(q,c*3+2,0) };

    delta2[c] = rs[c] - ls[c];
    delta1[c] = 2.0 * tk::dot(g1,vw) - delta2[c];
    delta3[c] = 2.0 * tk::dot(g2,vw) - delta2[c];

    // MUSCL extrapolation option 1:
    // ---------------------------------------------------------------------
    // Uncomment the following 3 blocks of code if this version is required.
    // this reconstruction is from the following paper:
    // Waltz, J., Morgan, N. R., Canfield, T. R., Charest, M. R.,
    // Risinger, L. D., & Wohlbier, J. G. (2014). A three-dimensional
    // finite element arbitrary Lagrangian–Eulerian method for shock
    // hydrodynamics on unstructured grids. Computers & Fluids, 92,
    // 172-187.

    //// form limiters
    //auto rcL = (delta2[c] + muscl_eps) / (delta1[c] + muscl_eps);
    //auto rcR = (delta2[c] + muscl_eps) / (delta3[c] + muscl_eps);
    //auto rLinv = (delta1[c] + muscl_eps) / (delta2[c] + muscl_eps);
    //auto rRinv = (delta3[c] + muscl_eps) / (delta2[c] + muscl_eps);

    //// van Leer limiter
    //// any other symmetric limiter could be used instead too
    //auto phiL = (std::abs(rcL) + rcL) / (std::abs(rcL) + 1.0);
    //auto phiR = (std::abs(rcR) + rcR) / (std::abs(rcR) + 1.0);
    //auto phi_L_inv = (std::abs(rLinv) + rLinv) / (std::abs(rLinv) + 1.0);
    //auto phi_R_inv = (std::abs(rRinv) + rRinv) / (std::abs(rRinv) + 1.0);

    //// update unknowns with reconstructed unknowns
    //const tk::real muscl_m1 = 1.0 - muscl_const;
    //const tk::real muscl_p1 = 1.0 + muscl_const;
    //url[c] += 0.25*(delta1[c]*muscl_m1*phiL + delta2[c]*muscl_p1*phi_L_inv);
    //urr[c] -= 0.25*(delta3[c]*muscl_m1*phiR + delta2[c]*muscl_p1*phi_R_inv);

    // ---------------------------------------------------------------------

    // MUSCL extrapolation option 2:
    // ---------------------------------------------------------------------
    // The following 2 blocks of code.
    // this reconstruction is from the following paper:
    // Luo, H., Baum, J. D., & Lohner, R. (1994). Edge-based finite element
    // scheme for the Euler equations. AIAA journal, 32(6), 1183-1190.
    // Van Leer, B. (1974). Towards the ultimate conservative difference
    // scheme. II. Monotonicity and conservation combined in a second-order
    // scheme. Journal of computational physics, 14(4), 361-370.

    // get Van Albada limiter
    // the following form is derived from the flux limiter phi as:
    // s = phi_inv - (1 - phi)
    auto sL = std::max(0.0, (2.0*delta1[c]*delta2[c] + muscl_eps)
      /(delta1[c]*delta1[c] + delta2[c]*delta2[c] + muscl_eps));
    auto sR = std::max(0.0, (2.0*delta3[c]*delta2[c] + muscl_eps)
      /(delta3[c]*delta3[c] + delta2[c]*delta2[c] + muscl_eps));

    // update unknowns with reconstructed unknowns
    url[c] += 0.25*sL*(delta1[c]*(1.0-muscl_const*sL)
      + delta2[c]*(1.0+muscl_const*sL));
    urr[c] -= 0.25*sR*(delta3[c]*(1.0-muscl_const*sR)
      + delta2[c]*(1.0+muscl_const*sR));

    // ---------------------------------------------------------------------
  }

  // force first order if the reconstructions for density or internal energy
  // would have allowed negative values
  if (ls[0] < delta1[0] || ls[4] < delta1[4]) url = ls;
  if (rs[0] < -delta3[0] || rs[4] < -delta3[4]) urr = rs;

  rL = url[0];
  uL = url[1];
  vL = url[2];
  wL = url[3];
  eL = url[4];

  rR = urr[0];
  uR = urr[1];
  vR = urr[2];
  wR = urr[3];
  eR = urr[4];
}

#pragma omp declare simd
static void
flux( tk::real nx, tk::real ny, tk::real nz,
      tk::real mx, tk::real my, tk::real mz,
      tk::real rL, tk::real ruL, tk::real rvL, tk::real rwL, tk::real reL,
      tk::real rR, tk::real ruR, tk::real rvR, tk::real rwR, tk::real reR,
      tk::real pL, tk::real pR,
      tk::real& fr, tk::real& fru, tk::real& frv, tk::real& frw, tk::real& fre )
// *****************************************************************************
//! Rusanov approximate Riemann solver flux function
//! \param[in] nx X component of the surface normal
//! \param[in] ny Y component of the surface normal
//! \param[in] nz Z component of the surface normal
//! \param[in] mx X component of the weighted surface normal on chare
//!   boundary, weighted by the number of contributions to the edge
//! \param[in] my Y component of the weighted surface normal on chare
//!   boundary, weighted by the number of contributions to the edge
//! \param[in] mz Z component of the weighted surface normal on chare
//!   boundary, weighted by the number of contributions to the edge
//! \param[in] rL Left density
//! \param[in] ruL Left X momentum
//! \param[in] rvL Left Y momentum
//! \param[in] rwL Left Z momentum
//! \param[in] reL Left total specific energy
//! \param[in] rR Right density
//! \param[in] ruR Right X momentum
//! \param[in] rvR Right Y momentum
//! \param[in] rwR Right Z momentum
//! \param[in] reR Right total specific energy
//! \param[in] pL Left pressure
//! \param[in] pR Right pressure
//! \param[in,out] fr Riemann solution for density according
//! \param[in,out] fru Riemann solution for X momenutm according
//! \param[in,out] frv Riemann solution for Y momenutm according
//! \param[in,out] frw Riemann solution for Z momenutm according
//! \param[in,out] fre Riemann solution for specific total energy
// *****************************************************************************
{
  auto ul = ruL/rL;
  auto vl = rvL/rL;
  auto wl = rwL/rL;

  auto ur = ruR/rR;
  auto vr = rvR/rR;
  auto wr = rwR/rR;

  auto al = eos_soundspeed( rL, pL );
  auto ar = eos_soundspeed( rR, pR );

  // dissipation
  tk::real len = tk::length( {mx,my,mz} );
  tk::real vml = ul*mx + vl*my + wl*mz;
  tk::real vmr = ur*mx + vr*my + wr*mz;
  auto sl = std::abs(vml) + al*len;
  auto sr = std::abs(vmr) + ar*len;
  auto smax = std::max( sl, sr );

  // face-normal velocities
  tk::real vnl = ul*nx + vl*ny + wl*nz;
  tk::real vnr = ur*nx + vr*ny + wr*nz;

  // numerical fluxes
  fr  = 0.5*(rL*vnl + rR*vnr - smax*(rR - rL));
  fru = 0.5*(ruL*vnl + pL*nx + ruR*vnr + pR*nx - smax*(ruR - ruL));
  frv = 0.5*(rvL*vnl + pL*ny + rvR*vnr + pR*ny - smax*(rvR - rvL));
  frw = 0.5*(rwL*vnl + pL*nz + rwR*vnr + pR*nz - smax*(rwR - rwL));
  fre = 0.5*(reL*vnl + reR*vnr
             + pL*(ruL*nx + rvL*ny + rwL*nz)/rL
             + pR*(ruR*nx + rvR*ny + rwR*nz)/rR
             - smax*(reR - reL));
}

static void
domainint( const std::array< std::vector< tk::real >, 3 >& coord,
           const std::vector< std::size_t >& gid,
           const std::vector< std::size_t >& edgenode,
           const std::vector< std::size_t >& edgeid,
           const std::pair< std::vector< std::size_t >,
                            std::vector< std::size_t > >& psup,
           const std::vector< tk::real >& dfn,
           const tk::Fields& U,
           const tk::Fields& Grad,
           const std::vector< tk::real >& spmult,
           tk::Fields& R )
// *****************************************************************************
//! Compute domain-edge integral
//! \param[in] coord Mesh node coordinates
//! \param[in] gid Local->glocal node ids
//! \param[in] edgenode Local node ids of edges
//! \param[in] edgeid Local node id pair -> edge id map
//! \param[in] psup Points surrounding points
//! \param[in] dfn Dual-face normals
//! \param[in] U Solution vector at recent time step
//! \param[in] Grad Nodal gradients
//! \param[in] spmult Sponge pressure multiplers at nodes, one per symBC set
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  // domain-edge integral: compute fluxes in edges
  std::vector< tk::real > dflux( edgenode.size()/2 * 5 );

  // access node coordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // number of side sets configured with sponge pressure multipliers
  std::size_t nset = spmult.size() / x.size();

  #pragma omp simd
  for (std::size_t e=0; e<edgenode.size()/2; ++e) {
    auto p = edgenode[e*2+0];
    auto q = edgenode[e*2+1];

    // compute primitive variables at edge-end points
    tk::real rL  = U(p,0,0);
    tk::real ruL = U(p,1,0) / rL;
    tk::real rvL = U(p,2,0) / rL;
    tk::real rwL = U(p,3,0) / rL;
    tk::real reL = U(p,4,0) / rL - 0.5*(ruL*ruL + rvL*rvL + rwL*rwL);
    tk::real rR  = U(q,0,0);
    tk::real ruR = U(q,1,0) / rR;
    tk::real rvR = U(q,2,0) / rR;
    tk::real rwR = U(q,3,0) / rR;
    tk::real reR = U(q,4,0) / rR - 0.5*(ruR*ruR + rvR*rvR + rwR*rwR);

    // apply stagnation BCs to primitive variables
    if ( !skipPoint(x[p],y[p],z[p]) && stagPoint(x[p],y[p],z[p]) )
      ruL = rvL = rwL = 0.0;
    if ( !skipPoint(x[q],y[q],z[q]) && stagPoint(x[q],y[q],z[q]) )
      ruR = rvR = rwR = 0.0;

    // compute MUSCL reconstruction in edge-end points
    muscl( p, q, coord, Grad,
           rL, ruL, rvL, rwL, reL, rR, ruR, rvR, rwR, reR );

    // convert back to conserved variables
    reL = (reL + 0.5*(ruL*ruL + rvL*rvL + rwL*rwL)) * rL;
    ruL *= rL;
    rvL *= rL;
    rwL *= rL;
    reR = (reR + 0.5*(ruR*ruR + rvR*rvR + rwR*rwR)) * rR;
    ruR *= rR;
    rvR *= rR;
    rwR *= rR;

    // evaluate pressure at edge-end points
    tk::real pL = eos_pressure( rL, ruL/rL, rvL/rL, rwL/rL, reL );
    tk::real pR = eos_pressure( rR, ruR/rR, rvR/rR, rwR/rR, reR );

    // apply sponge-pressure multipliers
    for (std::size_t s=0; s<nset; ++s) {
      pL -= pL*spmult[p*nset+s];
      pR -= pR*spmult[q*nset+s];
    }

    // compute Riemann flux using edge-end point states
    tk::real f[5];
    flux( dfn[e*6+0], dfn[e*6+1], dfn[e*6+2],
          dfn[e*6+3], dfn[e*6+4], dfn[e*6+5],
          rL, ruL, rvL, rwL, reL,
          rR, ruR, rvR, rwR, reR,
          pL, pR,
          f[0], f[1], f[2], f[3], f[4] );
    // store flux in edges
    for (std::size_t c=0; c<5; ++c) dflux[e*5+c] = f[c];
  }

  // access pointer to right hand side at component and offset
  std::array< const tk::real*, 5 > r;
  for (std::size_t c=0; c<5; ++c) r[c] = R.cptr( c, 0 );

  // domain-edge integral: sum flux contributions to points
  for (std::size_t p=0,k=0; p<U.nunk(); ++p)
    for (auto q : tk::Around(psup,p)) {
      auto s = gid[p] > gid[q] ? -1.0 : 1.0;
      auto e = edgeid[k++];
      // the 2.0 in the following expression is so that the RHS contribution
      // conforms with Eq 12 (Waltz et al. Computers & fluids (92) 2014);
      // The 1/2 in Eq 12 is extracted from the flux function (Rusanov).
      // However, Rusanov::flux computes the flux with the 1/2. This 2
      // cancels with the 1/2 in Rusanov::flux, so that the 1/2 can be
      // extracted out and multiplied as in Eq 12
      for (std::size_t c=0; c<5; ++c)
        R.var(r[c],p) -= 2.0*s*dflux[e*5+c];
    }

  tk::destroy(dflux);
}

static void
bndint( const std::array< std::vector< tk::real >, 3 >& coord,
        const std::vector< std::size_t >& triinpoel,
        const std::vector< int >& symbctri,
        const tk::Fields& U,
        const std::vector< tk::real >& spmult,
        tk::Fields& R )
// *****************************************************************************
//! Compute boundary integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] triinpoel Boundary triangle face connecitivity with local ids
//! \param[in] symbctri Vector with 1 at symmetry BC boundary triangles
//! \param[in] U Solution vector at recent time step
//! \param[in] spmult Sponge pressure multiplers at nodes, one per symBC set
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  // access node coordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // boundary integrals: compute fluxes in edges
  std::vector< tk::real > bflux( triinpoel.size() * 5 * 2 );

  // number of side sets configured with sponge pressure multipliers
  std::size_t nset = spmult.size() / x.size();

  #pragma omp simd
  for (std::size_t e=0; e<triinpoel.size()/3; ++e) {
    // access node IDs
    std::size_t N[3] =
      { triinpoel[e*3+0], triinpoel[e*3+1], triinpoel[e*3+2] };
    // access solution at element nodes
    tk::real rA  = U(N[0],0,0);
    tk::real rB  = U(N[1],0,0);
    tk::real rC  = U(N[2],0,0);
    tk::real ruA = U(N[0],1,0);
    tk::real ruB = U(N[1],1,0);
    tk::real ruC = U(N[2],1,0);
    tk::real rvA = U(N[0],2,0);
    tk::real rvB = U(N[1],2,0);
    tk::real rvC = U(N[2],2,0);
    tk::real rwA = U(N[0],3,0);
    tk::real rwB = U(N[1],3,0);
    tk::real rwC = U(N[2],3,0);
    tk::real reA = U(N[0],4,0);
    tk::real reB = U(N[1],4,0);
    tk::real reC = U(N[2],4,0);
    // apply stagnation BCs
    if ( !skipPoint(x[N[0]],y[N[0]],z[N[0]]) &&
         stagPoint(x[N[0]],y[N[0]],z[N[0]]) )
    {
      ruA = rvA = rwA = 0.0;
    }
    if ( !skipPoint(x[N[1]],y[N[1]],z[N[1]]) &&
         stagPoint(x[N[1]],y[N[1]],z[N[1]]) )
    {
      ruB = rvB = rwB = 0.0;
    }
    if ( !skipPoint(x[N[2]],y[N[2]],z[N[2]]) &&
         stagPoint(x[N[2]],y[N[2]],z[N[2]]) )
    {
      ruC = rvC = rwC = 0.0;
    }
    // compute face normal
    tk::real nx, ny, nz;
    tk::normal( x[N[0]], x[N[1]], x[N[2]],
                y[N[0]], y[N[1]], y[N[2]],
                z[N[0]], z[N[1]], z[N[2]],
                nx, ny, nz );
    // compute boundary flux
    tk::real f[5][3];
    tk::real p, vn;
    int sym = symbctri[e];
    p = eos_pressure( rA, ruA/rA, rvA/rA, rwA/rA, reA );
    for (std::size_t s=0; s<nset; ++s) p -= p*spmult[N[0]*nset+s];
    vn = sym ? 0.0 : (nx*(ruA/rA) + ny*(rvA/rA) + nz*(rwA/rA));
    f[0][0] = rA*vn;
    f[1][0] = ruA*vn + p*nx;
    f[2][0] = rvA*vn + p*ny;
    f[3][0] = rwA*vn + p*nz;
    f[4][0] = reA*vn + p*(sym ? 0.0 : (nx*ruA + ny*rvA + nz*rwA)/rA);
    p = eos_pressure( rB, ruB/rB, rvB/rB, rwB/rB, reB );
    for (std::size_t s=0; s<nset; ++s) p -= p*spmult[N[1]*nset+s];
    vn = sym ? 0.0 : (nx*(ruB/rB) + ny*(rvB/rB) + nz*(rwB/rB));
    f[0][1] = rB*vn;
    f[1][1] = ruB*vn + p*nx;
    f[2][1] = rvB*vn + p*ny;
    f[3][1] = rwB*vn + p*nz;
    f[4][1] = reB*vn + p*(sym ? 0.0 : (nx*ruB + ny*rvB + nz*rwB)/rB);
    p = eos_pressure( rC, ruC/rC, rvC/rC, rwC/rC, reC );
    for (std::size_t s=0; s<nset; ++s) p -= p*spmult[N[2]*nset+s];
    vn = sym ? 0.0 : (nx*(ruC/rC) + ny*(rvC/rC) + nz*(rwC/rC));
    f[0][2] = rC*vn;
    f[1][2] = ruC*vn + p*nx;
    f[2][2] = rvC*vn + p*ny;
    f[3][2] = rwC*vn + p*nz;
    f[4][2] = reC*vn + p*(sym ? 0.0 : (nx*ruC + ny*rvC + nz*rwC)/rC);
    // compute face area
    auto A6 = tk::area( x[N[0]], x[N[1]], x[N[2]],
                        y[N[0]], y[N[1]], y[N[2]],
                        z[N[0]], z[N[1]], z[N[2]] ) / 6.0;
    auto A24 = A6/4.0;
    // store flux in boundary elements
    for (std::size_t c=0; c<5; ++c) {
      auto eb = (e*5+c)*6;
      auto Bab = A24 * (f[c][0] + f[c][1]);
      bflux[eb+0] = Bab + A6 * f[c][0];
      bflux[eb+1] = Bab;
      Bab = A24 * (f[c][1] + f[c][2]);
      bflux[eb+2] = Bab + A6 * f[c][1];
      bflux[eb+3] = Bab;
      Bab = A24 * (f[c][2] + f[c][0]);
      bflux[eb+4] = Bab + A6 * f[c][2];
      bflux[eb+5] = Bab;
    }
  }

  // access pointer to right hand side at component and offset
  std::array< const tk::real*, 5 > r;
  for (std::size_t c=0; c<5; ++c) r[c] = R.cptr( c, 0 );

  // boundary integrals: sum flux contributions to points
  for (std::size_t e=0; e<triinpoel.size()/3; ++e)
    for (std::size_t c=0; c<5; ++c) {
      auto eb = (e*5+c)*6;
      R.var(r[c],triinpoel[e*3+0]) -= bflux[eb+0] + bflux[eb+5];
      R.var(r[c],triinpoel[e*3+1]) -= bflux[eb+1] + bflux[eb+2];
      R.var(r[c],triinpoel[e*3+2]) -= bflux[eb+3] + bflux[eb+4];
    }

  tk::destroy(bflux);
}

static void
src( const std::array< std::vector< tk::real >, 3 >& coord,
     const std::vector< std::size_t >& inpoel,
     tk::real t,
     const std::vector< tk::real >& tp,
     tk::Fields& R )
// *****************************************************************************
//! Compute optional source integral
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] t Physical time
//! \param[in] tp Physical time for each mesh node
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  auto src = SRC();
  if (!src) return;

  // access node coordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // access pointer to right hand side at component and offset
  std::array< const tk::real*, 5 > r;
  for (std::size_t c=0; c<5; ++c) r[c] = R.cptr( c, 0 );

  // source integral
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    std::size_t N[4] =
      { inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] };
    // compute element Jacobi determinant, J = 6V
    auto J24 = tk::triple(
      x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]],
      x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]],
      x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] ) / 24.0;
    // sum source contributions to nodes
    for (std::size_t a=0; a<4; ++a) {
      tk::real s[5];
      if (g_inputdeck.get< tag::discr, tag::steady_state >()) t = tp[N[a]];
      src( x[N[a]], y[N[a]], z[N[a]], t, s[0], s[1], s[2], s[3], s[4] );
      for (std::size_t c=0; c<5; ++c)
        R.var(r[c],N[a]) += J24 * s[c];
    }
  }
}

void
rhs( tk::real t,
     const std::array< std::vector< tk::real >, 3 >& coord,
     const std::vector< std::size_t >& inpoel,
     const std::vector< std::size_t >& triinpoel,
     const std::vector< std::size_t >& gid,
     const std::unordered_map< std::size_t, std::size_t >& bid,
     const std::unordered_map< std::size_t, std::size_t >& lid,
     const std::vector< tk::real >& dfn,
     const std::pair< std::vector< std::size_t >,
                      std::vector< std::size_t > >& esup,
     const std::pair< std::vector< std::size_t >,
                      std::vector< std::size_t > >& psup,
     const std::vector< int >& symbctri,
     const std::unordered_set< std::size_t >& spongenodes,
     const std::vector< tk::real >& vol,
     const std::vector< std::size_t >& edgenode,
     const std::vector< std::size_t >& edgeid,
     const std::vector< tk::real >& tp,
     const tk::Fields& bG,
     const tk::Fields& U,
     tk::Fields& R )
// *****************************************************************************
//! Compute right hand side
//! \param[in] t Physical time
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] triinpoel Boundary triangle face connecitivity with local ids
//! \param[in] bid Local chare-boundary node ids (value) associated to
//!    global node ids (key)
//! \param[in] gid Local->glocal node ids
//! \param[in] lid Global->local node ids
//! \param[in] dfn Dual-face normals
//! \param[in] esup Elements surrounding points
//! \param[in] psup Points surrounding points
//! \param[in] symbctri Vector with 1 at symmetry BC boundary triangles
//! \param[in] spongenodes Unique set of nodes at which to apply sponge
//               conditions
//! \param[in] vol Nodal volumes
//! \param[in] edgenode Local node IDs of edges
//! \param[in] edgeid Edge ids in the order of access
//! \param[in] bG Nodal gradients in chare boundaries
//! \param[in] U Solution vector at recent time step
//! \param[in] tp Physical time for each mesh node
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
          "vector at recent time step incorrect" );
  Assert( R.nunk() == coord[0].size(),
          "Number of unknowns and/or number of components in right-hand "
          "side vector incorrect" );

  // compute/assemble gradients in points
  auto G = grad( coord, inpoel, esup, lid, bid, vol, U, bG );

  // zero right hand side for all components
  R.fill( 0.0 );

  // compute sponge pressure multiplers at sponge side sets
  auto spmult = spongePressures( coord, spongenodes );

  // compute domain-edge integral
  domainint( coord, gid, edgenode, edgeid, psup, dfn, U, G, spmult, R );

  // compute boundary integrals
  bndint( coord, triinpoel, symbctri, U, spmult, R );

  // compute optional source integral
  src( coord, inpoel, t, tp, R );
}

tk::real
dt( const std::array< std::vector< tk::real >, 3 >& coord,
    const std::vector< std::size_t >& inpoel,
    tk::real t,
    const tk::Fields& U )
// *****************************************************************************
//! Compute the minimum time step size (for unsteady time stepping)
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] t Physical time
//! \param[in] U Solution vector at recent time step
//! \return Minimum time step size
// *****************************************************************************
{
  Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
          "vector at recent time step incorrect" );

  // energy source propagation time and velocity
  const auto& ic = g_inputdeck.get< tag::param, tag::compflow, tag::ic >();
  const auto& icbox = ic.get< tag::box >();

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // ratio of specific heats
  //auto g = g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[0][0];
  // compute the minimum dt across all elements we own
  tk::real mindt = std::numeric_limits< tk::real >::max();
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                           inpoel[e*4+2], inpoel[e*4+3] }};
    // compute cubic root of element volume as the characteristic length
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto L = std::cbrt( tk::triple( ba, ca, da ) / 6.0 );
    // access solution at element nodes at recent time step
    std::array< std::array< tk::real, 4 >, 5 > u;
    for (std::size_t c=0; c<5; ++c) u[c] = U.extract( c, 0, N );
    // compute the maximum length of the characteristic velocity (fluid
    // velocity + sound velocity) across the four element nodes
    tk::real maxvel = 0.0;
    for (std::size_t j=0; j<4; ++j) {
      auto& r  = u[0][j];    // rho
      auto& ru = u[1][j];    // rho * u
      auto& rv = u[2][j];    // rho * v
      auto& rw = u[3][j];    // rho * w
      auto& re = u[4][j];    // rho * e
      auto p = eos_pressure( r, ru/r, rv/r, rw/r, re );
      if (p < 0) p = 0.0;
      auto c = eos_soundspeed( r, p );
      auto v = std::sqrt((ru*ru + rv*rv + rw*rw)/r/r) + c; // char. velocity

      // energy source propagation velocity (in all IC boxes configured)
      if (icbox.size() > 0) {
        for (const auto& b : icbox[0]) {
          const auto& initiate = b.get< tag::initiate >();
          auto iv = initiate.get< tag::velocity >();
          auto inittype = initiate.get< tag::init >();
          if (inittype == ctr::InitiateType::LINEAR) {
            auto zmin = b.get< tag::zmin >();
            auto zmax = b.get< tag::zmax >();
            auto wFront = 0.08;
            auto tInit = 0.0;
            auto tFinal = tInit + (zmax - zmin - 2.0*wFront) /
              std::fabs(iv);
            if (t >= tInit && t <= tFinal)
              v = std::max(v, std::fabs(iv));
          }
        }
      }

      if (v > maxvel) maxvel = v;
    }
    // compute element dt for the Euler equations
    auto euler_dt = L / std::max( maxvel, 1.0e-8 );
    // compute element dt based on the viscous force
    //auto viscous_dt = m_physics.viscous_dt( L, u );
    // compute element dt based on thermal diffusion
    //auto conduct_dt = m_physics.conduct_dt( L, g, u );
    // compute minimum element dt
    auto elemdt = euler_dt;//std::min( euler_dt, std::min( viscous_dt, conduct_dt ) );
    // find minimum dt across all elements
    mindt = std::min( elemdt, mindt );
  }
  mindt *= g_inputdeck.get< tag::discr, tag::cfl >();

  return mindt;
}

void
dt( uint64_t,
    const std::vector< tk::real >& vol,
    const tk::Fields& U,
    std::vector< tk::real >& dtp )
// *****************************************************************************
//! Compute a time step size for each mesh node (for steady time stepping)
//! \param[in] U Solution vector at recent time step
//! \param[in] vol Nodal volume (with contributions from other chares)
//! \param[in,out] dtp Time step size for each mesh node
// *****************************************************************************
{
  for (std::size_t i=0; i<U.nunk(); ++i) {
    // compute cubic root of element volume as the characteristic length
    const auto L = std::cbrt( vol[i] );
    // access solution at node p at recent time step
    const auto u = U[i];
    // compute pressure
    auto p = eos_pressure( u[0], u[1]/u[0], u[2]/u[0], u[3]/u[0], u[4] );
    if (p < 0) p = 0.0;
    auto c = eos_soundspeed( u[0], p );
    // characteristic velocity
    auto v = std::sqrt((u[1]*u[1] + u[2]*u[2] + u[3]*u[3])/u[0]/u[0]) + c;
    // compute dt for node
    dtp[i] = L / v * g_inputdeck.get< tag::discr, tag::cfl >();
  }
}

std::map< std::size_t, std::vector< std::pair< bool, tk::real > > >
dirbc( tk::real t,
       tk::real deltat,
       const std::vector< tk::real >& tp,
       const std::vector< tk::real >& dtp,
       const std::pair< const int, std::vector< std::size_t > >& ss,
       const std::array< std::vector< tk::real >, 3 >& coord )
// *****************************************************************************
//! \brief Query Dirichlet boundary condition value on a given side set for
//!    all components in this PDE system
//! \param[in] t Physical time
//! \param[in] deltat Time step size
//! \param[in] tp Physical time for each mesh node
//! \param[in] dtp Time step size for each mesh node
//! \param[in] ss Pair of side set ID and (local) node IDs on the side set
//! \param[in] coord Mesh node coordinates
//! \return Vector of pairs of bool and boundary condition value associated
//!   to mesh node IDs at which Dirichlet boundary conditions are set.
// *****************************************************************************
{
  using tag::param; using tag::bcdir;
  using NodeBC = std::vector< std::pair< bool, tk::real > >;
  std::map< std::size_t, NodeBC > bc;
  const auto steady = g_inputdeck.get< tag::discr, tag::steady_state >();
  const auto& ubc =
    g_inputdeck.get< tag::param, tag::compflow, tag::bc, bcdir >();
  if (!ubc.empty()) {
    auto ic = IC();
    Assert( ubc.size() > 0, "Indexing out of Dirichlet BC vector" );
    const auto& x = coord[0];
    const auto& y = coord[1];
    const auto& z = coord[2];
    for (const auto& b : ubc[0])
      if (std::stoi(b) == ss.first)
        for (auto n : ss.second) {
          Assert( x.size() > n, "Indexing out of coordinate array" );
          if (steady) { t = tp[n]; deltat = dtp[n]; }
          auto s = ic( x[n], y[n], z[n], t+deltat );
          if ( !skipPoint(x[n],y[n],z[n]) && stagPoint(x[n],y[n],z[n]) ) {
            s[1] = s[2] = s[3] = 0.0;
          }
          bc[n] = {{ {true,s[0]}, {true,s[1]}, {true,s[2]}, {true,s[3]},
                     {true,s[4]} }};
        }
  }
  return bc;
}

void
symbc( tk::Fields& U,
       const std::array< std::vector< tk::real >, 3 >& coord,
       const std::unordered_map< int,
         std::unordered_map< std::size_t, std::array< tk::real, 4 > > >& bnorm,
       const std::unordered_set< std::size_t >& nodes )
// *****************************************************************************
//! Set symmetry boundary conditions at nodes
//! \param[in] U Solution vector at recent time step
//! \param[in] coord Mesh node coordinates
//! \param[in] bnorm Face normals in boundary points, key local node id,
//!   first 3 reals of value: unit normal, outer key: side set id
//! \param[in] nodes Unique set of node ids at which to set symmetry BCs
// *****************************************************************************
{
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];
  const auto& sbc =
    g_inputdeck.get< tag::param, tag::compflow, tag::bc, tag::bcsym >();
  if (sbc.size() > 0) {             // use symbcs for this system
    for (auto p : nodes) {                 // for all symbc nodes
      if (!skipPoint(x[p],y[p],z[p])) {
        // for all user-def symbc sets
        for (std::size_t s=0; s<sbc[0].size(); ++s) {
          // find nodes & normals for side
          auto j = bnorm.find(std::stoi(sbc[0][s]));
          if (j != end(bnorm)) {
            auto i = j->second.find(p);      // find normal for node
            if (i != end(j->second)) {
              std::array< tk::real, 3 >
                n{ i->second[0], i->second[1], i->second[2] },
                v{ U(p,1,0), U(p,2,0), U(p,3,0) };
              auto v_dot_n = tk::dot( v, n );
              // symbc: remove normal component of velocity
              U(p,1,0) -= v_dot_n * n[0];
              U(p,2,0) -= v_dot_n * n[1];
              U(p,3,0) -= v_dot_n * n[2];
            }
          }
        }
      }
    }
  }
}

void
farfieldbc( tk::Fields& U,
            const std::array< std::vector< tk::real >, 3 >& coord,
            const std::unordered_map< int,
              std::unordered_map< std::size_t, std::array<tk::real,4> > >& bnorm,
            const std::unordered_set< std::size_t >& nodes )
// *****************************************************************************
//! Set farfield boundary conditions at nodes
//! \param[in] U Solution vector at recent time step
//! \param[in] coord Mesh node coordinates
//! \param[in] bnorm Face normals in boundary points, key local node id,
//!   first 3 reals of value: unit normal, outer key: side set id
//! \param[in] nodes Unique set of node ids at which to set farfield BCs
// *****************************************************************************
{
  const auto& fbc =
    g_inputdeck.get< tag::param, tag::compflow, tag::bc, tag::bcfarfield >();

  if (fbc.empty()) return;

  tk::real fr =
    g_inputdeck.get< tag::param, tag::compflow, tag::farfield_density >().size() > 0 ?
    g_inputdeck.get< tag::param, tag::compflow, tag::farfield_density >()[0] : 1.0;

  tk::real fp =
    g_inputdeck.get< tag::param, tag::compflow, tag::farfield_pressure >().size() > 0 ?
    g_inputdeck.get< tag::param, tag::compflow, tag::farfield_pressure >()[0] : 1.0;

  auto fu =
    g_inputdeck.get< tag::param, tag::compflow, tag::farfield_velocity >().size() > 0 ?
    g_inputdeck.get< tag::param, tag::compflow, tag::farfield_velocity >()[0] :
    std::vector< tk::real >( 3, 0.0 );

  if (fbc.size() > 0) {
    const auto& x = coord[0];
    const auto& y = coord[1];
    const auto& z = coord[2];
    for (auto p : nodes)
      if (not skipPoint(x[p],y[p],z[p]))
        for (const auto& s : fbc[0]) {// for all user-def farbc sets
          auto j = bnorm.find(std::stoi(s));// find nodes & normals for side
          if (j != end(bnorm)) {
            auto i = j->second.find(p);      // find normal for node
            if (i != end(j->second)) {
              auto& r  = U(p,0,0);
              auto& ru = U(p,1,0);
              auto& rv = U(p,2,0);
              auto& rw = U(p,3,0);
              auto& re = U(p,4,0);
              auto vn =
                (ru*i->second[0] + rv*i->second[1] + rw*i->second[2]) / r;
              auto a = eos_soundspeed(r, eos_pressure(r, ru/r, rv/r, rw/r, re));
              auto M = vn / a;
              if (M <= -1.0) {                      // supersonic inflow
                r  = fr;
                ru = fr * fu[0];
                rv = fr * fu[1];
                rw = fr * fu[2];
                re = eos_totalenergy( fr, fu[0], fu[1], fu[2], fp );
              } else if (M > -1.0 && M < 0.0) {     // subsonic inflow
                r  = fr;
                ru = fr * fu[0];
                rv = fr * fu[1];
                rw = fr * fu[2];
                re =
                eos_totalenergy( fr, fu[0], fu[1],
                  fu[2], eos_pressure( r, ru/r, rv/r, rw/r, re ) );
              } else if (M >= 0.0 && M < 1.0) {     // subsonic outflow
                re = eos_totalenergy( r, ru/r, rv/r, rw/r, fp );
              }
            }
          }
        }
  }
}

void
sponge( tk::Fields& U,
        const std::array< std::vector< tk::real >, 3 >& coord,
        const std::unordered_set< std::size_t >& nodes )
// *****************************************************************************
//! Apply sponge conditions at sponge nodes
//! \param[in] U Solution vector at recent time step
//! \param[in] coord Mesh node coordinates
//! \param[in] nodes Unique set of node ids at which to apply sponge
//! \details This function applies a sponge-like parameter to nodes of a
//!   side set specified in the input file. We remove a user-specified
//!   percentage of the kinetic energy by reducing the tangential
//!   component of the velocity at a boundary and thereby modeling the
//!   effect of a solid wall on the fluid via fluid-structure interaction
//!   via a viscosity-like effect.
// *****************************************************************************
{
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  const auto& sponge =
    g_inputdeck.get< tag::param, tag::compflow, tag::sponge >();
  const auto& ss = sponge.get< tag::sideset >();

  if (ss.size() > 0) {          // sponge side set for this system
    const auto& spvel = sponge.get< tag::velocity >();
    for (auto p : nodes) {             // for all sponge nodes
      if (!skipPoint(x[p],y[p],z[p])) {
        std::vector< tk::real > sp( ss[0].size(), 0.0 );
        if (spvel.size() > 0) {
          sp = spvel[0];
          for (auto& s : sp) s = std::sqrt(s);
        }
        // sponge velocity: reduce kinetic energy by a user percentage
        for (std::size_t s=0; s<ss[0].size(); ++s) {
          U(p,1,0) -= U(p,1,0)*sp[s];
          U(p,2,0) -= U(p,2,0)*sp[s];
          U(p,3,0) -= U(p,3,0)*sp[s];
        }
      }
    }
  }
}

void
timedepbc( tk::real t,
  tk::Fields& U,
  const std::vector< std::unordered_set< std::size_t > >& nodes,
  const std::vector< tk::Table<5> >& timedepfn )
// *****************************************************************************
//! Apply user defined time dependent BCs
//! \param[in] t Physical time
//! \param[in,out] U Solution vector at recent time step
//! \param[in] nodes Vector of unique sets of node ids at which to apply BCs
//! \param[in] timedepfn Tables discretizing time dependence
//! \details This function applies user defined time dependent boundary
//!   conditions on groups of side sets specified in the input file.
//!   The user specifies pressure, density, and velocity as discrete
//!   functions of time, in the control file, associated with a group of
//!   side sets. Several such groups can be specified, each with their
//!   own discrete function: p(t), rho(t), vx(t), vy(t), vz(t).
// *****************************************************************************
{
  for (std::size_t ib=0; ib<nodes.size(); ++ib) {
    for (auto p:nodes[ib]) {
      // sample primitive vars from discrete data at time t
      auto unk = tk::sample<5>(t, timedepfn[ib]);

      // apply BCs after converting to conserved vars
      U(p,0,0) = unk[1];
      U(p,1,0) = unk[1]*unk[2];
      U(p,2,0) = unk[1]*unk[3];
      U(p,3,0) = unk[1]*unk[4];
      U(p,4,0) = eos_totalenergy( unk[1], unk[2], unk[3], unk[4], unk[0] );
    }
  }
}

} // inciter::
