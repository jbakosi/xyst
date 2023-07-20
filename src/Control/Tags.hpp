// *****************************************************************************
/*!
  \file      src/Control/Tags.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Tags
  \details   Tags are unique types, used for metaprogramming.
*/
// *****************************************************************************
#ifndef Tags_h
#define Tags_h

#include <string>

//! Tags used as unique-type labels for compile-time code-generation
namespace tag {

struct low {};
struct high {};
struct io { static std::string name() { return "io"; } };
struct quiescence { static std::string name() { return "quiescence"; } };
struct trace { static std::string name() { return "trace"; } };
struct version { static std::string name() { return "version"; } };
struct signal { static std::string name() { return "signal"; } };
struct input { static std::string name() { return "input"; } };
struct output { static std::string name() { return "output"; } };
struct restart { static std::string name() { return "restart"; } };
struct nrestart { static std::string name() { return "nrestart"; } };
struct diag { static std::string name() { return "diag"; } };
struct diag_iter { static std::string name() { return "diag_iter"; } };
struct diag_precision { static std::string name() { return "diag_precision"; } };
struct diag_format { static std::string name() { return "diag_format"; } };
struct histout { static std::string name() { return "histout"; } };
struct histout_iter { static std::string name() { return "histout_iter"; } };
struct histout_time { static std::string name() { return "histout_time"; } };
struct histout_range { static std::string name() { return "histout_range"; } };
struct histout_precision { static std::string name() { return "histout_precision"; } };
struct histout_format { static std::string name() { return "histout_format"; } };
struct integout { static std::string name() { return "integout"; } };
struct integout_iter { static std::string name() { return "integout_iter"; } };
struct integout_time { static std::string name() { return "integout_time"; } };
struct integout_range { static std::string name() { return "integout_range"; } };
struct integout_precision { static std::string name() { return "integout_precision"; } };
struct integout_format { static std::string name() { return "integout_format"; } };
struct evalLB {};
struct xminus { static std::string name() { return "xminus"; } };
struct xplus { static std::string name() { return "xplus"; } };
struct yminus { static std::string name() { return "yminus"; } };
struct yplus { static std::string name() { return "yplus"; } };
struct zminus { static std::string name() { return "zminus"; } };
struct zplus { static std::string name() { return "zplus"; } };
struct nonblocking { static std::string name() { return "nonblocking"; } };
struct benchmark { static std::string name() { return "benchmark"; } };
struct lboff {};
struct feedback { static std::string name() { return "feedback"; } };
struct reorder { static std::string name() { return "reorder"; } };
struct steady { static std::string name() { return "steady"; } };
struct residual { static std::string name() { return "residual"; } };
struct error { static std::string name() { return "error"; } };
struct lbfreq { static std::string name() { return "lbfreq"; } };
struct rsfreq { static std::string name() { return "rsfreq"; } };
struct dtfreq { static std::string name() { return "dtfreq"; } };
struct pdf { static std::string name() { return "pdf"; } };
struct nchare {};
struct bounds {};
struct maxit { static std::string name() { return "maxit"; } };
struct tolerance { static std::string name() { return "tolerace"; } };
struct mesh { static std::string name() { return "mesh"; } };
struct couple { static std::string name() { return "couple"; } };
struct filetype { static std::string name() { return "filetype"; } };
struct filename { static std::string name() { return "filename"; } };
struct binsize { static std::string name() { return "binsize"; } };
struct extent { static std::string name() { return "extent"; } };
struct problem { static std::string name() { return "problem"; } };
struct problem_ncomp { static std::string name() { return "problem_ncomp"; } };
struct problem_alpha { static std::string name() { return "problem_alpha"; } };
struct problem_kappa { static std::string name() { return "problem_kappa"; } };
struct problem_beta { static std::string name() { return "problem_beta"; } };
struct problem_r0 { static std::string name() { return "problem_r0"; } };
struct problem_p0 { static std::string name() { return "problem_p0"; } };
struct problem_ce { static std::string name() { return "problem_ce"; } };
struct mat_spec_heat_ratio { static std::string name() { return "mat_spec_heat_ratio"; } };
struct mat_spec_heat_const_vol { static std::string name() { return "mat_spec_heat_const_vol"; } };
struct mat_heat_conductivity { static std::string name() { return "mat_heat_conductivity"; } };
struct pstiff { static std::string name() { return "pstiff"; } };
struct pde { static std::string name() { return "pde"; } };
struct tolref { static std::string name() { return "tolref"; } };
struct href_t0 { static std::string name() { return "href_t0"; } };
struct href_dt { static std::string name() { return "href_dt"; } };
struct href_dtfreq { static std::string name() { return "href_dtfreq"; } };
struct href_maxlevels { static std::string name() { return "href_maxlevels"; } };
struct href_error { static std::string name() { return "href_error"; } };
struct href_init { static std::string name() { return "href_init"; } };
struct tolderef { static std::string name() { return "tolderef"; } };
struct t0ref { static std::string name() { return "t0ref"; } };
struct dtref { static std::string name() { return "dtref"; } };
struct dtref_uniform { static std::string name() { return "dtref_uniform"; } };
struct maxlevels { static std::string name() { return "maxlevels"; } };
struct part { static std::string name() { return "part"; } };
struct partitioned { static std::string name() { return "partitioned"; } };
struct scheme { static std::string name() { return "scheme"; } };
struct nstep { static std::string name() { return "nstep"; } };
struct term { static std::string name() { return "term"; } };
struct t0 { static std::string name() { return "t0"; } };
struct dt { static std::string name() { return "dt"; } };
struct cfl { static std::string name() { return "cfl"; } };
struct sys { static std::string name() { return "sys"; } };
struct refined { static std::string name() { return "refined output"; } };
struct matched {};
struct compatibility {};
struct bndint {};
struct particles { static std::string name() { return "particles"; } };
struct centroid {};
struct nmat { static std::string name() { return "nmat"; } };
struct ttyi { static std::string name() { return "ttyi"; } };
struct dump {};
struct plot {};
struct glob {};
struct control { static std::string name() { return "control"; } };
struct stat { static std::string name() { return "stat"; } };
struct field { static std::string name() { return "field"; } };
struct fieldout { static std::string name() { return "fieldout"; } };
struct fieldout_iter { static std::string name() { return "fieldout_iter"; } };
struct fieldout_time { static std::string name() { return "fieldout_time"; } };
struct fieldout_range { static std::string name() { return "fieldout_range"; } };
struct surface { static std::string name() { return "surface"; } };
struct rho2 { static std::string name() { return "rho2"; } };
struct rho { static std::string name() { return "rho"; } };
struct r { static std::string name() { return "r"; } };
struct c { static std::string name() { return "c"; } };
struct c0 { static std::string name() { return "c0"; } };
struct c1 {};
struct c2 {};
struct c3 { static std::string name() { return "c3"; } };
struct c4 { static std::string name() { return "c4"; } };
struct lambda { static std::string name() { return "lambda"; } };
struct theta { static std::string name() { return "theta"; } };
struct mu { static std::string name() { return "mu"; } };
struct timescale { static std::string name() { return "timescale"; } };
struct depvar { static std::string name() { return "depvar"; } };
struct href_refvar { static std::string name() { return "href_refvar"; } };
struct virtualization {static std::string name() { return "virtualization"; }};
struct normalization { static std::string name() { return "normalization"; } };
struct title { static std::string name() { return "title"; } };
struct selected { static std::string name() { return "selected"; } };
struct component { static std::string name() { return "component"; } };
struct rescomp { static std::string name() { return "residual component"; } };
struct iter { static std::string name() { return "iter"; } };
struct time { static std::string name() { return "time"; } };
struct data { static std::string name() { return "data"; } };
struct range { static std::string name() { return "range"; } };
struct cmd { static std::string name() { return "cmd"; } };
struct param { static std::string name() { return "param"; } };
struct init { static std::string name() { return "init"; } };
struct solve { static std::string name() { return "solve"; } };
struct chare { static std::string name() { return "chare"; } };
struct generator {};
struct help { static std::string name() { return "help"; } };
struct helpctr { static std::string name() { return "helpctr"; } };
struct cmdinfo {};
struct ctrinfo {};
struct group { static std::string name() { return "group"; } };
struct esup {};
struct psup {};
struct gid {};
struct transport { static std::string name() { return "transport"; } };
struct compflow { static std::string name() { return "compflow"; } };
struct physics { static std::string name() { return "physics"; } };
struct diffusivity { static std::string name() { return "diffusivity"; } };
struct u0 { static std::string name() { return "u0"; } };
struct bc_dir { static std::string name() { return "bc_dir"; } };
struct bc_sym { static std::string name() { return "bc_sym"; } };
struct bc_far { static std::string name() { return "bc_far"; } };
struct bc_far_density { static std::string name() { return "bc_far_density"; } };
struct bc_far_pressure { static std::string name() { return "bc_far_pressure"; } };
struct bc_far_velocity { static std::string name() { return "bc_far_velocity"; } };
struct bc_pre { static std::string name() { return "bc_pre"; } };
struct bc_pre_density { static std::string name() { return "bc_pre_density"; } };
struct bc_pre_pressure { static std::string name() { return "bc_pre_pressure"; } };
struct problem_source { static std::string name() { return "problem_source"; } };
struct location { static std::string name() { return "location"; } };
struct radius { static std::string name() { return "radius"; } };
struct release_time { static std::string name() { return "release_time"; } };
struct sideset { static std::string name() { return "sideset"; } };
struct eos { static std::string name() { return "eos"; } };
struct id { static std::string name() { return "id"; } };
struct edge { static std::string name() { return "edge"; } };
struct k { static std::string name() { return "k"; } };
struct com {};
struct queried {};
struct responded {};
struct coord {};
struct refinserted {};
struct discinserted {};
struct disccreated {};
struct workinserted {};
struct distributed {};
struct load {};
struct bcast {};
struct elem {};
struct avecost {};
struct stdcost {};
struct flux { static std::string name() { return "flux"; } };
struct ndof{ static std::string name() { return "ndof"; } };
struct limiter { static std::string name() { return "limiter"; } };
struct cweight { static std::string name() { return "cweight"; } };
struct update {};
struct ch {};
struct pe {};
struct fn { static std::string name() { return "fn"; } };
struct node {};
struct ic { static std::string name() { return "ic"; } };
struct ic_velocity { static std::string name() { return "ic_velocity"; } };
struct ic_density { static std::string name() { return "ic_density"; } };
struct ic_pressure { static std::string name() { return "ic_pressure"; } };
struct ic_energy { static std::string name() { return "ic_energy"; } };
struct ic_temperature { static std::string name() { return "ic_temperature"; } };
struct x { static std::string name() { return "x"; } };
struct y { static std::string name() { return "y"; } };
struct z { static std::string name() { return "z"; } };
struct xmin { static std::string name() { return "xmin"; } };
struct xmax { static std::string name() { return "xmax"; } };
struct ymin { static std::string name() { return "ymin"; } };
struct ymax { static std::string name() { return "ymax"; } };
struct zmin { static std::string name() { return "zmin"; } };
struct zmax { static std::string name() { return "zmax"; } };

} // tag::

#endif // Tags_h
