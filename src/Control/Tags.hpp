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
struct license { static std::string name() { return "license"; } };
struct input { static std::string name() { return "input"; } };
struct output { static std::string name() { return "output"; } };
struct screen { static std::string name() { return "screen"; } };
struct restart { static std::string name() { return "restart"; } };
struct nrestart { static std::string name() { return "nrestart"; } };
struct diag { static std::string name() { return "diag"; } };
struct history { static std::string name() { return "history"; } };
struct evalLB {};
struct xminus { static std::string name() { return "xminus"; } };
struct xplus { static std::string name() { return "xplus"; } };
struct yminus { static std::string name() { return "yminus"; } };
struct yplus { static std::string name() { return "yplus"; } };
struct zminus { static std::string name() { return "zminus"; } };
struct zplus { static std::string name() { return "zplus"; } };
struct verbose { static std::string name() { return "verbose"; } };
struct nonblocking { static std::string name() { return "nonblocking"; } };
struct benchmark { static std::string name() { return "benchmark"; } };
struct lboff {};
struct feedback { static std::string name() { return "feedback"; } };
struct reorder { static std::string name() { return "reorder"; } };
struct pelocal_reorder {
  static std::string name() { return "pelocal_reorder"; } };
struct operator_reorder {
  static std::string name() { return "operator_reorder"; } };
struct steady_state {
  static std::string name() { return "steady_state"; } };
struct residual { static std::string name() { return "residual"; } };
struct error { static std::string name() { return "error"; } };
struct lbfreq { static std::string name() { return "lbfreq"; } };
struct rsfreq { static std::string name() { return "rsfreq"; } };
struct dtfreq { static std::string name() { return "dtfreq"; } };
struct pdf { static std::string name() { return "pdf"; } };
struct ordpdf {};
struct cenpdf {};
struct nchare {};
struct bounds {};
struct meshvelocity { static std::string name() { return "meshvelocity"; } };
struct smoother { static std::string name() { return "smoother"; } };
struct meshforce { static std::string name() { return "meshforce"; } };
struct mesh_motion{ static std::string name() { return "mesh_motion"; } };
struct vortmult { static std::string name() { return "vortmult"; } };
struct maxit { static std::string name() { return "maxit"; } };
struct tolerance { static std::string name() { return "tolerace"; } };
struct mesh { static std::string name() { return "mesh"; } };
struct couple { static std::string name() { return "couple"; } };
struct transfer { static std::string name() { return "transfer"; } };
struct filetype { static std::string name() { return "filetype"; } };
struct filename { static std::string name() { return "filename"; } };
struct location { static std::string name() { return "location"; } };
struct orientation { static std::string name() { return "orientation"; } };
struct reference { static std::string name() { return "reference"; } };
struct pdfpolicy { static std::string name() { return "pdfpolicy"; } };
struct pdfctr { static std::string name() { return "pdfctr"; } };
struct pdfnames { static std::string name() { return "pdfnames"; } };
struct flformat { static std::string name() { return "flformat"; } };
struct prec { static std::string name() { return "precision"; } };
struct binsize { static std::string name() { return "binsize"; } };
struct extent { static std::string name() { return "extent"; } };
struct p0 { static std::string name() { return "p0"; } };
struct beta { static std::string name() { return "beta"; } };
struct betax { static std::string name() { return "betax"; } };
struct betay { static std::string name() { return "betay"; } };
struct betaz { static std::string name() { return "betaz"; } };
struct r0 { static std::string name() { return "r0"; } };
struct ce { static std::string name() { return "ce"; } };
struct alpha { static std::string name() { return "alpha"; } };
struct gamma { static std::string name() { return "gamma"; } };
struct pstiff { static std::string name() { return "pstiff"; } };
struct pde { static std::string name() { return "pde"; } };
struct pref { static std::string name() { return "pref"; } };
struct tolref { static std::string name() { return "tolref"; } };
struct ndofmax { static std::string name() { return "ndofmax"; } };
struct indicator{ static std::string name() { return "indicator"; } };
struct amr { static std::string name() { return "amr"; } };
struct ale { static std::string name() { return "ale"; } };
struct move { static std::string name() { return "move"; } };
struct tolderef { static std::string name() { return "tolderef"; } };
struct t0ref { static std::string name() { return "t0ref"; } };
struct dtref { static std::string name() { return "dtref"; } };
struct dtref_uniform { static std::string name() { return "dtref_uniform"; } };
struct maxlevels { static std::string name() { return "maxlevels"; } };
struct partitioner { static std::string name() { return "partitioner"; } };
struct partitioned { static std::string name() { return "partitioned"; } };
struct scheme { static std::string name() { return "scheme"; } };
struct nstep { static std::string name() { return "nstep"; } };
struct term { static std::string name() { return "term"; } };
struct t0 { static std::string name() { return "t0"; } };
struct dt { static std::string name() { return "dt"; } };
struct cfl { static std::string name() { return "cfl"; } };
struct dvcfl { static std::string name() { return "dvcfl"; } };
struct fct { static std::string name() { return "fct"; } };
struct fctclip { static std::string name() { return "fctclip"; } };
struct sys { static std::string name() { return "sys"; } };
struct sysfct { static std::string name() { return "sysfct"; } };
struct sysfctvar { static std::string name() { return "sysfctvar"; } };
struct fcteps { static std::string name() { return "fcteps"; } };
struct ctau { static std::string name() { return "ctau"; } };
struct refined { static std::string name() { return "refined output"; } };
struct matched {};
struct compatibility {};
struct bndint {};
struct part { static std::string name() { return "part"; } };
struct particles { static std::string name() { return "particles"; } };
struct centroid {};
struct ncomp { static std::string name() { return "ncomp"; } };
struct nmat { static std::string name() { return "nmat"; } };
struct tty { static std::string name() { return "tty"; } };
struct dump {};
struct plot {};
struct glob {};
struct control { static std::string name() { return "control"; } };
struct stat { static std::string name() { return "stat"; } };
struct field { static std::string name() { return "field"; } };
struct surface { static std::string name() { return "surface"; } };
struct integral { static std::string name() { return "integral"; } };
struct kappa { static std::string name() { return "kappa"; } };
struct rho2 { static std::string name() { return "rho2"; } };
struct rho { static std::string name() { return "rho"; } };
struct mean_gradient { static std::string name() { return "mean_gradient"; } };
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
struct refvar { static std::string name() { return "refvar"; } };
struct virtualization {static std::string name() { return "virtualization"; }};
struct normalization { static std::string name() { return "normalization"; } };
struct mass { static std::string name() { return "mass"; } };
struct title { static std::string name() { return "title"; } };
struct selected { static std::string name() { return "selected"; } };
struct discr { static std::string name() { return "discr"; } };
struct bc { static std::string name() { return "bc"; } };
struct farfield_pressure {
  static std::string name() { return "farfield_pressure"; } };
struct farfield_density {
  static std::string name() { return "farfield_density"; } };
struct farfield_velocity {
  static std::string name() { return "farfield_velocity"; } };
struct farfield { static std::string name() { return "farfield"; } };
struct pressure_pressure { static std::string name() { return "pressure_pressure"; } };
struct pressure_density { static std::string name() { return "pressure_density"; } };
struct component { static std::string name() { return "component"; } };
struct rescomp { static std::string name() { return "residual component"; } };
struct iter { static std::string name() { return "iter"; } };
struct time { static std::string name() { return "time"; } };
struct range { static std::string name() { return "range"; } };
struct cmd { static std::string name() { return "cmd"; } };
struct param { static std::string name() { return "param"; } };
struct init { static std::string name() { return "init"; } };
struct solve { static std::string name() { return "solve"; } };
struct chare { static std::string name() { return "chare"; } };
struct generator {};
struct help { static std::string name() { return "help"; } };
struct helpctr { static std::string name() { return "helpctr"; } };
struct helpkw {};
struct cmdinfo {};
struct ctrinfo {};
struct group { static std::string name() { return "group"; } };
struct esup {};
struct psup {};
struct gid {};
struct transport { static std::string name() { return "transport"; } };
struct compflow { static std::string name() { return "compflow"; } };
struct problem { static std::string name() { return "problem"; } };
struct physics { static std::string name() { return "physics"; } };
struct diffusivity { static std::string name() { return "diffusivity"; } };
struct u0 { static std::string name() { return "u0"; } };
struct dirichlet { static std::string name() { return "dirichlet"; } };
struct symmetry { static std::string name() { return "symmetry"; } };
struct point { static std::string name() { return "point"; } };
struct source { static std::string name() { return "source"; } };
struct radius { static std::string name() { return "radius"; } };
struct sideset { static std::string name() { return "sideset"; } };
struct material { static std::string name() { return "material"; } };
struct eos { static std::string name() { return "eos"; } };
struct id { static std::string name() { return "id"; } };
struct edge { static std::string name() { return "edge"; } };
struct cv { static std::string name() { return "cv"; } };
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
struct it {};
struct fn { static std::string name() { return "fn"; } };
struct fntype { static std::string name() { return "fntype"; } };
struct node {};
struct ic { static std::string name() { return "ic"; } };
struct velocity { static std::string name() { return "velocity"; } };
struct density { static std::string name() { return "density"; } };
struct pressure { static std::string name() { return "pressure"; } };
struct energy { static std::string name() { return "energy"; } };
struct energy_content { static std::string name() { return "energy_content"; } };
struct temperature { static std::string name() { return "temperature"; } };
struct box { static std::string name() { return "box"; } };
struct xmin { static std::string name() { return "xmin"; } };
struct xmax { static std::string name() { return "xmax"; } };
struct ymin { static std::string name() { return "ymin"; } };
struct ymax { static std::string name() { return "ymax"; } };
struct zmin { static std::string name() { return "zmin"; } };
struct zmax { static std::string name() { return "zmax"; } };

} // tag::

#endif // Tags_h
