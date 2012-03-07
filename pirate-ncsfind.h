/*! \file pirate-ncsfind.h pirate ncs library */
/* Copyright 2003-2006 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

#include "simplex-lib.h"

/*! Class for a candidate operator. */
class Local_rtop
{
 public:
  Local_rtop() {}
  Local_rtop( const clipper::Rotation& rot, const clipper::Coord_orth& src, const clipper::Coord_orth& tgt ) : rot_(rot), src_(src), tgt_(tgt) {}
  const clipper::Rotation& rot() const { return rot_; }
  const clipper::Coord_orth& src() const { return src_; }
  const clipper::Coord_orth& tgt() const { return tgt_; }
  clipper::Rotation& rot() { return rot_; }
  clipper::Coord_orth& src() { return src_; }
  clipper::Coord_orth& tgt() { return tgt_; }
  bool is_null() const { return rot().is_null(); }
  Local_rtop inverse() const { return Local_rtop( rot().inverse(), tgt(), src() ); }
  clipper::RTop_orth rtop_orth() const { return clipper::RTop_orth( rot_.matrix(), tgt_ - rot_.matrix()*src_ ); }
  Local_rtop transform( const clipper::RTop_orth& r1, const clipper::RTop_orth& r2 ) const;
  std::pair<double,Local_rtop> symm_match( const Local_rtop& other, const clipper::Spacegroup& spgr, const clipper::Cell& cell, const double& tol_dst, const double& tol_ang ) const;
  Local_rtop proper( const clipper::Spacegroup& spgr, const clipper::Cell& cell ) const;
  std::vector<Local_rtop> proper_closed_group( const clipper::Spacegroup& spgr, const clipper::Cell& cell, const std::vector<int>& orders, double tol_dst = 5.0 ) const;
  static std::vector<Local_rtop> include_inverse( const std::vector<Local_rtop>& ncsops, const clipper::Spacegroup& spgr, const clipper::Cell& cell, double tol_dst = 3.0, double tol_ang = 0.1 );
  static std::vector<Local_rtop> exclude_inverse( const std::vector<Local_rtop>& ncsops, const clipper::Spacegroup& spgr, const clipper::Cell& cell, double tol_dst = 3.0, double tol_ang = 0.1 );
  static std::vector<Local_rtop> tidy( const std::vector<Local_rtop>& ncsops, const clipper::Spacegroup& spgr, const clipper::Cell& cell, double tol_dst = 3.0, double tol_ang = 0.1 );
  bool operator < ( const Local_rtop& other ) const { return rot().w() < other.rot().w(); }
 private:
  clipper::Rotation rot_; clipper::Coord_orth src_, tgt_;
};


/*! Gaussian with width and weight. */
class Atom_overlap {
 public:
  Atom_overlap() {}
  Atom_overlap( const clipper::Coord_orth coord_orth, const double& half_width, const double& weight ) : coord_orth_(coord_orth), half_width_(half_width), weight_(weight) {}
  const clipper::Coord_orth& coord_orth() const { return coord_orth_; }
  const double& half_width()              const { return half_width_; }
  const double& weight()                  const { return weight_; }
  double evaluate( const double& radsq )  const { double z = 0.5*radsq/(half_width_*half_width_); return weight_/(1+z+z*z); }
  //double evaluate( const double& radsq )  const { return weight_*exp(-0.5*radsq/(half_width_*half_width_)); }
 private:
  clipper::Coord_orth coord_orth_;
  double half_width_, weight_;
};


/*! List of Gaussians with widths and weights. */
typedef std::vector<Atom_overlap> Atom_overlap_list;


/*! Target function for nxmap rotation optimisation. */
class Target_fn_coord_rotation : public Target_fn_order_zero {
 public:
  Target_fn_coord_rotation() {}
  Target_fn_coord_rotation( const Atom_overlap_list& source, const Atom_overlap_list& target, const double& rot_step );
  ~Target_fn_coord_rotation() {}
  int num_params() const { return 3; }
  //! evaluate target function for given rotation
  double operator() ( const clipper::Rotation& rot ) const;
  //! \internal evaluate target function for EulerXYZr offset from rot_
  double operator() ( const std::vector<double>& args ) const;
  //! \internal convert params to rotation
  clipper::Rotation rotation( const std::vector<double>& args ) const;
  //! refine rotation
  clipper::Rotation refine( const clipper::Rotation& rot );
 private:
  const Atom_overlap_list* source_;
  const Atom_overlap_list* target_;
  double rot_step_;
  clipper::Rotation rot_;
};


/*! Return a scored list of rotations mapping one map to another. */
class Search_coord_rotation {
 public:
  //! Null constructor
  Search_coord_rotation() {}
  Search_coord_rotation( const double& rot_step, const int& max_match );
  std::vector<std::pair<double,clipper::Rotation> > operator() ( const Atom_overlap_list& source, const Atom_overlap_list& target ) const;
 private:
  double rot_step_;
  int max_match_;
};


/*! Target function for nxmap rotation optimisation. */
class Target_fn_xmap_rtop : public Target_fn_order_zero {
 public:
  enum TYPE { CORREL, RMSD };
  Target_fn_xmap_rtop() {}
  Target_fn_xmap_rtop( const clipper::Xmap<float>& source, const clipper::Xmap<float>& target, const double& rad, const double& rot_step, const double& trn_step, const TYPE type=CORREL );
  ~Target_fn_xmap_rtop() {}
  int num_params() const { return 6; }
  //! evaluate target function for given rotation
  double operator() ( const Local_rtop& rot ) const;
  //! \internal evaluate target function for EulerXYZr+uvw offset from rot_
  double operator() ( const std::vector<double>& args ) const;
  //! \internal convert params to rotation
  Local_rtop local_rtop( const std::vector<double>& args ) const;
  //! refine rotation
  Local_rtop refine( const Local_rtop& rot );
 private:
  const clipper::Xmap<float>* source_;
  const clipper::Xmap<float>* target_;
  double rad_, rot_step_, trn_step_;
  TYPE type_;
  Local_rtop rot_;
};


/*! Target function for nxmap rotation optimisation. */
class Target_fn_xmap_mask_rtop : public Target_fn_order_zero {
 public:
  enum TYPE { CORREL, RMSD };
  Target_fn_xmap_mask_rtop() {}
  Target_fn_xmap_mask_rtop( const clipper::Xmap<float>& source, const clipper::Xmap<float>& target, const clipper::NXmap<float>& mask, const double& rot_step, const double& trn_step );
  ~Target_fn_xmap_mask_rtop() {}
  int num_params() const { return 6; }
  //! evaluate target function for given rotation
  double operator() ( const Local_rtop& rot ) const;
  //! \internal evaluate target function for EulerXYZr+uvw offset from rot_
  double operator() ( const std::vector<double>& args ) const;
  //! \internal convert params to rotation
  Local_rtop local_rtop( const std::vector<double>& args ) const;
  //! refine rotation
  Local_rtop refine( const Local_rtop& rot );
 private:
  const clipper::Xmap<float>* source_;
  const clipper::Xmap<float>* target_;
  const clipper::NXmap<float>* mask_;
  double rad_, rot_step_, trn_step_;
  Local_rtop rot_;
  double mmin, mmax;
};


/*! Method object for heavy-atom+map -> NCS calculation */
class NCSfind {
 public:
  NCSfind() { debug_ = false; }
  NCSfind( const double& tol_dst, const double& tol_ang ) : tol_dst_(tol_dst), tol_ang_(tol_ang) { debug_ = false; }
  std::vector<Local_rtop> operator() ( const clipper::Atom_list& atoms, const clipper::Xmap<float>& xmap ) const;
  std::vector<Local_rtop> find_ncs_candidates( const clipper::Atom_list& atoms, const clipper::Spacegroup& spgr, const clipper::Cell& cell ) const;
  std::vector<Local_rtop> filter_ncs_candidates( const std::vector<Local_rtop>& rtops, const clipper::Atom_list& atoms, const clipper::Xmap<float>& xmap, double sigma ) const;
  void debug() { debug_ = true; }
  static std::vector<double> local_variance_from_map( clipper::NXmap<float>& nxmap0, clipper::NXmap<float>& nxmap1, clipper::NXmap<float>& correl, clipper::NXmap<float>& rmsdev, clipper::NXmap<float>& mskdev, const clipper::Xmap<float>& xmap, const Local_rtop& nxop, const double& map_radius, const double& local_radius, const double& pncs0 );
  static Local_rtop masked_refine_from_map( const clipper::Xmap<float>& xmap, const Local_rtop& nxop, const double& map_radius, const double& local_radius );
 private:
  static double minimum_distance_between_lines( const clipper::Coord_orth& a0, const clipper::Coord_orth& b0, const clipper::Coord_orth& a1, const clipper::Coord_orth& b1 ); 
  double tol_ang_, tol_dst_;
  bool debug_;
};
