/*! \file pirancslib.h pirate ncs library */
/* Copyright 2003-2005 Kevin Cowtan & University of York all rights reserved */

//L   This code is distributed under the terms and conditions of the
//L   CCP4 Program Suite Licence Agreement as a CCP4 Application.
//L   A copy of the CCP4 licence can be obtained by writing to the
//L   CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>


/*! Abstract base class for zero-th order function. */
class Target_fn_order_zero {
 public:
  Target_fn_order_zero() {}
  virtual ~Target_fn_order_zero() {}
  virtual int num_params() const = 0;
  virtual double operator() ( const std::vector<double>& args ) const = 0;
};


/*! Simplex optimiser. */
class Optimiser_simplex {
 public:
  enum TYPE { NORMAL, GRADIENT };
  Optimiser_simplex( double tolerance = 0.001, int max_cycles = 50, TYPE type = NORMAL );
  std::vector<double> operator() ( const Target_fn_order_zero& target_fn, const std::vector<std::vector<double> >& args ) const;
  void debug() { debug_mode = true; }
 private:
  double tolerance_, max_cycles_;
  TYPE type_;
  bool debug_mode;
  std::vector<double> params_;
};


/*! Create and store an NXmap from an orthogonal sub-region of an Xmap.
    This will usually be a sphere of density for use in a rotation fn. */
class NXmap_cube : public clipper::NXmap<float> {
 public:
  //! Null constructor
  NXmap_cube() {}
  NXmap_cube( const clipper::Xmap<float>& xmap, const clipper::Coord_orth& origin, const double& spacing, const double& extent );
  const clipper::Coord_orth& origin() const { return origin_; }
  const double& spacing() const { return spacing_; }
  const double& extent() const { return extent_; }

  inline const float& operator [] ( const clipper::Coord_orth& coord ) const
    { return get_data( clipper::Coord_grid( clipper::Util::intr(rx_*coord.x()+tx_), clipper::Util::intr(rx_*coord.y()+tx_), clipper::Util::intr(rx_*coord.z()+tx_) ) ); }
 private:
  clipper::Coord_orth origin_;
  double spacing_, extent_, rx_, tx_;
};


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
  Local_rtop inverse() const { return Local_rtop( rot().inverse(), tgt(), src() ); }
  clipper::RTop_orth rtop_orth() const { return clipper::RTop_orth( rot_.matrix(), tgt_ - rot_.matrix()*src_ ); }
  Local_rtop transform( const clipper::RTop_orth& r1, const clipper::RTop_orth& r2 ) const;
  std::pair<double,Local_rtop> symm_match( const Local_rtop& other, const clipper::Spacegroup& spgr, const clipper::Cell& cell, const double& tol_dst, const double& tol_ang ) const;
  std::pair<int,double> atom_match( const clipper::Spacegroup& spgr, const clipper::Cell& cell, const std::vector<clipper::Coord_orth>& coords, const double& tol ) const;
  std::pair<int,int>    atom_loop ( const clipper::Spacegroup& spgr, const clipper::Cell& cell, const std::vector<clipper::Coord_orth>& coords, const double& tol ) const;
  Local_rtop proper( const clipper::Spacegroup& spgr, const clipper::Cell& cell ) const;
  bool operator < ( const Local_rtop& other ) const { return rot().w() < other.rot().w(); }
 private:
  clipper::Rotation rot_; clipper::Coord_orth src_, tgt_;
};


/*! Target function for nxmap rotation optimisation. */
class Target_fn_nxmap_rotation : public Target_fn_order_zero {
 public:
  enum TYPE { CORREL, RMSD };
  Target_fn_nxmap_rotation() {}
  Target_fn_nxmap_rotation( const NXmap_cube& source, const NXmap_cube& target, const double& rot_step, const double& trn_step, const TYPE type=CORREL );
  ~Target_fn_nxmap_rotation() {}
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
  const NXmap_cube* source_;
  const NXmap_cube* target_;
  double rot_step_, trn_step_;
  TYPE type_;
  clipper::Rotation rot_;
};


/*! Return a scored list of rotations mapping one map to another. */
class Search_nxmap_rotation {
 public:
  //! Null constructor
  Search_nxmap_rotation() {}
  Search_nxmap_rotation( const double& rot_step, const double& trn_step, const Target_fn_nxmap_rotation::TYPE type=Target_fn_nxmap_rotation::CORREL );
  std::vector<std::pair<double,clipper::Rotation> > operator() ( const NXmap_cube& source, const NXmap_cube& target ) const;
 private:
  double rot_step_, trn_step_;
  Target_fn_nxmap_rotation::TYPE type_;
};


/*! Gaussian with width and weight. */
class Gaussian_orth {
 public:
  Gaussian_orth() {}
  Gaussian_orth( const clipper::Coord_orth coord_orth, const double& half_width, const double& weight ) : coord_orth_(coord_orth), half_width_(half_width), weight_(weight) {}
  const clipper::Coord_orth& coord_orth() const { return coord_orth_; }
  const double& half_width()              const { return half_width_; }
  const double& weight()                  const { return weight_; }
  double evaluate( const double& radsq )  const { return weight_*exp(-0.5*radsq/(half_width_*half_width_)); }
 private:
  clipper::Coord_orth coord_orth_;
  double half_width_, weight_;
};


/*! List of Gaussians with widths and weights. */
typedef std::vector<Gaussian_orth> Gaussian_orth_list;


/*! Target function for nxmap rotation optimisation. */
class Target_fn_coord_rotation : public Target_fn_order_zero {
 public:
  Target_fn_coord_rotation() {}
  Target_fn_coord_rotation( const Gaussian_orth_list& source, const Gaussian_orth_list& target, const double& rot_step );
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
  const Gaussian_orth_list* source_;
  const Gaussian_orth_list* target_;
  double rot_step_;
  clipper::Rotation rot_;
};


/*! Return a scored list of rotations mapping one map to another. */
class Search_coord_rotation {
 public:
  //! Null constructor
  Search_coord_rotation() {}
  Search_coord_rotation( const double& rot_step, const int& max_match );
  std::vector<std::pair<double,clipper::Rotation> > operator() ( const Gaussian_orth_list& source, const Gaussian_orth_list& target ) const;
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
class Search_NCS_from_atom_map {
 public:
  Search_NCS_from_atom_map() {}
  Search_NCS_from_atom_map( const double& tol_dst, const double& tol_ang ) : tol_dst_(tol_dst), tol_ang_(tol_ang) {}
  std::vector<Local_rtop> operator() ( const clipper::Atom_list& atoms, const clipper::Xmap<float>& xmap ) const;
  std::vector<Local_rtop> include_inverse( const std::vector<Local_rtop>& ncsops, const clipper::Spacegroup& spgr, const clipper::Cell& cell ) const;
  std::vector<Local_rtop> exclude_inverse( const std::vector<Local_rtop>& ncsops, const clipper::Spacegroup& spgr, const clipper::Cell& cell ) const;
  static Local_rtop masked_refine_from_map( const clipper::Xmap<float>& xmap, const Local_rtop& nxop, const double& map_radius, const double& local_radius );
 private:
  double tol_ang_, tol_dst_;
};


/*! 1-dimensional Guassian distribution. */
class Gaussian_probability_1d {
 public:
  //! null constructor
  Gaussian_probability_1d() :
    mean_(0.0), sigma_(0.0) {}
  //! constructor
  Gaussian_probability_1d( const double& mean, const double& sigma ) :
    mean_(mean), sigma_(sigma) {}
  const double& mean()    const { return mean_; }   //!< get mean
  const double& std_dev() const { return sigma_; }  //!< get standard deviation
  //! return value
  double operator() ( const double& x ) 
    { return (0.3989422804/sigma_)*exp(-0.5*clipper::Util::sqr(0.5*(x-mean_)/sigma_)); }
  //! multiply two Gaussians
  friend Gaussian_probability_1d operator *(const Gaussian_probability_1d& g1, const Gaussian_probability_1d& g2) {
    double m1 = g1.mean();
    double m2 = g2.mean();
    double w1 = g1.std_dev()*g1.std_dev();
    double w2 = g2.std_dev()*g2.std_dev();
    return Gaussian_probability_1d( (m1*w2+m2*w1)/(w1+w2), sqrt((w1*w2)/(w1+w2)) );
  }
  //! multiply two Gaussians (+ is synonym for *)
  friend Gaussian_probability_1d operator +(const Gaussian_probability_1d& g1, const Gaussian_probability_1d& g2) { return g1+g2; }
 private:
  double mean_, sigma_;
};


/*! Xmap of 1-dimensional Guassian distributions.
  Includes methods for constraining density using NCS operators. */
class Xmap_ncs : public clipper::Xmap<Gaussian_probability_1d> {
 public:
  void restrain_ncs( const clipper::Xmap<float>& xmap, const Local_rtop& nxop, const double& map_radius, const double& local_radius );
  const Gaussian_probability_1d& operator= (const Gaussian_probability_1d& value) { return static_cast<clipper::Xmap<Gaussian_probability_1d>&>(*this) = value; }
};
