/*! \file piralib.h pirate library */
/* Copyright 2003-2005 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

#include <map>

#include "pirancslib.h"


//! Reflection data type: Argand gradient
typedef clipper::data32::F_phi ArgGrad;

//! Reflection data type: Argand curvature
class ArgCurv : private clipper::Datatype_base
{
  typedef clipper::Util Util;
  typedef clipper::ftype ftype;
  typedef clipper::xtype xtype;
  typedef clipper::String String;
public:
  ArgCurv() : caa_(Util::nan()), cab_(Util::nan()), cbb_(Util::nan()) {}
  void set_null() { Util::set_null(caa_); Util::set_null(cab_);
                    Util::set_null(cbb_); }
  static String type() { return "ArgCurv"; }
  void friedel() { clipper::Message::message( clipper::Message_fatal( "ArgCurv: transform not implemented" ) ); }
  void shift_phase(const ftype& dphi) { friedel(); }
  bool missing() const { return (Util::is_nan(caa_) || Util::is_nan(cab_) || Util::is_nan(cbb_)); }
  static int data_size() { return 3; }
  static String data_names() { return "Caa Cab Cbb"; }
  void data_export( xtype array[] ) const
    { array[0] = curv_aa(); array[1] = curv_ab(); array[2] = curv_bb(); }
  void data_import( const xtype array[] )
    { curv_aa() = array[0]; curv_ab() = array[1]; curv_bb() = array[2]; }
  //! resolve along phase direction
  float resolve(const float phi) { const float va=cos(phi); const float vb=sin(phi); return va*va*caa_+2.0*va*vb*cab_+vb*vb*cbb_; }
  // accessors
  const float& curv_aa() const { return caa_; }  //<! read access
  const float& curv_ab() const { return cab_; }  //<! read access
  const float& curv_bb() const { return cbb_; }  //<! read access
  float& curv_aa() { return caa_; }  //<! write access
  float& curv_ab() { return cab_; }  //<! write access
  float& curv_bb() { return cbb_; }  //<! write access
private:
  float caa_,cab_,cbb_;
};


//! Base class for map likelihood functions
/*!
 */
class Xmap_target_base {
 public:
  //! calculate value of the target function
  virtual double func( const clipper::HKL_data<clipper::data32::F_phi>& fphi ) const = 0;
  //! calculate (diagonal) derivatives of the target w.r.t. structure factors
  virtual void curv( const clipper::HKL_data<clipper::data32::F_phi>& fphi, clipper::HKL_data<ArgGrad>& grad, clipper::HKL_data<ArgCurv>& curv ) const = 0;

  void debug( const clipper::HKL_data<clipper::data32::F_phi>& fphi ) const;
 protected:
  static double degrees_of_freedom( const clipper::HKL_data_base& data );
};


//! Class for compact representation of a probability histogram
/*! This class represents probability distributions for unitary
  variables, and in particular unitary electron densities, using the
  sun of a set of Gaussians of equal width and spacing. The Gaussian
  coefficients are represented by bytes holding the square-root of the
  height (i.e. quadratic compression to increase dynamic range). The
  width of the Gaussians is determined such the the inflections of
  neighbouring Gaussians intersect, i.e. spacing is 2 sigma. The
  number of Gaussians is determined by the template parameter. The
  overall scale is lost, so these are more useful for storing
  log-likelioohs.
*/
template<int N> class MiniHist {
 public:
  //! null constructor
  MiniHist() {}
  //! constructor: initialise from data
  MiniHist( const clipper::Histogram& hist );
  //! return value and derivatives of function
  void func_curv( const double& x, double& func, double& grad, double& curv ) const;
  //! return value and derivatives of log-likelihood of function
  void llk_curv( const double& x, double& func, double& grad, double& curv ) const;
  //! return mean of function
  double mean() const;
  //! read/write member data directly
  unsigned char& operator[] ( const int& i ) { return data[i]; }
 private:
  unsigned char data[N];
};


//! Class for representation of a probability histogram
/*! This class represents probability distributions for unitary
  variables, and in particular unitary electron densities, using the
  sun of a set of Gaussians of equal width and spacing. The
  width of the Gaussians is determined such the the inflections of
  neighbouring Gaussians intersect, i.e. spacing is 2 sigma. The
  number of Gaussians is determined by the template parameter. The
  overall scale is lost, so these are more useful for storing
  log-likelihoods.
*/
template<int N> class GaussianHistogram {
 public:
  //! null constructor
  GaussianHistogram() {}
  //! constructor: initialise from data
  GaussianHistogram( const clipper::Histogram& hist );
  //! constructor: initialise from range and Gaussian params
  GaussianHistogram( const clipper::Range<double>& range, const double& mean, const double& std_dev );
  //! read member data directly
  const double& operator[] ( const int& i ) const { return data[i]; }
  //! combine two histograms (with same range)
  friend GaussianHistogram<N> operator *(const GaussianHistogram<N>& g1, const GaussianHistogram<N>& g2) {
    GaussianHistogram<N> result = g1;
    for ( int i = 0; i < N; i++ ) result.data[i] *= g2.data[i];
    double f = 0.0;
    for ( int i = 0; i < N; i++ ) f = clipper::Util::max( f, result.data[i] );
    for ( int i = 0; i < N; i++ ) result.data[i] /= f;
    return result;
  }
 private:
  double data[N];
  clipper::Range<double> range_;
};


//! Class for compact representation of a probability histogram
/*! This class represents probability distributions for unitary
  variables, and in particular unitary electron densities, using the
  sun of a set of Gaussians of equal width and spacing. The Gaussian
  coefficients are represented by bytes holding the square-root of the
  height (i.e. quadratic compression to increase dynamic range). The
  width of the Gaussians is determined such the the inflections of
  neighbouring Gaussians intersect, i.e. spacing is 2 sigma. The
  number of Gaussians is determined by the template parameter. The
  overall scale is lost, so these are more useful for storing
  log-likelihoods.
*/
template<int N> class GaussianHistogramCompressed {
 public:
  //! null constructor
  GaussianHistogramCompressed() {}
  //! constructor: initialise from GaussianHistogram
  GaussianHistogramCompressed( const GaussianHistogram<N>& ghist );
  //! constructor: initialise from data
  GaussianHistogramCompressed( const clipper::Histogram& hist );
  //! return value and derivatives of function
  void func_curv( const double& x, double& func, double& grad, double& curv ) const;
  //! return value and derivatives of log-likelihood of function
  void llk_curv( const double& x, double& func, double& grad, double& curv ) const;
  //! return mean of function
  double mean() const;
  //! read member data directly
  const unsigned char& operator[] ( const int& i ) const { return data[i]; }
 private:
  unsigned char data[N];
};


typedef GaussianHistogram<16> TargetHist;
typedef GaussianHistogramCompressed<16> TargetHistCompr;


//! Class for calculating a probability histogram based Xmap target
/*!
 */
class Xmap_target_minihist : public Xmap_target_base {
 public:
  //! initialise the target function
  void init( const clipper::Xmap<TargetHistCompr>& target, const clipper::Range<double>& range );
  //! calculate value of the target function
  double func( const clipper::HKL_data<clipper::data32::F_phi>& fphi ) const;
  //! calculate (diagonal) derivatives of the target w.r.t. structure factors
  void curv( const clipper::HKL_data<clipper::data32::F_phi>& fphi, clipper::HKL_data<ArgGrad>& grad, clipper::HKL_data<ArgCurv>& curv ) const;

  //! return the mean density of these distributions
  double mean() const;
  //! normalise the range of the histograms to match given mean density
  void set_mean( const double& newmean );
 private:
  const clipper::Xmap<TargetHistCompr>* target_;
  clipper::Range<double> range_;
};


//! Class for calculating a gaussian Xmap target
/*!
 */
class Xmap_target_gaussian : public Xmap_target_base {
 public:
  //! initialise the target function
  void init( const clipper::Xmap<float>& target, const clipper::Xmap<float>& weight );
  //! calculate value of the target function
  double func( const clipper::HKL_data<clipper::data32::F_phi>& fphi ) const;
  //! calculate (diagonal) derivatives of the target w.r.t. structure factors
  void curv( const clipper::HKL_data<clipper::data32::F_phi>& fphi, clipper::HKL_data<ArgGrad>& grad, clipper::HKL_data<ArgCurv>& curv ) const;
 private:
  const clipper::Xmap<float>* target_;
  const clipper::Xmap<float>* weight_;
};


//! Class for refining HL coeffs from a Gaussian Xmap target
/*!
 */
class Refine_HL_coeff {
 public:
  Refine_HL_coeff() {}
  Refine_HL_coeff( const std::vector<std::pair<double,double> >& w_cyc, const double llkscale = 0.1 );
  void init( const std::vector<std::pair<double,double> >& w_cyc, const double llkscale = 0.1 );
  bool operator() ( clipper::HKL_data<clipper::data32::F_phi>& fphi,
		    clipper::HKL_data<clipper::data32::ABCD>& abcd_new,
		    const clipper::HKL_data<clipper::data32::ABCD>& abcd,
		    const clipper::HKL_data<clipper::data32::F_sigF>& fsig,
		    const clipper::HKL_data<clipper::data32::F_sigF>& fobs,
		    const Xmap_target_base& xtgt );
  double r_factor_work() { return rfac_w; }
  double r_factor_free() { return rfac_f; }
  double e_correl_work() { return ecor_w; }
  double e_correl_free() { return ecor_f; }
  double f_correl_work() { return fcor_w; }
  double f_correl_free() { return fcor_f; }
  double llk_gain_work() { return llkg_w; }
  double llk_gain_free() { return llkg_f; }
  static void debug( const clipper::HKL_data<clipper::data32::ABCD>& abcd,
		     const clipper::HKL_data<clipper::data32::F_sigF>& fsig,
		     const Xmap_target_base& xtgt );
 private:
  std::vector<std::pair<double,double> > w_cyc_;
  double llkscale_, rfac_w, rfac_f, ecor_w, ecor_f, fcor_w, fcor_f, llkg_w, llkg_f;
};


//! Smooth sphere or shell radial map filter
/*! This function implements a radial shell function with an optional
  hollow core and smooth drop-off. */
class MapFilterFn_shell : public clipper::MapFilterFn_base {
 public:
  MapFilterFn_shell() {}
  //! constructor: takes radius for step function cutoff
  MapFilterFn_shell( const clipper::ftype& inner, const clipper::ftype& outer, const bool hollow=false );
  //! evaluate radial shell function
  clipper::ftype operator() ( const clipper::ftype& radius ) const;
 private:
  clipper::ftype inner_, outer_;
};


//! Local map stats helper class
/*! This class calculates local map statistics in the form or
  ordinalised local mean and variance. */
class Map_local_moment_ordinal {
 public:
  Map_local_moment_ordinal( const clipper::Xmap<float>& xmap, const clipper::MapFilterFn_base& fn );
  double ord_moment_1( const clipper::Xmap<float>::Map_reference_index& ix ) const { return ord_mom1.ordinal( lmom1[ix] ); }
  double ord_moment_2( const clipper::Xmap<float>::Map_reference_index& ix ) const { return ord_mom2.ordinal( lmom2[ix] ); }
 private:
  clipper::Xmap<float> lmom1, lmom2;
  clipper::Generic_ordinal ord_mom1, ord_mom2;
};


//! Refine_HL_simulate phase improvement class
/*! This class is a method object which will refine HL coeffs on the
  work structure using map coefficients from a reference structure to
  construct density targets. */
class Refine_HL_simulate {
 public:
  Refine_HL_simulate() {}
  Refine_HL_simulate( bool un_bias, double rad_inner, double rad_outer, double ncs_radius = 6.0, double ncs_volume = 4.0, double weight_llk = 0.1, double weight_ramp = 2.0, double skew_moment1 = 0.0, double skew_moment2 = 0.0, int nbins_moment1 = 9, int nbins_moment2 = 9, int ncyc_int = 10, double oversampling = 1.5 );
  void init( bool un_bias, double rad_inner, double rad_outer, double ncs_radius = 6.0, double ncs_volume = 4.0, double weight_llk = 0.1, double weight_ramp = 2.0, double skew_moment1 = 0.0, double skew_moment2 = 0.0, int nbins_moment1 = 9, int nbins_moment2 = 9, int ncyc_int = 10, double oversampling = 1.5 );
  bool operator() ( clipper::HKL_data<clipper::data32::F_phi>& fphi,
		    clipper::HKL_data<clipper::data32::ABCD>& abcd_new,
		    const clipper::HKL_data<clipper::data32::F_sigF>& fsig,
		    const clipper::HKL_data<clipper::data32::F_sigF>& fobs,
		    const clipper::HKL_data<clipper::data32::ABCD>& abcd,
		    const clipper::HKL_data<clipper::data32::F_sigF>& ref_f,
		    const clipper::HKL_data<clipper::data32::ABCD>& ref_hlcal,
		    const clipper::HKL_data<clipper::data32::ABCD>& ref_hlsim );
  bool operator() ( clipper::HKL_data<clipper::data32::F_phi>& fphi,
		    clipper::HKL_data<clipper::data32::ABCD>& abcd_new,
		    const clipper::HKL_data<clipper::data32::F_sigF>& fsig,
		    const clipper::HKL_data<clipper::data32::F_sigF>& fobs,
		    const clipper::HKL_data<clipper::data32::ABCD>& abcd,
		    const clipper::HKL_data<clipper::data32::F_sigF>& ref_f,
		    const clipper::HKL_data<clipper::data32::ABCD>& ref_hlcal,
		    const clipper::HKL_data<clipper::data32::ABCD>& ref_hlsim,
		    const std::vector<Local_rtop>& ncsops );
  double r_factor_work() const { return rfac_w; }  //! return stats
  double r_factor_free() const { return rfac_f; }  //! return stats
  double e_correl_work() const { return ecor_w; }  //! return stats
  double e_correl_free() const { return ecor_f; }  //! return stats
  double f_correl_work() const { return fcor_w; }  //! return stats
  double f_correl_free() const { return fcor_f; }  //! return stats
  double llk_gain_work() const { return llkg_w; }  //! return stats
  double llk_gain_free() const { return llkg_f; }  //! return stats
  const std::vector<Local_rtop>& ncs_operators() const { return ncsops_; }
  const clipper::Array2d<clipper::Histogram>& hist_raw() const { return hist0; }
  const clipper::Array2d<clipper::Histogram>& hist_mod() const { return hist1; }
  const double& rmsd_calc() const { return rms_cal; }  //! verbose info
  const double& rmsd_simu() const { return rms_sim; }  //! verbose info
  const double& rmsd_work() const { return rms_wrk; }  //! verbose info
  void debug() { debug_mode = true; }
 private:
  // control parameters
  bool unbias;
  double rad1, rad2, ncsrad, ncsvol, wtllk, wtrmp, skew1, skew2, oversam;
  int ncyc, nbins1, nbins2;
  // results
  double rfac_w, rfac_f, ecor_w, ecor_f, fcor_w, fcor_f, llkg_w, llkg_f;
  std::vector<Local_rtop> ncsops_;
  // extra information for debugging
  bool debug_mode;
  clipper::Array2d<clipper::Histogram> hist0, hist1;
  double rms_cal, rms_sim, rms_wrk;
};


//! Map simulation class
/* This class simulates a set of HL coeffs for a reference structure
   (A & B only) matching the properties of the coefficients of a work
   structure. */
class MapSimulate {
 public:
  MapSimulate( int nresbins = 100, int binmin = 20 );
  bool operator() ( clipper::HKL_data<clipper::data32::F_sigF>& sim_f,
		    clipper::HKL_data<clipper::data32::ABCD>& sim_hl,
		    const clipper::HKL_data<clipper::data32::F_sigF>& ref_f,
		    const clipper::HKL_data<clipper::data32::ABCD>& ref_hl,
		    const clipper::HKL_data<clipper::data32::F_sigF>& wrk_f,
		    const clipper::HKL_data<clipper::data32::ABCD>& wrk_hl ) const;
 private:
  class EMagCompare {
  public:
    EMagCompare(const clipper::HKL_data<clipper::data32::E_sigE>& m ) { p = &m; }
    bool operator() ( const int& i1, const int& i2 ) const { return (*p)[i1].E() < (*p)[i2].E(); }
    const clipper::HKL_data<clipper::data32::E_sigE>* p;
  };
  int n_res_bins, bin_min;
};

