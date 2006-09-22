/*! \file pirate-ncsaver.h pirate ncs library */
/* Copyright 2003-2006 Kevin Cowtan & University of York all rights reserved */

#include "pirate-ncsfind.h"


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
    { return (0.3989422804/sigma_)*exp(-0.5*clipper::Util::sqr((x-mean_)/sigma_)); }
  //! multiply two Gaussians
  friend Gaussian_probability_1d operator *(const Gaussian_probability_1d& g1, const Gaussian_probability_1d& g2) {
    double m1 = g1.mean();
    double m2 = g2.mean();
    double w1 = g1.std_dev()*g1.std_dev();
    double w2 = g2.std_dev()*g2.std_dev();
    return Gaussian_probability_1d( (m1*w2+m2*w1)/(w1+w2), sqrt((w1*w2)/(w1+w2)) );
  }
  //! multiply two Gaussians (+ is synonym for *)
  friend Gaussian_probability_1d operator +(const Gaussian_probability_1d& g1, const Gaussian_probability_1d& g2) { return g1*g2; }
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
