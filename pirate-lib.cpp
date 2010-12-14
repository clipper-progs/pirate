/*! \file pirate-lib.cpp pirate library */
/* Copyright 2003-2005 Kevin Cowtan & University of York all rights reserved */

#include "pirate-lib.h"

#include <clipper/clipper-ccp4.h>

extern "C" {
#include <stdlib.h>
}

// Histogram smoother

class HistogramSmoother {
 public:
  static clipper::Histogram smooth_region( const clipper::Histogram src, int i1, int i2 );
  static clipper::Histogram smooth_curvature( clipper::Histogram src ); 
};

clipper::Histogram HistogramSmoother::smooth_region( const clipper::Histogram src, int i1, int i2 )
{
  if ( i1 < 0 )            i1 = 0;
  if ( i2 > src.size()-1 ) i2 = src.size()-1;
  clipper::Histogram tgt( src, src.size() );
  for ( int i = 0; i < i1; i++ )
    tgt.accumulate( src.x(i), src.y(i) );
  tgt.accumulate( src.x(i1), 0.75*src.y(i1)+0.25*src.y(i1+1) );
  for ( int i = i1+1; i <= i2-1; i++ )
    tgt.accumulate( src.x(i), 0.25*src.y(i-1)+0.5*src.y(i)+0.25*src.y(i+1) );
  tgt.accumulate( src.x(i2), 0.25*src.y(i2-1)+0.75*src.y(i2) );
  for ( int i = i2+1; i < src.size(); i++ )
    tgt.accumulate( src.x(i), src.y(i) );
  return tgt;
}

clipper::Histogram HistogramSmoother::smooth_curvature( clipper::Histogram src )
{
  // initial smoothing
  clipper::Histogram tgt = smooth_region( src, 0, src.size()-1 );

  // specfic smoothing of negative curvature regions
  for ( int c = 0; c < 100; c++ ) {
    /*
    // locate maximum to preserve from smoothing
    int imax = 0;
    for ( int i = 0; i < tgt.size(); i++ )
      if ( tgt.y(i) > tgt.y(imax) ) imax = i;
    // calc curvatures
    std::vector<double> crv( tgt.size(), 0.0 );
    for ( int i = 1; i < crv.size()-1; i++ )
      crv[i] = tgt.y(i-1) - 2.0*tgt.y(i) + tgt.y(i+1);
    // tabulate curvature regions
    std::vector<std::pair<int,int> > regions;
    int i1, i2;
    for ( i1 = 1; i1 < crv.size(); i1++ )
      if ( crv[i1-1] >= 0.0 && crv[i1] < 0.0 ) {
	for ( i2 = i1; i2 < crv.size()-1; i2++ )
	  if ( crv[i2+1] >= 0.0 ) break;
	if ( i1 > imax || i2 < imax )
	  regions.push_back( std::pair<int,int>( i1, i2 ) );
      }
    // now smooth them
    if ( regions.size() > 0 )
      tgt = smooth_region( tgt, regions[0].first-1, regions[regions.size()-1].second+1 );
    for ( int r = 0; r < regions.size(); r++ )
      tgt = smooth_region( tgt, regions[r].first-1, regions[r].second+1 );
    */
    int nmax = 0;
    for ( int i = 1; i < tgt.size()-1; i++ )
      if ( tgt.y(i) > tgt.y(i-1) && tgt.y(i) > tgt.y(i+1) ) nmax++;
    if ( nmax > 1 ) tgt = smooth_region( tgt, 0, tgt.size()-1 );
  }

  // done
  return tgt;
}


// GaussianHistogram

template<int N> GaussianHistogram<N>::GaussianHistogram( const clipper::Histogram& hist )
{
  // fit params
  std::vector<double> params( N, 1.0 );
  std::vector<double> grad( N ), step( N ), w( N );
  clipper::Matrix<double> curv( N, N );
  double x0, x1, x2, f, val;
  const double d = 0.5;  // edge offset
  const double s = double(params.size()) + 2.0*d;  // histogram scale
  for ( int i = 0; i < 3; i++ ) {
    // calculate function, gradients and curvatures
    val = 0.0;  // zero everything
    for ( int i = 0; i < params.size(); i++ ) grad[i] = 0.0;
    for ( int i = 0; i < params.size(); i++ )
      for ( int j = 0; j < params.size(); j++ ) curv(i,j) = 0.0;
    for ( int i = 0; i < hist.size(); i++ ) {
      x0 = s * ( hist.x(i) - hist.min()) / hist.range();
      for ( int j = 0; j < params.size(); j++ ) {
	x1 = double(j) + 0.5 + d;
	w[j] = exp( -2.0 * clipper::Util::sqr( ((1.0))* (x0-x1) ) );
      }
      f = 0.0;
      for ( int j = 0; j < params.size(); j++ ) f += params[j] * w[j];
      val += clipper::Util::sqr( f - hist.y(i) );
      for ( int j = 0; j < params.size(); j++ ) {
	grad[j] += 2.0 * ( f - hist.y(i) ) * w[j];
	for ( int k = 0; k < params.size(); k++ )
	  curv(j,k) += 2.0 * w[j] * w[k];
      }
    }
    // solve for Newton Raphson step and apply
    step = curv.solve( grad );
    for ( int j = 0; j < step.size(); j++ )
      params[j] = clipper::Util::max( params[j]-step[j], 0.0 );
  }

  // store
  f = 0.0;
  for ( int i = 0; i < N; i++ ) f = clipper::Util::max( f, params[i] );
  for ( int i = 0; i < N; i++ ) data[i] = params[i]/f;

  // store range
  range_ = hist;
}


template<int N> GaussianHistogram<N>::GaussianHistogram( const clipper::Range<double>& range, const clipper::Histogram& hist )
{
  // fit params
  std::vector<double> params( N, 1.0 );
  std::vector<double> grad( N ), step( N ), w( N );
  clipper::Matrix<double> curv( N, N );
  double x0, x1, x2, f, val;
  const double d = 0.5;  // edge offset
  const double s = double(params.size()) + 2.0*d;  // histogram scale
  int i0 = int( hist.size() * (range.min()-hist.min())/hist.range() ) - 1;
  int i1 = int( hist.size() * (range.max()-hist.min())/hist.range() ) + 1;
  //std::cout << i0 << " " << i1 << "  i0i1\n";
  for ( int c = 0; c < 3; c++ ) {
    // calculate function, gradients and curvatures
    val = 0.0;  // zero everything
    for ( int i = 0; i < params.size(); i++ ) grad[i] = 0.0;
    for ( int i = 0; i < params.size(); i++ )
      for ( int j = 0; j < params.size(); j++ ) curv(i,j) = 0.0;
    for ( int i = i0; i < i1; i++ ) {
      // fit points over a region wider than the histogram
      double y = 0.0;
      if ( i >= 0 && i < hist.size() ) y = hist.y(i);
      // calc weights
      x0 = s * ( hist.x(i) - range.min() ) / range.range();
      for ( int j = 0; j < params.size(); j++ ) {
	x1 = double(j) + 0.5 + d;
	w[j] = exp( -2.0 * clipper::Util::sqr( ((1.0))* (x0-x1) ) );
      }
      f = 0.0;
      for ( int j = 0; j < params.size(); j++ ) f += params[j] * w[j];
      // calc func/grad/curv
      val += clipper::Util::sqr( f - y );
      for ( int j = 0; j < params.size(); j++ ) {
	grad[j] += 2.0 * ( f - y ) * w[j];
	for ( int k = 0; k < params.size(); k++ )
	  curv(j,k) += 2.0 * w[j] * w[k];
      }
    }
    // solve for Newton Raphson step and apply
    step = curv.solve( grad );
    for ( int j = 0; j < step.size(); j++ )
      params[j] = clipper::Util::max( params[j]-step[j], 0.0 );
  }

  // store
  f = 0.0;
  for ( int i = 0; i < N; i++ ) f = clipper::Util::max( f, params[i] );
  for ( int i = 0; i < N; i++ ) data[i] = params[i]/f;

  //std::cout << "GH: "; for ( int i = 0; i < N; i++ ) std::cout << clipper::String(int(99.9*data[i]),2) << " "; std::cout << "\n";

  // store range
  range_ = range;
}


template<int N> GaussianHistogram<N>::GaussianHistogram( const clipper::Range<double>& range, const double& mean, const double& std_dev )
{
  // calculate Gaussian
  for ( int i = 0; i < N; i++ ) {
    double x = range.min() + range.range() * double(i+1)/double(N+1);
    data[i] = exp( -0.5*clipper::Util::sqr((x-mean)/std_dev) );
  }

  // store
  double f = 0.0;
  for ( int i = 0; i < N; i++ ) f = clipper::Util::max( f, data[i] );
  for ( int i = 0; i < N; i++ ) data[i] /= f;

  if (clipper::Util::isnan(data[0])) { std::cout << "GH " << range.min() << " " << range.max() << " " << mean << " " << std_dev << " " << f << std::endl; for ( int i = 0; i < N; i++ ) std::cout << i << " " << data[i] << "\n"; }
  // store range
  range_ = range;
}

template class GaussianHistogram<16>;
template class GaussianHistogram<24>;
template class GaussianHistogram<32>;
template class GaussianHistogram<48>;


// GaussianHistogramCompressed

template<int N> GaussianHistogramCompressed<N>::GaussianHistogramCompressed( const GaussianHistogram<N>& ghist )
{
  // compress and store
  for ( int i = 0; i < N; i++ )
    data[i] = clipper::Util::intr( 255 * sqrt( ghist[i] ) );
}

template<int N> GaussianHistogramCompressed<N>::GaussianHistogramCompressed( const clipper::Histogram& hist )
{
  // make GuassianHistogram
  GaussianHistogram<N> ghist( hist );
  // compress and store
  for ( int i = 0; i < N; i++ )
    data[i] = clipper::Util::intr( 255 * sqrt( ghist[i] ) );
}

template<int N> void GaussianHistogramCompressed<N>::func_curv( const double& x, double& func, double& grad, double& curv ) const
{
  const double s(N+1);
  const double b(4.0*s*s);
  double f, dx;
  func = grad = curv = 0.0;
  for ( int i = 0; i < N; i++ ) {
    dx = ((1.0))* ( x - double(i+1)/s );
    if ( fabs( dx * s ) < 8.0 ) {  // optimisation
      f = double( data[i] );
      f = f * f * exp( -0.5 * b * dx * dx );
      func += f;
      grad += -b * dx * f;
      curv += (b * f) * ( b * dx * dx - 1.0 );
    }
  }
}

template<int N> void GaussianHistogramCompressed<N>::llk_curv( const double& x, double& func, double& grad, double& curv ) const
{
  const double sml = 1.0e-50;
  func_curv( x, func, grad, curv );
  func = clipper::Util::max( func, sml );
  grad = grad / func;
  curv = curv / func - grad * grad;
  func = log( func );
}

template<int N> double GaussianHistogramCompressed<N>::mean() const
{
  const double s(N+1);
  double s0, s1, f, dx;
  s0 = s1 = 0.0;
  for ( int i = 0; i < N; i++ ) {
    dx = double(i+1)/s;
    f = double( data[i] );
    f = f * f;
    s0 += f;
    s1 += dx * f;
  }
  return s1 / s0;
}

template class GaussianHistogramCompressed<16>;
template class GaussianHistogramCompressed<24>;
template class GaussianHistogramCompressed<32>;
template class GaussianHistogramCompressed<48>;



// generic map likelihood target

double Xmap_target_base::degrees_of_freedom( const clipper::HKL_data_base& data )
{
  clipper::HKL_info::HKL_reference_index ih;
  double d = 0.0;
  for ( ih = data.first(); !ih.last(); ih.next() ) {
    if ( ih.hkl_class().centric() ) d += 1.0;
    else                            d += 2.0;
  }
  return d;
}

void Xmap_target_base::debug( const clipper::HKL_data<clipper::data32::F_phi>& fphi ) const
{
  clipper::HKL_data<clipper::data32::F_phi> fphid( fphi.base_hkl_info() );
  clipper::HKL_data<ArgGrad> dfphi( fphi.base_hkl_info() );
  clipper::HKL_data<ArgGrad> dfphia( fphi.base_hkl_info() );
  clipper::HKL_data<ArgGrad> dfphib( fphi.base_hkl_info() );
  clipper::HKL_data<ArgCurv> ddfphi( fphi.base_hkl_info() );
  clipper::HKL_data<ArgCurv> dummy( fphi.base_hkl_info() );
  double t, ta, tb;
  t = func( fphi );
  curv( fphi, dfphi, ddfphi );
  clipper::HKL_info::HKL_reference_index ih;
  for ( ih = fphi.first(); !ih.last(); ih.next() )
    if ( (abs(ih.hkl().h())<=4&&abs(ih.hkl().k())<=2&&abs(ih.hkl().l())<=1 ) ||
	 ih.index() % 1000 == 0 ) {
      double d = 0.1;
      double phic = 0.0;
      fphid = fphi;
      fphid[ih] = std::complex<float>( fphid[ih] ) +
                  std::complex<float>( std::polar( d, phic ) );
      ta = func( fphid );
      curv( fphid, dfphia, dummy );
      fphid = fphi;
      fphid[ih] = std::complex<float>( fphid[ih] ) +
                  std::complex<float>( std::polar( d, phic + 1.570796 ) );
      tb = func( fphid );
      curv( fphid, dfphib, dummy );
      std::complex<float> da = ( std::complex<float>( dfphia[ih] ) -
				 std::complex<float>( dfphi[ih] ) );
      std::complex<float> db = ( std::complex<float>( dfphib[ih] ) -
				 std::complex<float>( dfphi[ih] ) );
      float caa, cbb, cab;
      caa = ddfphi[ih].curv_aa();
      cab = ddfphi[ih].curv_ab();
      cbb = ddfphi[ih].curv_bb();

      float tol = 0.00002;
      std::cout << ih.hkl().h() << " " << ih.hkl().k() << " " << ih.hkl().l() << ", " << ih.hkl_class().epsilon();
      if ( !ih.hkl_class().centric() ) {
	std::cout << "\t(" << (ta-t)/d << ",\t" << (tb-t)/d << ")\t(" << std::complex<float>( dfphi[ih] ).real() << ",\t" << std::complex<float>( dfphi[ih] ).imag() << ")\n";
	std::cout << "[ " << da.real()/d << ",\t" << da.imag()/d << ",\t" << db.real()/d << ",\t" << db.imag()/d << " ]\n";
	std::cout << "[ " << caa << ",\t" << cab << ",\t" << cab << ",\t" << cbb << " ]\n";
	if ( fabs( caa - da.real()/d ) > tol || fabs( cbb - db.imag()/d ) > tol )
	  std::cout << "ERR1\n";
	if ( fabs( cab ) > tol )
	  if ( fabs( cab - db.real()/d ) > tol || fabs( cab - da.imag()/d ) > tol ) std::cout << "ERR2\n";
      } else {
	float va = cos( ih.hkl_class().allowed() );
	float vb = sin( ih.hkl_class().allowed() );
	float g1 = ((ta-t)*va+(tb-t)*vb)/d;
	float g2 = dfphi[ih].resolve( ih.hkl_class().allowed() );
	std::cout << "\t(" << g1 << ")\t(" << g2 << ")\n";
	float c1 =
	  (da.real()*va*va+da.imag()*va*vb+db.real()*va*vb+db.imag()*vb*vb)/d;
	float c2 = ddfphi[ih].resolve( ih.hkl_class().allowed() );
	std::cout << "[ " << c1 << " ]\n";
	std::cout << "[ " << c2 << " ]\n";
	if ( fabs( c1 - c2 ) > tol  ) std::cout << "ERRC\n";
      }
    }
}



// Minihist target

void Xmap_target_minihist::init( const clipper::Xmap<TargetHistCompr>& target, const clipper::Range<double>& range )
{
  target_ = &target;
  range_ = range;
}

double Xmap_target_minihist::mean() const
{
  const clipper::Xmap<TargetHistCompr>& target = *target_;
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  double fn = 0.0;
  for ( MRI ix = target.first(); !ix.last(); ix.next() )
    fn += target[ix].mean()
      * target.spacegroup().num_symops() / target.multiplicity(ix.coord());
  fn /= target.grid_sampling().size();
  return fn * range_.range() + range_.min();
}

void Xmap_target_minihist::set_mean( const double& newmean )
{
  double shift = newmean - mean();
  range_ = clipper::Range<double>( range_.min() + shift, range_.max() + shift );
}

double Xmap_target_minihist::func( const clipper::HKL_data<clipper::data32::F_phi>& fphi ) const
{
  const clipper::Xmap<TargetHistCompr>& target = *target_;

  double scale = degrees_of_freedom(fphi);

  // prepare gradient and curvature maps
  clipper::Xmap<float> xmap( target.spacegroup(), target.cell(),
			     target.grid_sampling() );
  // calculate target map
  xmap.fft_from( fphi );

  // calculate target function and derivs
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  double f, df, ddf;
  double s = range_.range();
  double fn = 0.0;
  for ( MRI ix = xmap.first(); !ix.last(); ix.next() ) {
    target[ix].llk_curv( (xmap[ix]-range_.min())/s, f, df, ddf );
    fn += f * xmap.spacegroup().num_symops() / xmap.multiplicity(ix.coord());
  }

  // scale function
  fn *= scale / xmap.grid_sampling().size();
  return fn;
}

void Xmap_target_minihist::curv( const clipper::HKL_data<clipper::data32::F_phi>& fphi, clipper::HKL_data<ArgGrad>& grad, clipper::HKL_data<ArgCurv>& curv ) const
{
  const clipper::Xmap<TargetHistCompr>& target = *target_;

  typedef clipper::Xmap<float>::Map_reference_index MRI;
  typedef clipper::HKL_info::HKL_reference_index HRI;

  double scale = degrees_of_freedom(fphi);

  // prepare gradient and curvature maps
  clipper::Xmap<float> xmap( target.spacegroup(), target.cell(),
			     target.grid_sampling() );
  clipper::FFTmap    fftmap( target.spacegroup(), target.cell(),
			     target.grid_sampling() );

  // calculate target map
  xmap.fft_from( fphi );

  // calculate target function and derivs
  double f, df, ddf;
  double s = range_.range();
  for ( MRI ix = xmap.first(); !ix.last(); ix.next() ) {
    target[ix].llk_curv( (xmap[ix]-range_.min())/s, f, df, ddf );
    xmap[ix] = df/s;
    fftmap.set_real_data( ix.coord(), 0.5*ddf/(s*s) );
  }

  // calc and scale gradient
  xmap.fft_to( grad );
  scale *= 2.0 / clipper::Util::sqr( xmap.cell().volume() );
  for ( HRI ih = fphi.first(); !ih.last(); ih.next() )
    grad[ih].f() *= scale * xmap.spacegroup().num_symops() /
      ih.hkl_class().epsilonc();

  // calc reflection curvatures
  fftmap.fft_x_to_h();
  clipper::HKL hkl, hi, hj;
  clipper::data32::F_phi fp;
  double caa, cab, cbb, di, dj, c;
  scale *= 2.0 / target.cell().volume();
  for ( HRI ih = curv.first(); !ih.last(); ih.next() ) {
    hkl = ih.hkl();
    caa = cab = cbb = 0.0;
    for ( int i = 0; i < target.spacegroup().num_symops(); i++ )
      for ( int j = 0; j < target.spacegroup().num_symops(); j++ ) {
	hi = hkl.transform( target.spacegroup().symop(i) );
	hj = hkl.transform( target.spacegroup().symop(j) );
	di = hkl.sym_phase_shift( target.spacegroup().symop(i) );
	dj = hkl.sym_phase_shift( target.spacegroup().symop(j) );
	fp = fftmap.get_recip_data( hi + hj );
	c = cos(fp.phi()-di-dj);
	caa +=  fp.f() * c;
	cbb += -fp.f() * c;
	cab +=  fp.f() * sin(fp.phi()-di-dj);
	fp = fftmap.get_recip_data( -hi + hj );
	c = cos(fp.phi()+di-dj);
	caa +=  fp.f() * c;
	cbb +=  fp.f() * c;
      }
    curv[ih].curv_aa()=caa*scale/clipper::Util::sqr(ih.hkl_class().epsilonc());
    curv[ih].curv_bb()=cbb*scale/clipper::Util::sqr(ih.hkl_class().epsilonc());
    curv[ih].curv_ab()=cab*scale/clipper::Util::sqr(ih.hkl_class().epsilonc());
  }
}



// Gaussuan target

void Xmap_target_gaussian::init( const clipper::Xmap<float>& target, const clipper::Xmap<float>& weight )
{
  target_ = &target;
  weight_ = &weight;
}

double Xmap_target_gaussian::func( const clipper::HKL_data<clipper::data32::F_phi>& fphi ) const
{
  const clipper::Xmap<float>& target = *target_;
  const clipper::Xmap<float>& weight = *weight_;

  // calculate target map
  clipper::Xmap<float> xmap( target.spacegroup(), target.cell(),
			     target.grid_sampling() );
  xmap.fft_from( fphi );

  double scale = -degrees_of_freedom(fphi);

  // calculate target function
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  double fn = 0.0;
  for ( MRI ix = xmap.first(); !ix.last(); ix.next() )
    fn += weight[ix] * clipper::Util::sqr( xmap[ix] - target[ix] )
      * xmap.spacegroup().num_symops() / xmap.multiplicity( ix.coord() );
  fn *= scale / xmap.grid_sampling().size();
  return fn;
}

void Xmap_target_gaussian::curv( const clipper::HKL_data<clipper::data32::F_phi>& fphi, clipper::HKL_data<ArgGrad>& grad, clipper::HKL_data<ArgCurv>& curv ) const
{
  const clipper::Xmap<float>& target = *target_;
  const clipper::Xmap<float>& weight = *weight_;

  typedef clipper::Xmap<float>::Map_reference_index MRI;
  typedef clipper::HKL_info::HKL_reference_index HRI;

  // calculate target map
  clipper::Xmap<float> xmap( target.spacegroup(), target.cell(),
			     target.grid_sampling() );
  xmap.fft_from( fphi );

  double scale = -degrees_of_freedom(fphi);

  // calculate gradient map
  for ( MRI ix = xmap.first(); !ix.last(); ix.next() )
    xmap[ix] = 2.0 * weight[ix]*( xmap[ix] - target[ix] );
  xmap.fft_to( grad );
  scale *= 2.0 / clipper::Util::sqr( xmap.cell().volume() );
  for ( HRI ih = fphi.first(); !ih.last(); ih.next() )
    grad[ih].f() *= scale * xmap.spacegroup().num_symops() /
      ih.hkl_class().epsilonc();

  // calc fft of weight map
  clipper::FFTmap fftmap( weight.spacegroup(), weight.cell(),
			  weight.grid_sampling() );
  for ( MRI ix = weight.first();	!ix.last(); ix.next() )
    fftmap.set_real_data( ix.coord(), weight[ix] );

  // calc reflection curvatures
  fftmap.fft_x_to_h();
  clipper::HKL hkl, hi, hj;
  clipper::data32::F_phi fp;
  double caa, cab, cbb, di, dj, c;
  scale *= 2.0 / weight.cell().volume();
  for ( HRI ih = curv.first(); !ih.last(); ih.next() ) {
    hkl = ih.hkl();
    caa = cab = cbb = 0.0;
    for ( int i = 0; i < weight.spacegroup().num_symops(); i++ )
      for ( int j = 0; j < weight.spacegroup().num_symops(); j++ ) {
	hi = hkl.transform( weight.spacegroup().symop(i) );
	hj = hkl.transform( weight.spacegroup().symop(j) );
	di = hkl.sym_phase_shift( weight.spacegroup().symop(i) );
	dj = hkl.sym_phase_shift( weight.spacegroup().symop(j) );
	fp = fftmap.get_recip_data( hi + hj );
	c = cos(fp.phi()-di-dj);
	caa +=  fp.f() * c;
	cbb += -fp.f() * c;
	cab +=  fp.f() * sin(fp.phi()-di-dj);
	fp = fftmap.get_recip_data( -hi + hj );
	c = cos(fp.phi()+di-dj);
	caa +=  fp.f() * c;
	cbb +=  fp.f() * c;
      }
    curv[ih].curv_aa()=caa*scale/clipper::Util::sqr(ih.hkl_class().epsilonc());
    curv[ih].curv_bb()=cbb*scale/clipper::Util::sqr(ih.hkl_class().epsilonc());
    curv[ih].curv_ab()=cab*scale/clipper::Util::sqr(ih.hkl_class().epsilonc());
  }
}



// Refinement function

Refine_HL_coeff::Refine_HL_coeff( const std::vector<std::pair<double,double> >& w_cyc, const double llkscale ) { init( w_cyc, llkscale ); }

void Refine_HL_coeff::init( const std::vector<std::pair<double,double> >& w_cyc, const double llkscale )
{
  w_cyc_ = w_cyc;
  llkscale_ = llkscale;
}


// helper subroutine
// calculate maximum of llk distribution (with Wilson contrib)
float llk_max( const clipper::data32::F_phi& fphi, const ArgGrad& grad, const ArgCurv& curv, const float& f2, const float& llkscale )
{
  // get f's
  float a = fphi.a();
  float b = fphi.b();
  // get grad and curvature
  float da  = grad.a()*llkscale;
  float db  = grad.b()*llkscale;
  float daa = curv.curv_aa()*llkscale;
  float dab = curv.curv_ab()*llkscale;
  float dbb = curv.curv_bb()*llkscale;
  // offset to get gradients at origin by subtracting current location
  da -= daa * a + dab * b;
  db -= dab * a + dbb * b;
  // now apply wilson probability distribution
  float dww = -1.0/f2;
  daa += dww;
  dbb += dww;
  // calculate most probable F
  return sqrt( ( clipper::Util::sqr(dbb*da-dab*db) +
		 clipper::Util::sqr(daa*db-dab*da) )
	       / clipper::Util::sqr(daa*dbb-dab*dab) );
}

// calculate llk distribution for a given F (with Wilson constrib)
clipper::LogPhaseProb<72> llk_dist( const clipper::data32::F_phi& fphi, const ArgGrad& grad, const ArgCurv& curv, const float& f2, const float& llkscale, const clipper::HKL_class& cls, const float& f )
{
  // get f's
  float a = fphi.a();
  float b = fphi.b();
  // get grad and curvature
  float da  = grad.a()*llkscale;
  float db  = grad.b()*llkscale;
  float daa = curv.curv_aa()*llkscale;
  float dab = curv.curv_ab()*llkscale;
  float dbb = curv.curv_bb()*llkscale;
  // offset to get gradients at origin by subtracting current location
  da -= daa * a + dab * b;
  db -= dab * a + dbb * b;
  // now apply wilson probability distribution
  float dww = -1.0/f2;
  daa += dww;
  dbb += dww;
  // get normalisation factor to set integral to 1
  float s;
  if ( !cls.centric() ) {
    s = log( clipper::Util::twopi()/sqrt(daa*dbb-dab*dab) );
  } else {
    a = cos( cls.allowed() );
    b = sin( cls.allowed() );
    s = log( clipper::Util::twopi()/fabs(a*a*daa+b*b*dbb+2.0*a*b*dab) );
  }
  // calculate offset to scale factor for non-origin centring of Gaussian
  float a0 = -(dbb*da-dab*db) / (daa*dbb-dab*dab);
  float b0 = -(daa*db-dab*da) / (daa*dbb-dab*dab);
  float c0 = a0*da + b0*db + 0.5*a0*a0*daa + a0*b0*dab + 0.5*b0*b0*dbb;
  // calculate phase probability distribution
  clipper::LogPhaseProb<72> q( cls );
  for ( int p = 0; p < q.size(); p ++ ) {
    a = f * cos( q.phase(p) );
    b = f * sin( q.phase(p) );
    q[p] = (a*da + b*db + 0.5*a*a*daa + a*b*dab + 0.5*b*b*dbb - c0) - s;
    //float a1 = a - a0; float b1 = b - b0;
    //float x = (a1*da + b1*db + 0.5*a1*a1*daa + a1*b1*dab + 0.5*b1*b1*dbb) - s;
    //if ( fabs(q[p]-x) > 0.1 ) std::cout << "ERR " << q[p] << " " << x << "\n";
  }
  return q;
}



bool Refine_HL_coeff::operator()
		  ( clipper::HKL_data<clipper::data32::F_phi>& fphi,
		    clipper::HKL_data<clipper::data32::ABCD>& abcd_new,
		    const clipper::HKL_data<clipper::data32::ABCD>& abcd,
		    const clipper::HKL_data<clipper::data32::F_sigF>& fsig,
		    const clipper::HKL_data<clipper::data32::F_sigF>& fobs,
		    const Xmap_target_base& xtgt )
{
  // set up derivative lists
  clipper::HKL_data<ArgGrad> grad( fobs.base_hkl_info() );
  clipper::HKL_data<ArgCurv> curv( fobs.base_hkl_info() );

  // first get phase and weight for first cycle
  clipper::data32::Phi_fom phifom;
  clipper::HKL_info::HKL_reference_index ih;
  for ( ih = fphi.first(); !ih.last(); ih.next() ) {
    clipper::LogPhaseProb<72> q( ih.hkl_class() );
    q.set_abcd( abcd[ih] );
    q.get_phi_fom( phifom );
    if ( !fsig[ih].missing() ) {
      fphi[ih].f()   = fsig[ih].f() * phifom.fom();
      fphi[ih].phi() = phifom.phi();
    } else {
      fphi[ih].f()   = fphi[ih].phi() = 0.0;
    }
  }

  // and get mean F^2
  const int n_scl_param = 10;
  std::vector<double> params( n_scl_param, 1.0 );
  clipper::BasisFn_spline basisfn( fobs, n_scl_param );
  clipper::TargetFn_meanFnth<clipper::data32::F_sigF> targetfn( fobs, 2.0 );
  clipper::ResolutionFn meanf2( fobs.base_hkl_info(), basisfn, targetfn, params );

  // main loop
  float f, f2, fm;
  for ( int cyc = 0; cyc < w_cyc_.size(); cyc++ ) {
    // calc function value and derivs
    xtgt.curv( fphi, grad, curv );

    // loop over reflections and make new probability distributions
    for ( ih = fphi.first(); !ih.last(); ih.next() )
      if ( ih.hkl() != clipper::HKL::zero() ) {
	clipper::LogPhaseProb<72> qexp(ih.hkl_class()), qmap(ih.hkl_class());
	qexp.set_abcd( abcd[ih] );
	f2 = ih.hkl_class().epsilonc() * meanf2.f(ih);
	if ( !fsig[ih].missing() ) {
	  // non-missing reflections - use measured f
	  f = fsig[ih].f();
	  qmap = llk_dist( fphi[ih], grad[ih], curv[ih], f2, llkscale_,
			   ih.hkl_class(), f );
	  for ( int p = 0; p < qmap.size(); p ++ )
	    qmap[p] = w_cyc_[cyc].first*qexp[p] + w_cyc_[cyc].second*qmap[p];
	  qmap.get_phi_fom( phifom );
	  fphi[ih].f() = f * phifom.fom();
	  fphi[ih].phi() = phifom.phi();
	} else {
	  // non-missing reflections - use best estmate of f
	  f = llk_max( fphi[ih], grad[ih], curv[ih], f2, llkscale_ );
	  qmap = llk_dist( fphi[ih], grad[ih], curv[ih], f2, llkscale_,
			   ih.hkl_class(), f );
	  for ( int p = 0; p < qmap.size(); p ++ )
	    qmap[p] = w_cyc_[cyc].first*qexp[p] + w_cyc_[cyc].second*qmap[p];
	  qmap.get_phi_fom( phifom );
	  fphi[ih].f() = f;
	  fphi[ih].phi() = phifom.phi();
	}
	qmap.get_abcd( abcd_new[ih] );
      }
  }

  // now do a statistics run
  clipper::HKL_data<clipper::data32::E_sigE>
    eo( fobs.base_hkl_info() ), ec( fobs.base_hkl_info() );
  ArgGrad gradnull; ArgCurv curvnull;
  gradnull.f() = gradnull.phi() = 0.0;
  curvnull.curv_aa() = curvnull.curv_bb() = curvnull.curv_ab() = 0.0;
  double nw, nf, lw, lf, rw, rf, dw, df, s1, s2;
  double fo1, fc1, eo1, ec1;
  double eow, eof, ecw, ecf, eoow, eoof, eccw, eccf, eocw, eocf; 
  double fow, fof, fcw, fcf, foow, foof, fccw, fccf, focw, focf; 
  nw = nf = lw = lf = rw = rf = dw = df = 0.0;
  eow = eof = ecw = ecf = eoow = eoof = eccw = eccf = eocw = eocf = 0.0; 
  fow = fof = fcw = fcf = foow = foof = fccw = fccf = focw = focf = 0.0; 
  xtgt.curv( fphi, grad, curv );
  // first get most probable F's to predict E's
  for ( ih = fphi.first(); !ih.last(); ih.next() )
    if ( ih.hkl() != clipper::HKL::zero() && !fobs[ih].missing() ) {
      fm = llk_max( fphi[ih], grad[ih], curv[ih], f2, llkscale_ );
      eo[ih].E() = fobs[ih].f() / sqrt(ih.hkl_class().epsilon());
      ec[ih].E() = fm           / sqrt(ih.hkl_class().epsilon());
      eo[ih].sigE() = ec[ih].sigE() = 1.0;
    }
  // now calc E's
  clipper::TargetFn_scaleEsq<clipper::data32::E_sigE> target_fo( eo );
  clipper::ResolutionFn scalefosq( fobs.base_hkl_info(), basisfn, target_fo, params );
  clipper::TargetFn_scaleEsq<clipper::data32::E_sigE> target_fc( ec );
  clipper::ResolutionFn scalefcsq( fobs.base_hkl_info(), basisfn, target_fc, params );
  // now do stats
  for ( ih = fphi.first(); !ih.last(); ih.next() )
    if ( ih.hkl() != clipper::HKL::zero() && !fobs[ih].missing() ) {
      // calculate most probable F
      f2 = ih.hkl_class().epsilonc() * meanf2.f(ih);
      fm = llk_max( fphi[ih], grad[ih], curv[ih], f2, llkscale_ );
      fo1 = eo[ih].E();
      fc1 = ec[ih].E();
      eo1 = fo1 * sqrt( clipper::Util::max( scalefosq.f(ih), 0.0 ) );
      ec1 = fc1 * sqrt( clipper::Util::max( scalefcsq.f(ih), 0.0 ) );
      // calculate phase probability distributions
      clipper::LogPhaseProb<72> qexp(ih.hkl_class()),
	  qnul(ih.hkl_class()), qmap(ih.hkl_class());
      qexp.set_abcd( abcd[ih] );
      qnul = llk_dist( gradnull, gradnull, curvnull, f2, llkscale_,
		       ih.hkl_class(), fobs[ih].f() );
      qmap = llk_dist( fphi[ih], grad[ih], curv[ih], f2, llkscale_,
		       ih.hkl_class(), fobs[ih].f() );
      // calc llg
      s1 = s2 = 0.0;  // sum exp, map, posterior probs
      for ( int p = 0; p < qmap.size(); p ++ ) {
	s1 += exp( qexp[p] + qnul[p] );
	s2 += exp( qexp[p] + qmap[p] );
      }
      double llg = log( s2 / s1 );
      // calc R, R-free, LLK, LLK-free
      /*
      f = fobs[ih].f()+0.05*sqrt(f2);
      l1 = 0.0;
      for ( int p = 0; p < qmap.size(); p ++ )
	l1 += clipper::Util::twopi()/double(qmap.size()) * f * exp( qmap[p] );
      l2 = f * exp( - fobs[ih].f()*fobs[ih].f()/f2 ) / f2;
      */
      if ( !fsig[ih].missing() ) {
	nw += 1.0;
	lw += llg;
	rw += fabs( fm - fobs[ih].f() );
	dw += fobs[ih].f();
	eow += eo1;
	ecw += ec1;
	eoow += eo1 * eo1;
	eccw += ec1 * ec1;
	eocw += eo1 * ec1;
	fow += fo1;
	fcw += fc1;
	foow += fo1 * fo1;
	fccw += fc1 * fc1;
	focw += fo1 * fc1;
      } else {
	nf += 1.0;
	lf += llg;
	rf += fabs( fm - fobs[ih].f() );
	df += fobs[ih].f();
	eof += eo1;
	ecf += ec1;
	eoof += eo1 * eo1;
	eccf += ec1 * ec1;
	eocf += eo1 * ec1;
	fof += fo1;
	fcf += fc1;
	foof += fo1 * fo1;
	fccf += fc1 * fc1;
	focf += fo1 * fc1;
      }
    }

  // store free results
  rfac_w = rw/dw;
  rfac_f = (nf>0.0)?(rf/df):0.0;
  ecor_w = (eocw*nw-eow*ecw) / sqrt((eoow*nw-eow*eow)*(eccw*nw-ecw*ecw));
  ecor_f = (nf>0.0)?((eocf*nf-eof*ecf) / sqrt((eoof*nf-eof*eof)*(eccf*nf-ecf*ecf))):0.0;
  fcor_w = (focw*nw-fow*fcw) / sqrt((foow*nw-fow*fow)*(fccw*nw-fcw*fcw));
  fcor_f = (nf>0.0)?((focf*nf-fof*fcf) / sqrt((foof*nf-fof*fof)*(fccf*nf-fcf*fcf))):0.0;
  llkg_w = lw/nw;
  llkg_f = (nf>0.0)?(lf/nf):0.0;

  // done
  return true;
}

void Refine_HL_coeff::debug( const clipper::HKL_data<clipper::data32::ABCD>& abcd, const clipper::HKL_data<clipper::data32::F_sigF>& fsig, const Xmap_target_base& xtgt )
{
  // set up first cycle
  clipper::HKL_data<clipper::data32::F_phi> fphi( fsig.base_hkl_info() );
  clipper::HKL_data<clipper::data32::F_phi> fphid( fsig.base_hkl_info() );
  // set up derivative lists
  clipper::HKL_data<ArgGrad> grad( fsig.base_hkl_info() );
  clipper::HKL_data<ArgCurv> curv( fsig.base_hkl_info() );

  // first get phase and weight for first cycle
  clipper::data32::Phi_fom phifom;
  clipper::HKL_info::HKL_reference_index ih;
  for ( ih = fphi.first(); !ih.last(); ih.next() ) {
    clipper::LogPhaseProb<72> q( ih.hkl_class() );
    q.set_abcd( abcd[ih] );
    q.get_phi_fom( phifom );
    if ( !fsig[ih].missing() ) {
      fphi[ih].f()   = fsig[ih].f() * phifom.fom();
      fphi[ih].phi() = phifom.phi();
    } else {
      fphi[ih].f()   = fphi[ih].phi() = 0.0;
    }
  }

  // and get mean F^2
  std::vector<double> params( 10, 1.0 );
  clipper::BasisFn_spline basisfn( fsig, 10 );
  clipper::TargetFn_meanFnth<clipper::data32::F_sigF> targetfn( fsig, 2.0 );
  clipper::ResolutionFn meanf2( fsig.base_hkl_info(), basisfn, targetfn, params );

  // main loop
  for ( int cyc = 0; cyc < 1; cyc++ ) {
    // calc function value and derivs
    xtgt.curv( fphi, grad, curv );

    // loop over reflections and make new probability distributions
    float f, f2, fm, fw, qmap;
    double llkscale = 0.1;
    for ( ih = fphi.first(); !ih.last(); ih.next() )
     if ( ih.hkl() != clipper::HKL::zero() && ih.index() % 100 == 0 ) {
      // get prior abcd
      clipper::LogPhaseProb<72> qmap( ih.hkl_class() );
      f2 = ih.hkl_class().epsilonc() * meanf2.f(ih);
      // calculate most probable F
      fm = llk_max( fphi[ih], grad[ih], curv[ih], f2, llkscale );
      // Now pick real F, or estimate missing F
      if ( !fsig[ih].missing() ) f = fsig[ih].f();
      else                       f = fm;
      // calc R, R-free, LLK, LLK-free
      qmap = llk_dist( fphi[ih], grad[ih], curv[ih], f2, llkscale,
		       ih.hkl_class(), f );
      // get f's
      float a = fphi[ih].a();
      float b = fphi[ih].b();
      // get grad and curvature
      float da  = grad[ih].a()*llkscale;
      float db  = grad[ih].b()*llkscale;
      float daa = curv[ih].curv_aa()*llkscale;
      float dab = curv[ih].curv_ab()*llkscale;
      float dbb = curv[ih].curv_bb()*llkscale;
      // offset to get gradients at origin by subtracting current location
      da -= daa * a + dab * b;
      db -= dab * a + dbb * b;
      // now apply wilson probability distribution
      float dww = -1.0/f2;
      daa += dww;
      dbb += dww;
      // calc q's
      clipper::LogPhaseProb<72> q1( ih.hkl_class() );
      clipper::LogPhaseProb<72> q2( ih.hkl_class() );
      int dp = ( qmap.size() == 2 ) ? 1 : 12;
      for ( int p = 0; p < qmap.size(); p += dp ) {
	a = f * cos( qmap.phase(p) );
	b = f * sin( qmap.phase(p) );
	q1[p] = ( a*da + b*db + 0.5*a*a*daa + a*b*dab + 0.5*b*b*dbb );
	fphid = fphi;
	fphid[ih].f() = f;
	fphid[ih].phi() = qmap.phase(p);
	q2[p] = xtgt.func( fphid )*llkscale;
	std::cout << ih.hkl().format() << "  " << p << "   \t" << q1[p] - q1[0] << "   \t" << q2[p] - q2[0] << "   \t" << qmap[p] - qmap[0] << "\n";
      }
      if ( fsig[ih].missing() ) {
	double t, t1, t2, t3, t4;
	a = -(dbb*da-dab*db)/(daa*dbb-dab*dab);
	b = -(daa*db-dab*da)/(daa*dbb-dab*dab);
	double s1, s2, s3, s4, s;
	s = a*da + b*db + 0.5*a*a*daa + a*b*dab + 0.5*b*b*dbb;
	s1 = (a+5)*da + b*db + 0.5*(a+5)*(a+5)*daa + (a+5)*b*dab + 0.5*b*b*dbb;
	s2 = (a-5)*da + b*db + 0.5*(a-5)*(a-5)*daa + (a-5)*b*dab + 0.5*b*b*dbb;
	s3 = a*da + (b+5)*db + 0.5*a*a*daa + a*(b+5)*dab + 0.5*(b+5)*(b+5)*dbb;
	s4 = a*da + (b-5)*db + 0.5*a*a*daa + a*(b-5)*dab + 0.5*(b-5)*(b-5)*dbb;
	std::cout << "Missing " << s1 - s << " \t" << s2 - s << " \t"
                                << s3 - s << " \t" << s4 - s << "\n";
      }
    }
  }
}


MapFilterFn_shell::MapFilterFn_shell( const clipper::ftype& inner, const clipper::ftype& outer, const bool hollow )
{
  outer_ = outer;
  inner_ = (hollow)?(inner):(-inner);
}


clipper::ftype MapFilterFn_shell::operator() ( const clipper::ftype& radius ) const
{
  if      ( fabs( radius - inner_ ) <= inner_ )
    return 0.5+0.5*sin(0.5*clipper::Util::pi()*(radius-inner_)/fabs(inner_));
  else if ( fabs( radius - outer_ ) <= fabs(inner_) )
    return 0.5-0.5*sin(0.5*clipper::Util::pi()*(radius-outer_)/fabs(inner_));
  else if ( radius > inner_ && radius < outer_ )
    return 1.0;
  else
    return 0.0;
}


Map_local_moment_ordinal::Map_local_moment_ordinal( const clipper::Xmap<float>& xmap, const clipper::MapFilterFn_base& fn )
{
  // make squared map
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  clipper::Xmap<float> xmap2( xmap );
  for ( MRI ix = xmap2.first(); !ix.last(); ix.next() )
    xmap2[ix] = clipper::Util::sqr( xmap2[ix] );

  // now calculate local mom1, local mom1 squared
  clipper::MapFilter_fft<float> fltr( fn, 1.0, clipper::MapFilter_fft<float>::Relative );
  fltr( lmom1, xmap );
  fltr( lmom2, xmap2 );

  // calculate std deviation
  for ( MRI ix = lmom1.first(); !ix.last(); ix.next() )
    lmom2[ix] = sqrt( lmom2[ix] - clipper::Util::sqr( lmom1[ix] ) );

  // work out ordinals for map distributions
  std::vector<double> vals;
  for ( MRI ix = lmom1.first(); !ix.last(); ix.next() )
    vals.push_back( lmom1[ix] );
  ord_mom1.init( vals );
  vals.clear();
  for ( MRI ix = lmom2.first(); !ix.last(); ix.next() )
    vals.push_back( lmom2[ix] );
  ord_mom2.init( vals );
  vals.clear();
}


Refine_HL_simulate::Refine_HL_simulate( bool un_bias, double rad_inner, double rad_outer, double ncs_radius, double ncs_volume, double weight_llk, double weight_ramp, double skew_moment1, double skew_moment2, int nbins_moment1, int nbins_moment2, int ncyc_int, double oversampling ) { init( un_bias, rad_inner, rad_outer, ncs_radius, ncs_volume, weight_llk, weight_ramp, skew_moment1, skew_moment2, nbins_moment1, nbins_moment2, ncyc_int, oversampling ); }

void Refine_HL_simulate::init( bool un_bias, double rad_inner, double rad_outer, double ncs_radius, double ncs_volume, double weight_llk, double weight_ramp, double skew_moment1, double skew_moment2, int nbins_moment1, int nbins_moment2, int ncyc_int, double oversampling )
{
  unbias = un_bias;
  rad1 = rad_inner; rad2 = rad_outer;
  ncsrad = ncs_radius; ncsvol = ncs_volume;
  wtllk = weight_llk;
  wtrmp = weight_ramp;
  skew1  = exp(skew_moment1); skew2  = exp(skew_moment2);
  nbins1 = nbins_moment1;     nbins2 = nbins_moment2;
  ncyc = ncyc_int;
  oversam = oversampling;
  debug_mode = false;
}


bool Refine_HL_simulate::operator() (
  clipper::HKL_data<clipper::data32::F_phi>& fphi,
  clipper::HKL_data<clipper::data32::ABCD>& abcd_new,
  const clipper::HKL_data<clipper::data32::F_sigF>& fsig,
  const clipper::HKL_data<clipper::data32::F_sigF>& fobs,
  const clipper::HKL_data<clipper::data32::ABCD>& abcd,
  const clipper::HKL_data<clipper::data32::F_sigF>& ref_f,
  const clipper::HKL_data<clipper::data32::ABCD>& ref_hlcal,
  const clipper::HKL_data<clipper::data32::ABCD>& ref_hlsim,
  const std::vector<Local_rtop>& ncsops )
{
  // shared data used in both sections:
  clipper::Range<double> range_ref;
  clipper::Array2d<clipper::Histogram> hist_lib( nbins1, nbins2 );
  MapFilterFn_shell fn( rad1, rad2 );

  typedef clipper::Xmap<float>::Map_reference_index MRI;

  // FIRST STAGE:
  // compile the target histograms using the reference data
  {
    // constants
    const int n_hist = 200;     // source histogram sampling
    const int bin_min = 1000;  // minimum data in source hist

    // get source info
    const clipper::HKL_info& hkls = ref_f.base_hkl_info();

    // get map coeffs
    clipper::HKL_data<clipper::data32::Phi_fom> phiw( hkls );
    clipper::HKL_data<clipper::data32::F_phi> fcal( hkls );
    clipper::HKL_data<clipper::data32::F_phi> fsim( hkls );
    phiw.compute( ref_hlcal, clipper::data32::Compute_phifom_from_abcd() );
    fcal.compute(ref_f,phiw,clipper::data32::Compute_fphi_from_fsigf_phifom());
    phiw.compute( ref_hlsim, clipper::data32::Compute_phifom_from_abcd() );
    fsim.compute(ref_f,phiw,clipper::data32::Compute_fphi_from_fsigf_phifom());

    // make electron density map    
    clipper::Grid_sampling grid( hkls.spacegroup(), hkls.cell(), hkls.resolution()/*, oversam*/ );
    clipper::Xmap<float> xmap( hkls.spacegroup(), hkls.cell(), grid );
    xmap.fft_from( fsim );
    rms_sim = clipper::Map_stats( xmap ).std_dev();

    // calculate local stats and ordinals
    Map_local_moment_ordinal local_ord( xmap, fn );

    // now get fcalc map
    xmap.fft_from( fcal );
    rms_cal = clipper::Map_stats( xmap ).std_dev();

    // find the range of the map
    clipper::Range<double> range_map;
    for ( MRI ix = xmap.first(); !ix.last(); ix.next() )
      range_map.include( xmap[ix] );
    // cut off the most extreme 0.1%/0.01% each way
    clipper::Generic_ordinal mapord;
    mapord.init( range_map );
    for ( MRI ix = xmap.first(); !ix.last(); ix.next() )
      mapord.accumulate( xmap[ix] );
    mapord.prep_ordinal();
    mapord.invert();
    range_ref  = clipper::Range<double>( mapord.ordinal(0.0001),
					 mapord.ordinal(0.9999) );

    // now build the histograms
    clipper::Histogram hist_null( range_map, n_hist );
    hist1 = hist0 =
      clipper::Array2d<clipper::Histogram>( nbins1, nbins2, hist_null );

    // accumulate initial histograms
    for ( MRI ix = xmap.first(); !ix.last(); ix.next() )
      if ( range_map.contains( xmap[ix] ) ) {
	int i1 = clipper::Util::intf( nbins1*local_ord.ord_moment_1(ix) );
	int i2 = clipper::Util::intf( nbins2*local_ord.ord_moment_2(ix) );
	hist0(i1,i2).accumulate( xmap[ix] );
      }

    // merge data from adjacent bins if data too sparse
    for ( int i = 0; i < nbins1; i++ )
      for ( int j = 0; j < nbins2; j++ ) {
	hist1(i,j) = hist_null;
	int box = 0;
	while ( hist1(i,j).sum() < bin_min ) {
	  for ( int k = i - box; k <= i + box; k++ )
	    for ( int l = j - box; l <= j + box; l++ )
	      if ( k >= 0 && k < nbins1 && l >= 0 && l < nbins2 )
		hist1(i,j) += hist0(k,l);
	  box++;
	}
      }

    // smooth the histograms
    for ( int i = 0; i < nbins1; i++ )
      for ( int j = 0; j < nbins2; j++ ) {
	hist0(i,j) = hist1(i,j);
	hist1(i,j) = HistogramSmoother::smooth_curvature( hist0(i,j) );
      }
    //for ( int k = 0; k < n_hist; k++ ) std::cout << k << "\t" << hist0(nbins1-1,nbins2-1).y(k) << "\t" << hist1(nbins1-1,nbins2-1).y(k) << "\n";

    hist_lib = hist1;
  }

  // SECOND STAGE:
  // compile the target histograms using the reference data
  {
    // get source info
    const clipper::HKL_info& hkls = fphi.base_hkl_info();
    const clipper::Spacegroup& spgr = hkls.spacegroup();
    const clipper::Cell&       cell = hkls.cell();

    // make map
    clipper::HKL_data<clipper::data32::Phi_fom> phiw( hkls );
    phiw.compute( abcd, clipper::data32::Compute_phifom_from_abcd() );
    fphi.compute( fobs, phiw, clipper::data32::Compute_fphi_from_fsigf_phifom() );
    clipper::Grid_sampling grid( spgr, cell, hkls.resolution(), oversam );
    clipper::Xmap<float> xmap( spgr, cell, grid );
    xmap.fft_from( fphi );
    rms_wrk = clipper::Map_stats( xmap ).std_dev();

    // calculate local stats and ordinals
    Map_local_moment_ordinal local_ord( xmap, fn );

    // now do the NCS search
    NCSfind ncsref( 3.0, 0.1 );
    // refine
    double mskrad = ncsvol*pow(cell.volume()/double(spgr.num_symops()),0.333);
    ncsops_ = Local_rtop::exclude_inverse( ncsops, spgr, cell );
    
    Local_rtop rtop;
    std::vector<Local_rtop> ncsops_tmp;
    for ( int i = 0; i < ncsops_.size(); i++ ) {
      rtop = ncsref.masked_refine_from_map( xmap, ncsops_[i], mskrad, ncsrad );
      if ( !rtop.is_null() ) ncsops_tmp.push_back( rtop );
    }
    // expand
    ncsops_ = Local_rtop::include_inverse( ncsops_tmp, spgr, cell );
    ncsops_ = Local_rtop::tidy( ncsops_, spgr, cell );

    // now generate NCS restrains
    Xmap_ncs xncs;
    xncs.init( xmap.spacegroup(), xmap.cell(), xmap.grid_sampling() );
    xncs = Gaussian_probability_1d( 0.0, 1.0e6 );
    for ( int i = 0; i < ncsops_.size(); i++ )
      xncs.restrain_ncs( xmap, ncsops_[i], mskrad, ncsrad );

    // find the range of the ncs map
    clipper::Range<double> range_map, range_wrk;
    for ( MRI ix = xncs.first(); !ix.last(); ix.next() )
      range_map.include( xncs[ix].mean() );
    // cut off the most extreme 0.1%/0.01% each way
    clipper::Generic_ordinal mapord;
    mapord.init( range_map );
    for ( MRI ix = xncs.first(); !ix.last(); ix.next() )
      mapord.accumulate( xncs[ix].mean() );
    mapord.prep_ordinal();
    mapord.invert();
    range_wrk  = clipper::Range<double>( mapord.ordinal(0.001),
					 mapord.ordinal(0.999) );
    //std::cout << "Ranges: " << range_ref.min() << " " << range_ref.max() << " \t" << range_wrk.min() << " " << range_wrk.max() << " \t" << range_map.min() << " " << range_map.max() << "\n";
    // and extend overall range to include it
    range_ref.include( range_wrk.min() );
    range_ref.include( range_wrk.max() );
    range_ref.include( range_map.min() );
    range_ref.include( range_map.max() );

    clipper::Xmap<float> x1,x2;
    x1.init( xmap.spacegroup(), xmap.cell(), xmap.grid_sampling() );
    x2.init( xmap.spacegroup(), xmap.cell(), xmap.grid_sampling() );
    for ( MRI ix = xncs.first(); !ix.last(); ix.next() ) {
      x1[ix] = xncs[ix].mean();
      x2[ix] = xncs[ix].std_dev();
    }
    //clipper::CCP4MAPfile mapout;
    //mapout.open_write( "x1.map" );
    //mapout.export_xmap( x1 );
    //mapout.close_write();
    //mapout.open_write( "x2.map" );
    //mapout.export_xmap( x2 );
    //mapout.close_write();

    // make minihists from hists
    clipper::Array2d<TargetHist> hist_bins( nbins1, nbins2 );
    for ( int i = 0; i < nbins1; i++ )
      for ( int j = 0; j < nbins2; j++ )
	hist_bins(i,j) = TargetHist( range_ref, hist_lib(i,j) );

    // make the target map
    clipper::Xmap<TargetHistCompr> target( spgr, cell, xmap.grid_sampling() );
    for ( MRI ix = xmap.first(); !ix.last(); ix.next() ) {
      double x1 = local_ord.ord_moment_1(ix);
      double x2 = local_ord.ord_moment_2(ix);
      x1 = skew1*x1/(1.0-(1.0-skew1)*x1);
      x2 = skew2*x2/(1.0-(1.0-skew2)*x2);
      int i1 = clipper::Util::intf( nbins1*x1 );
      int i2 = clipper::Util::intf( nbins2*x2 );
      TargetHist h1 = hist_bins( i1, i2 );
      TargetHist h2( range_ref, xncs[ix].mean(), xncs[ix].std_dev() );
      target[ix] = TargetHistCompr( h1*h2 );
    }

    // make target fn
    Xmap_target_minihist xtarget;
    xtarget.init( target, range_ref );

    // match origin term to zero be shifting histograms
    xtarget.set_mean( 0.0 );

    // make refinement weights
    std::vector<std::pair<double,double> > wcyc(ncyc);
    for ( int i = 0; i < ncyc; i++ ) {
      double x = double( i + 1 ) / double( ncyc );
      x = clipper::Util::min( wtrmp*x, 1.0 );
      if ( !unbias ) {
	wcyc[i].first  = 1.0;
	wcyc[i].second = x;
      } else {
	wcyc[i].first = 1.0 - x;
	wcyc[i].second = x;
      }
    }

    // refine HL coeffs
    Refine_HL_coeff refine( wcyc, wtllk );
    refine( fphi, abcd_new, abcd, fsig, fobs, xtarget );

    // get results
    rfac_w = refine.r_factor_work();
    rfac_f = refine.r_factor_free();
    ecor_w = refine.e_correl_work();
    ecor_f = refine.e_correl_free();
    fcor_w = refine.f_correl_work();
    fcor_f = refine.f_correl_free();
    llkg_w = refine.llk_gain_work();
    llkg_f = refine.llk_gain_free();

    // debugging code
    if ( debug_mode ) {
      xtarget.debug( fphi );
      refine.debug( abcd, fsig, xtarget );
    }
  }

  return true;
}

bool Refine_HL_simulate::operator() (
  clipper::HKL_data<clipper::data32::F_phi>& fphi,
  clipper::HKL_data<clipper::data32::ABCD>& abcd_new,
  const clipper::HKL_data<clipper::data32::F_sigF>& fsig,
  const clipper::HKL_data<clipper::data32::F_sigF>& fobs,
  const clipper::HKL_data<clipper::data32::ABCD>& abcd,
  const std::vector<Local_rtop>& ncsops )
{
  clipper::Range<double> range_ref( -1.0, 1.0 );
  typedef clipper::Xmap<float>::Map_reference_index MRI;

  {
    // get source info
    const clipper::HKL_info& hkls = fphi.base_hkl_info();
    const clipper::Spacegroup& spgr = hkls.spacegroup();
    const clipper::Cell&       cell = hkls.cell();

    // make map
    clipper::HKL_data<clipper::data32::Phi_fom> phiw( hkls );
    phiw.compute( abcd, clipper::data32::Compute_phifom_from_abcd() );
    fphi.compute( fobs, phiw, clipper::data32::Compute_fphi_from_fsigf_phifom() );
    clipper::Grid_sampling grid( spgr, cell, hkls.resolution(), oversam );
    clipper::Xmap<float> xmap( spgr, cell, grid );
    xmap.fft_from( fphi );
    rms_wrk = clipper::Map_stats( xmap ).std_dev();

    // expand
    ncsops_ = Local_rtop::include_inverse( ncsops, spgr, cell );
    ncsops_ = Local_rtop::tidy( ncsops_, spgr, cell );

    // now generate NCS restrains
    double mskrad = ncsvol*pow(cell.volume()/double(spgr.num_symops()),0.333);
    Xmap_ncs xncs;
    xncs.init( xmap.spacegroup(), xmap.cell(), xmap.grid_sampling() );
    xncs = Gaussian_probability_1d( 0.0, 1.0e6 );
    for ( int i = 0; i < ncsops_.size(); i++ )
      xncs.restrain_ncs( xmap, ncsops_[i], mskrad, ncsrad );

    // find the range of the ncs map
    clipper::Range<double> range_map, range_wrk;
    for ( MRI ix = xncs.first(); !ix.last(); ix.next() )
      range_map.include( xncs[ix].mean() );
    // cut off the most extreme 0.1%/0.01% each way
    clipper::Generic_ordinal mapord;
    mapord.init( range_map );
    for ( MRI ix = xncs.first(); !ix.last(); ix.next() )
      mapord.accumulate( xncs[ix].mean() );
    mapord.prep_ordinal();
    mapord.invert();
    range_wrk  = clipper::Range<double>( mapord.ordinal(0.001),
					 mapord.ordinal(0.999) );
    //std::cout << "Ranges: " << range_ref.min() << " " << range_ref.max() << " \t" << range_wrk.min() << " " << range_wrk.max() << " \t" << range_map.min() << " " << range_map.max() << "\n";
    // and extend overall range to include it
    range_ref.include( range_wrk.min() );
    range_ref.include( range_wrk.max() );
    range_ref.include( range_map.min() );
    range_ref.include( range_map.max() );

    clipper::Xmap<float> x1,x2;
    x1.init( xmap.spacegroup(), xmap.cell(), xmap.grid_sampling() );
    x2.init( xmap.spacegroup(), xmap.cell(), xmap.grid_sampling() );
    for ( MRI ix = xncs.first(); !ix.last(); ix.next() ) {
      x1[ix] = xncs[ix].mean();
      x2[ix] = xncs[ix].std_dev();
    }
    //clipper::CCP4MAPfile mapout;
    //mapout.open_write( "x1.map" );
    //mapout.export_xmap( x1 );
    //mapout.close_write();
    //mapout.open_write( "x2.map" );
    //mapout.export_xmap( x2 );
    //mapout.close_write();

    // make the target map
    clipper::Xmap<TargetHistCompr> target( spgr, cell, xmap.grid_sampling() );
    for ( MRI ix = xmap.first(); !ix.last(); ix.next() ) {
      TargetHist h2( range_ref, xncs[ix].mean(), xncs[ix].std_dev() );
      target[ix] = TargetHistCompr( h2 );
    }

    // make target fn
    Xmap_target_minihist xtarget;
    xtarget.init( target, range_ref );

    // match origin term to zero be shifting histograms
    xtarget.set_mean( 0.0 );

    // make refinement weights
    std::vector<std::pair<double,double> > wcyc(ncyc);
    for ( int i = 0; i < ncyc; i++ ) {
      double x = double( i + 1 ) / double( ncyc );
      x = clipper::Util::min( wtrmp*x, 1.0 );
      if ( !unbias ) {
	wcyc[i].first  = 1.0;
	wcyc[i].second = x;
      } else {
	wcyc[i].first = 1.0 - x;
	wcyc[i].second = x;
      }
    }

    // refine HL coeffs
    Refine_HL_coeff refine( wcyc, wtllk );
    refine( fphi, abcd_new, abcd, fsig, fobs, xtarget );

    // get results
    rfac_w = refine.r_factor_work();
    rfac_f = refine.r_factor_free();
    ecor_w = refine.e_correl_work();
    ecor_f = refine.e_correl_free();
    fcor_w = refine.f_correl_work();
    fcor_f = refine.f_correl_free();
    llkg_w = refine.llk_gain_work();
    llkg_f = refine.llk_gain_free();

    // debugging code
    if ( debug_mode ) {
      xtarget.debug( fphi );
      refine.debug( abcd, fsig, xtarget );
    }
  }

  return true;
}




double Refine_HL_simulate::correl_moment_1(
  clipper::HKL_data<clipper::data32::F_phi>& fphi,
  clipper::HKL_data<clipper::data32::ABCD>& abcd_new,
  const clipper::HKL_data<clipper::data32::F_sigF>& fsig,
  const clipper::HKL_data<clipper::data32::F_sigF>& fobs,
  const clipper::HKL_data<clipper::data32::ABCD>& abcd,
  const clipper::HKL_data<clipper::data32::F_sigF>& ref_f,
  const clipper::HKL_data<clipper::data32::ABCD>& ref_hlcal,
  const clipper::HKL_data<clipper::data32::ABCD>& ref_hlsim,
  const std::vector<Local_rtop>& ncsops )
{
  // shared data used in both sections:
  clipper::Range<double> range_ref;
  clipper::Generic_ordinal ord_ref, ord_wrk;
  MapFilterFn_shell fn( rad1, rad2 );
  double correl;

  typedef clipper::Xmap<float>::Map_reference_index MRI;

  // FIRST STAGE:
  // compile the target histograms using the reference data
  {
    // constants
    const int n_hist = 200;     // source histogram sampling
    const int bin_min = 1000;  // minimum data in source hist

    // get source info
    const clipper::HKL_info& hkls = ref_f.base_hkl_info();

    // get map coeffs
    clipper::HKL_data<clipper::data32::Phi_fom> phiw( hkls );
    clipper::HKL_data<clipper::data32::F_phi> fsim( hkls );
    phiw.compute( ref_hlsim, clipper::data32::Compute_phifom_from_abcd() );
    fsim.compute( ref_f, phiw, clipper::data32::Compute_fphi_from_fsigf_phifom());
    clipper::Grid_sampling grid( hkls.spacegroup(), hkls.cell(), hkls.resolution()/*, oversam*/ );
    clipper::Xmap<float> xmap( hkls.spacegroup(), hkls.cell(), grid );
    xmap.fft_from( fsim );
    rms_sim = clipper::Map_stats( xmap ).std_dev();

    // calculate local stats and ordinals
    Map_local_moment_ordinal local_ord( xmap, fn );
    ord_ref = local_ord.ordinal_fn_moment_1();
  }
  ord_ref.invert();

  // SECOND STAGE:
  // compile the target histograms using the reference data
  {
    // get source info
    const clipper::HKL_info& hkls = fphi.base_hkl_info();
    const clipper::Spacegroup& spgr = hkls.spacegroup();
    const clipper::Cell&       cell = hkls.cell();

    // make map
    clipper::HKL_data<clipper::data32::Phi_fom> phiw( hkls );
    phiw.compute( abcd, clipper::data32::Compute_phifom_from_abcd() );
    fphi.compute( fobs, phiw, clipper::data32::Compute_fphi_from_fsigf_phifom() );
    clipper::Grid_sampling grid( spgr, cell, hkls.resolution(), oversam );
    clipper::Xmap<float> xmap( spgr, cell, grid );
    xmap.fft_from( fphi );
    rms_wrk = clipper::Map_stats( xmap ).std_dev();

    // calculate local stats and ordinals
    Map_local_moment_ordinal local_ord( xmap, fn );
    ord_wrk = local_ord.ordinal_fn_moment_1();

    // make the target map
    clipper::Xmap<TargetHistCompr> target( spgr, cell, xmap.grid_sampling() );
    double sn, sx, sy, sxx, syy, sxy;
    sn = sx = sy = sxx = syy = sxy = 0.0;
    for ( MRI ix = xmap.first(); !ix.last(); ix.next() ) {
      double x = local_ord.local_moment_1()[ix];
      double o1 = ord_wrk.ordinal(x);
      double o2 = skew1*o1/(1.0-(1.0-skew1)*o1);
      double y = ord_ref.ordinal( o2 );
      sn += 1.0;
      sx += x;
      sy += y;
      sxx += x*x;
      syy += y*y;
      sxy += x*y;
    }
    correl = (sn*sxy-sx*sy)/sqrt((sn*sxx-sx*sx)*(sn*syy-sy*sy));
  }
  return correl;
}

double Refine_HL_simulate::correl_moment_2(
  clipper::HKL_data<clipper::data32::F_phi>& fphi,
  clipper::HKL_data<clipper::data32::ABCD>& abcd_new,
  const clipper::HKL_data<clipper::data32::F_sigF>& fsig,
  const clipper::HKL_data<clipper::data32::F_sigF>& fobs,
  const clipper::HKL_data<clipper::data32::ABCD>& abcd,
  const clipper::HKL_data<clipper::data32::F_sigF>& ref_f,
  const clipper::HKL_data<clipper::data32::ABCD>& ref_hlcal,
  const clipper::HKL_data<clipper::data32::ABCD>& ref_hlsim,
  const std::vector<Local_rtop>& ncsops )
{
  // shared data used in both sections:
  clipper::Range<double> range_ref;
  clipper::Generic_ordinal ord_ref, ord_wrk;
  MapFilterFn_shell fn( rad1, rad2 );
  double correl;

  typedef clipper::Xmap<float>::Map_reference_index MRI;

  // FIRST STAGE:
  // compile the target histograms using the reference data
  {
    // constants
    const int n_hist = 200;     // source histogram sampling
    const int bin_min = 1000;  // minimum data in source hist

    // get source info
    const clipper::HKL_info& hkls = ref_f.base_hkl_info();

    // get map coeffs
    clipper::HKL_data<clipper::data32::Phi_fom> phiw( hkls );
    clipper::HKL_data<clipper::data32::F_phi> fsim( hkls );
    phiw.compute( ref_hlsim, clipper::data32::Compute_phifom_from_abcd() );
    fsim.compute( ref_f, phiw, clipper::data32::Compute_fphi_from_fsigf_phifom());
    clipper::Grid_sampling grid( hkls.spacegroup(), hkls.cell(), hkls.resolution()/*, oversam*/ );
    clipper::Xmap<float> xmap( hkls.spacegroup(), hkls.cell(), grid );
    xmap.fft_from( fsim );
    rms_sim = clipper::Map_stats( xmap ).std_dev();

    // calculate local stats and ordinals
    Map_local_moment_ordinal local_ord( xmap, fn );
    ord_ref = local_ord.ordinal_fn_moment_2();
  }
  ord_ref.invert();

  // SECOND STAGE:
  // compile the target histograms using the reference data
  {
    // get source info
    const clipper::HKL_info& hkls = fphi.base_hkl_info();
    const clipper::Spacegroup& spgr = hkls.spacegroup();
    const clipper::Cell&       cell = hkls.cell();

    // make map
    clipper::HKL_data<clipper::data32::Phi_fom> phiw( hkls );
    phiw.compute( abcd, clipper::data32::Compute_phifom_from_abcd() );
    fphi.compute( fobs, phiw, clipper::data32::Compute_fphi_from_fsigf_phifom() );
    clipper::Grid_sampling grid( spgr, cell, hkls.resolution(), oversam );
    clipper::Xmap<float> xmap( spgr, cell, grid );
    xmap.fft_from( fphi );
    rms_wrk = clipper::Map_stats( xmap ).std_dev();

    // calculate local stats and ordinals
    Map_local_moment_ordinal local_ord( xmap, fn );
    ord_wrk = local_ord.ordinal_fn_moment_2();

    // make the target map
    clipper::Xmap<TargetHistCompr> target( spgr, cell, xmap.grid_sampling() );
    double sn, sx, sy, sxx, syy, sxy;
    sn = sx = sy = sxx = syy = sxy = 0.0;
    for ( MRI ix = xmap.first(); !ix.last(); ix.next() ) {
      double x = local_ord.local_moment_2()[ix];
      double o1 = ord_wrk.ordinal(x);
      double o2 = skew2*o1/(1.0-(1.0-skew2)*o1);
      double y = ord_ref.ordinal( o2 );
      sn += 1.0;
      sx += x;
      sy += y;
      sxx += x*x;
      syy += y*y;
      sxy += x*y;
    }
    correl = (sn*sxy-sx*sy)/sqrt((sn*sxx-sx*sx)*(sn*syy-sy*sy));
  }
  return correl;
}

