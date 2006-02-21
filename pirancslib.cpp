/*! \file pirancslib.cpp pirate library */
/* Copyright 2003-2005 Kevin Cowtan & University of York all rights reserved */

//L   This code is distributed under the terms and conditions of the
//L   CCP4 Program Suite Licence Agreement as a CCP4 Application.
//L   A copy of the CCP4 licence can be obtained by writing to the
//L   CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.

#include "pirancslib.h"


struct Compare_pair_score_rotation{ bool operator() ( const std::pair<double,clipper::Rotation>& r1, const std::pair<double,clipper::Rotation>& r2 ) const { return ( r1.first < r2.first ); } };


Optimiser_simplex::Optimiser_simplex( double tolerance, int max_cycles, TYPE type )
{
  tolerance_ = tolerance;
  max_cycles_ = max_cycles;
  type_ = type;
  debug_mode = false;
}


std::vector<double> Optimiser_simplex::operator() ( const Target_fn_order_zero& target_fn, const std::vector<std::vector<double> >& args ) const
{
  // check input
  int size = target_fn.num_params();
  if ( args.size() != size+1 )
    clipper::Message::message( clipper::Message_fatal( "Optimiser_simplex: parameter size mismatch" ) );
  for ( int i = 0; i < args.size(); i++ )
    if ( args[i].size() != size )
      clipper::Message::message( clipper::Message_fatal( "Optimiser_simplex: parameter size mismatch" ) );

  // make arrays
  std::vector<std::vector<double> > params( size + 1 );
  std::vector<double> fn( size + 1 );
  // calc target fn
  for ( int i = 0; i < args.size(); i++ ) {
    params[i] = args[i];
    fn[i] = target_fn( args[i] );
  }

  // simplex loop
  int iw, ib;
  double f1(0.0), f2(0.0);
  for ( int cyc = 0; cyc < max_cycles_; cyc++ ) {
    if ( debug_mode )  // DEBUG OUTPUT
      for ( int i = 0; i < params.size(); i++ ) {
	std::cout << i << "  " << clipper::String( fn[i] ) << "  ";
	for ( int j = 0; j < size; j++ )
	  std::cout << clipper::String( params[i][j] ) << "\t";
	std::cout << "\n";
      }
    // find worst point
    iw = ib = 0;
    for ( int i = 0; i < params.size(); i++ ) {
      if ( fn[i] > fn[iw] ) iw = i;
      if ( fn[i] < fn[ib] ) ib = i;
    }
    // termination condition
    if ( fn[iw] - fn[ib] < tolerance_ ) break;
    // find centroid of the rest
    std::vector<double> centr( size, 0.0 ),
      shift( size ), t0( size ), t1( size ), t2( size );
    double sumweight = 0.0;
    for ( int i = 0; i < params.size(); i++ )
      if ( i != iw ) {
	double weight = 1.0;
	if ( type_ == GRADIENT ) {
	  double r2 = 0.0;
	  for ( int j = 0; j < size; j++ )
	    r2 += clipper::Util::sqr( params[i][j] - params[iw][j] );
	  weight = ( fn[iw] - fn[i] ) / sqrt( r2 );
	}
	for ( int j = 0; j < size; j++ )
	  centr[j] += weight * params[i][j];
	sumweight += weight;
      }
    for ( int j = 0; j < size; j++ ) centr[j] = centr[j] / sumweight;
    for ( int j = 0; j < size; j++ ) shift[j] = centr[j] - params[iw][j];
    // calculate first trial point
    for ( int j = 0; j < size; j++ ) {
      t0[j] = centr[j] - 0.5 * shift[j];
      t1[j] = centr[j] + 1.0 * shift[j];
      t2[j] = centr[j] + 2.0 * shift[j];
    }
    f1 = target_fn( t1 );
    // simplex conditions
    if ( !clipper::Util::is_nan( f1 ) ) {
      if ( f1 < fn[iw] ) {  // new point is better than worst
	if ( f1 < fn[ib] ) {  // new point is better than best
	  f2 = target_fn( t2 );
	  if ( !clipper::Util::is_nan( f2 ) ) {
	    if ( f2 < f1 ) {  // extended step is best
	      params[iw] = t2; fn[iw] = f2;         // extended step
	    } else {  // normal step better than extended step
	      params[iw] = t1; fn[iw] = f1;         // normal step
	    }
	  } else {  // extended step is invalid
	    params[iw] = t1; fn[iw] = f1;           // normal step
	  }
	} else {  // normal step is neither best nor worst (default)
	  params[iw] = t1; fn[iw] = f1;             // normal step
	}
      } else {  // normal step is worse
	params[iw] = t0; fn[iw] = target_fn( t0 );  // contraction
      }
    } else {  // normal step is invalid
      params[iw] = t0; fn[iw] = target_fn( t0 );    // contraction
    }
    if ( debug_mode )  // DEBUG OUTPUT
      if      ( fn[iw] == f2 ) std::cout << "Extn-step\n";
      else if ( fn[iw] == f1 ) std::cout << "Nrml-step\n";
      else                     std::cout << "Ctrn-step\n";
  }
  return params[ib];
}


NXmap_cube::NXmap_cube( const clipper::Xmap<float>& xmap, const clipper::Coord_orth& origin, const double& spacing, const double& extent )
{
  // basic properties
  origin_  = origin;
  spacing_ = spacing;
  extent_  = extent;

  // calc NXmap properties
  double gby2 = rint( 0.5 * extent_ / spacing_ );
  double g = 2.0 * gby2;
  rx_ = g / extent;
  tx_ = g / 2.0;
  int ngrid = 2 * int(gby2) + 1;
  clipper::Grid grid( ngrid, ngrid, ngrid );

  clipper::Mat33<> mat( rx_, 0.0, 0.0, 0.0, rx_, 0.0, 0.0, 0.0, rx_ );
  clipper::Vec3<>  trn( tx_, tx_, tx_ );

  clipper::RTop<> rtop( mat, trn );

  // make the nxmap
  clipper::NXmap<float>::init( grid, rtop );
  clipper::NXmap<float>& nxmap_ = *this;

  // and fill it
  typedef clipper::NXmap<double>::Map_reference_index MRI;
  clipper::Coord_orth co;
  clipper::Coord_frac cf;
  for ( MRI ix = nxmap_.first(); !ix.last(); ix.next() ) {
    // Note: can optimise this transformation to a set of grid steps
    cf = ( ix.coord_orth() + origin ).coord_frac( xmap.cell() );
    nxmap_[ix] = xmap.interp<clipper::Interp_cubic>( cf );
  }
}


Local_rtop Local_rtop::transform( const clipper::RTop_orth& r1, const clipper::RTop_orth& r2 ) const
{
  Local_rtop res;
  res.src() = r1 * src();
  res.tgt() = r2 * tgt();
  res.rot() = clipper::Rotation(r2.rot() * rot().matrix() * r1.rot().inverse());
  return res;
}


std::pair<double,Local_rtop> Local_rtop::symm_match( const Local_rtop& other, const clipper::Spacegroup& spgr, const clipper::Cell& cell, const double& tol_dst, const double& tol_ang ) const
{
  // prepare orthogonal symops
  std::vector<clipper::RTop_orth> symop_orth( spgr.num_symops() );
  for ( int i = 0; i < symop_orth.size(); i++ )
    symop_orth[i] = spgr.symop(i).rtop_orth(cell);

  // find scored symmetry match between two ops
  double ang_rad = 10.0;
  Local_rtop rcbest;
  double scbest = 1.0e9;
  clipper::Coord_frac cf;
  for ( int sym1 = 0; sym1 < spgr.num_symops(); sym1++ )
    for ( int sym2 = 0; sym2 < spgr.num_symops(); sym2++ ) {
      Local_rtop rc1 = other;
      Local_rtop rc2 = transform( symop_orth[sym1], symop_orth[sym2] );
      //rc2.src() = rc2.src().coord_frac(cell).lattice_copy_near(rc1.src().coord_frac(cell)).coord_orth(cell);
      //rc2.tgt() = rc2.tgt().coord_frac(cell).lattice_copy_near(rc1.tgt().coord_frac(cell)).coord_orth(cell);
      cf = ( rc2.src() - rc1.src() ).coord_frac(cell);
      cf = clipper::Coord_frac(rint(cf.u()),rint(cf.v()),rint(cf.w()));
      rc2.src() = rc2.src() - cf.coord_orth(cell);
      cf = ( rc2.tgt() - rc1.tgt() ).coord_frac(cell);
      cf = clipper::Coord_frac(rint(cf.u()),rint(cf.v()),rint(cf.w()));
      rc2.tgt() = rc2.tgt() - cf.coord_orth(cell);
      clipper::Rotation rot = rc1.rot().inverse() * rc2.rot();
      double a = 2.0*acos(rot.w());
      double r = sqrt( ( rc2.rtop_orth()*rc1.src() - rc1.tgt() ).lengthsq() );
      double a2 = (a*a)/(tol_ang*tol_ang);
      double r2 = (r*r)/(tol_dst*tol_dst);
      double s2 = r2 + a2;
      if ( s2 < scbest ) {
	scbest = s2;
	rcbest = rc2;
      }
    }
  return std::pair<double,Local_rtop>( scbest, rcbest );
}

std::pair<int,double> Local_rtop::atom_match( const clipper::Spacegroup& spgr, const clipper::Cell& cell, const std::vector<clipper::Coord_orth>& coords, const double& tol ) const
{
  // prepare orthogonal symops
  std::vector<clipper::RTop_orth> symop_orth( spgr.num_symops() );
  for ( int i = 0; i < symop_orth.size(); i++ )
    symop_orth[i] = spgr.symop(i).rtop_orth(cell);

  // score how many atoms are related by this operator
  double radsq = tol*tol;
  clipper::RTop_orth rtop = rtop_orth();
  clipper::Coord_orth a1, a2;
  clipper::Coord_frac cf;
  int nmatch = 0;
  double smatch = 0.0;
  for ( int atm1 = 0; atm1 < coords.size(); atm1++ ) {
    double r2min = 1.0e9;
    for ( int sym1 = 0; sym1 < spgr.num_symops(); sym1++ ) {
      a1 = symop_orth[sym1] * coords[atm1];
      //a1 = a1.coord_frac(cell).lattice_copy_near(src().coord_frac(cell)).coord_orth(cell);
      cf = (a1-src()).coord_frac(cell);
      cf = clipper::Coord_frac(rint(cf.u()),rint(cf.v()),rint(cf.w()));
      a1 = a1 - cf.coord_orth(cell);
      for ( int atm2 = 0; atm2 < coords.size(); atm2++ )
	for ( int sym2 = 0; sym2 < spgr.num_symops(); sym2++ ) {
	  a2 = symop_orth[sym2] * coords[atm2];
	  //a2 = a2.coord_frac(cell).lattice_copy_near(tgt().coord_frac(cell)).coord_orth(cell);
	  cf = (a2-tgt()).coord_frac(cell);
	  cf = clipper::Coord_frac(rint(cf.u()),rint(cf.v()),rint(cf.w()));
	  a2 = a2 - cf.coord_orth(cell);
	  double r2 = (rtop*a1 - a2).lengthsq();
	  if ( r2 < r2min ) r2min = r2;
	}
    }
    nmatch += ( r2min < radsq ) ? 1 : 0;
    smatch += exp( -r2min / radsq );
  }
  return std::pair<int,double>(nmatch,smatch);
}

std::pair<int,int> Local_rtop::atom_loop( const clipper::Spacegroup& spgr, const clipper::Cell& cell, const std::vector<clipper::Coord_orth>& coords, const double& tol ) const
{
  // prepare orthogonal symops
  std::vector<clipper::RTop_orth> symop_orth( spgr.num_symops() );
  for ( int i = 0; i < symop_orth.size(); i++ )
    symop_orth[i] = spgr.symop(i).rtop_orth(cell);

  // identify proper ncs loops related by this operator
  double radsq = tol*tol;
  int nncs, nmis, maxncs = 20;
  clipper::RTop_orth rtop = rtop_orth();
  clipper::Coord_orth next = tgt();
  nncs = nmis = 0;
  for ( nncs = 2; nncs < maxncs; nncs++ ) {
    next = rtop * next;
    if ( (next-src()).lengthsq() < radsq ) break;
    clipper::Coord_orth trial, best = coords[0];
    for ( int atm = 0; atm < coords.size(); atm++ )
      for ( int sym = 0; sym < spgr.num_symops(); sym++ ) {
	trial = symop_orth[sym] * coords[atm];
	//trial = trial.coord_frac(cell).lattice_copy_near(next.coord_frac(cell).coord_orth(cell);
	clipper::Coord_frac cf = (trial-next).coord_frac(cell);
	cf = clipper::Coord_frac(rint(cf.u()),rint(cf.v()),rint(cf.w()));
	trial = trial - cf.coord_orth(cell);
	if ( (next-trial).lengthsq() < (next-best).lengthsq() )
	  best = trial;
      }
    if ( (next-best).lengthsq() > radsq ) { best = next; nmis++; }
  }
  if ( nncs == maxncs ) nncs = 0;
  return std::pair<int,int>(nncs,nncs-nmis);
}

Local_rtop Local_rtop::proper( const clipper::Spacegroup& spgr, const clipper::Cell& cell ) const
{
  // prepare orthogonal symops
  std::vector<clipper::RTop_orth> symop_orth( spgr.num_symops() );
  for ( int i = 0; i < symop_orth.size(); i++ )
    symop_orth[i] = spgr.symop(i).rtop_orth(cell);

  // most proper version
  int maxncs = 20;
  double wloop( 1.0 ), wscrw( 0.001 ), wtran( 0.001 );
  // make list of candidate operators
  std::vector<std::pair<double,Local_rtop> > resultsym;
  clipper::RTop_orth rtid = clipper::RTop_orth::identity();
  for ( int j = 0; j < spgr.num_symops(); j++ ) {
    Local_rtop rsym, rcel;
    rsym = transform( rtid, symop_orth[j] );
    clipper::Coord_frac cf = (rsym.tgt()-rsym.src()).coord_frac(cell);
    clipper::Coord_frac df( floor(cf.u()), floor(cf.v()), floor(cf.w()) );
    for ( double u = df.u(); u <= df.u()+1.01; u += 1.0 )
      for ( double v = df.v(); v <= df.v()+1.01; v += 1.0 )
	for ( double w = df.w(); w <= df.w()+1.01; w += 1.0 ) {
	  rcel = rsym;
	  rcel.tgt() -= clipper::Coord_frac(u,v,w).coord_orth(cell);
	  resultsym.push_back(std::pair<double,Local_rtop>(0.0,rcel));
	}
  }
  // score them for properness
  for ( int k = 0; k < resultsym.size(); k++ ) {
    clipper::Mat33<> mat = resultsym[k].second.rot().matrix();
    clipper::Vec3<> v0( mat(0,0) - 1.0, mat(1,0), mat(2,0) );
    clipper::Vec3<> v1( mat(0,1), mat(1,1) - 1.0, mat(2,1) );
    clipper::Vec3<> v2( mat(0,2), mat(1,2), mat(2,2) - 1.0 );
    clipper::Vec3<> v3 = clipper::Vec3<>::cross( v1, v2 );
    clipper::Vec3<> v4 = clipper::Vec3<>::cross( v0, v2 );
    clipper::Vec3<> v5 = clipper::Vec3<>::cross( v0, v1 );
    if ( v3*v3 > v4*v4 && v3*v3 > v5*v5 )
      v0 = v3.unit();
    else if ( v4*v4 > v5*v5 )
      v0 = v4.unit();
    else
      v0 = v5.unit();
    double scrwshift = clipper::Util::sqr( v0 * resultsym[k].second.rtop_orth().trn() );
    double transhift = ( resultsym[k].second.tgt() - resultsym[k].second.src() ).lengthsq();
    double loopshift = 1.0e9;
    clipper::RTop_orth rtop;
    clipper::Coord_orth init, next;
    init = next = resultsym[k].second.src();
    for ( int r = 0; r < maxncs; r++ ) {
      next = rtop * next;
      loopshift = clipper::Util::min( loopshift, (next-init).lengthsq() );
    }
    resultsym[k].first = wloop*loopshift+wscrw*scrwshift+wtran*transhift;
  }
  // and pick the best
  std::sort( resultsym.begin(), resultsym.end() );
  return resultsym[0].second;
}


Target_fn_nxmap_rotation::Target_fn_nxmap_rotation( const NXmap_cube& source, const NXmap_cube& target, const double& rot_step, const double& trn_step, const TYPE type )
{
  source_ = &source;
  target_ = &target;
  rot_step_ = rot_step;
  trn_step_ = trn_step;
  type_ = type;
}


double Target_fn_nxmap_rotation::operator() ( const clipper::Rotation& rot ) const
{
  const NXmap_cube& source( *source_ );
  const NXmap_cube& target( *target_ );
  clipper::Mat33<> mat = rot.matrix();
  double step = trn_step_;
  double radius = 0.5 * clipper::Util::min( source.extent(), target.extent() );
  double radlim = step * floor( radius/step );
  double x, y, z;
  double score, r1, r2;
  double sn, sd, s1, s2, s11, s22, s12;
  clipper::Coord_orth c1, c2;
  sn = sd = s1 = s2 = s11 = s22 = s12 = 0.0;
  for ( z = -radlim; z <= radius; z += step )
    for ( y = -radlim; y <= radius; y += step )
      for ( x = -radlim; x <= radius; x += step )
        if ( x*x+y*y+z*z <= radius*radius ) {
          c1 = clipper::Coord_orth( x, y, z );
          c2 = clipper::Coord_orth( mat * c1 );
          r1 = source[c1];
          r2 = target[c2];
          sn += 1.0;
          sd += clipper::Util::sqr( r1 - r2 );
          s1 += r1;
          s2 += r2;
          s11 += r1*r1;
          s22 += r2*r2;
          s12 += r1*r2;
        }
  if ( type_ == Target_fn_nxmap_rotation::CORREL ) {
    score = -(sn*s12-s1*s2) / sqrt((sn*s11-s1*s1)*(sn*s22-s2*s2));
  } else if ( type_ == Target_fn_nxmap_rotation::RMSD ) {
    score = sqrt( sd / sn );
  } else {
    score = -1.0;
  }
  return score; 
}


double Target_fn_nxmap_rotation::operator() ( const std::vector<double>& args ) const
{
  return (*this)( rotation( args ) );
}


clipper::Rotation Target_fn_nxmap_rotation::rotation( const std::vector<double>& args ) const
{
  return clipper::Rotation(clipper::Euler<clipper::Rotation::EulerXYZs>(args[0],args[1],args[2])) * rot_;
}


clipper::Rotation Target_fn_nxmap_rotation::refine( const clipper::Rotation& rot )
{
  // store initial rotation
  rot_ = rot;

  // calculate initial params
  std::vector<std::vector<double> > args;
  std::vector<double> arg(3);
  double step = 0.5 * rot_step_;
  clipper::Euler<clipper::Rotation::EulerXYZs> euler( 0.0, 0.0, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( step, 0.0, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( 0.0, step, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( 0.0, 0.0, step );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );

  // simple refinement
  Optimiser_simplex os( 0.0001, 50 );
  return rotation( os( *this, args ) );
}


Search_nxmap_rotation::Search_nxmap_rotation( const double& rot_step, const double& trn_step, const Target_fn_nxmap_rotation::TYPE type )
{
  rot_step_ = rot_step;
  trn_step_ = trn_step;
  type_ = type;
}


std::vector<std::pair<double,clipper::Rotation> > Search_nxmap_rotation::operator() ( const NXmap_cube& source, const NXmap_cube& target ) const
{
  std::vector<std::pair<double,clipper::Rotation> > result;

  // make a list of rotation angles
  {
    clipper::Rotation rot;
    double pi = clipper::Util::pi();
    double step = rot_step_;
    double alpha, beta, gamma, thpl, thmi, spl, smi, opl, omi;
    opl = 0.25;
    for ( beta = 0.5*step; beta <= pi; beta += step ) {
      opl = 1.0 - opl;
      // step along theta+/-
      spl = 2.0*pi/ceil(2.0*pi*0.9428*cos(0.5*beta)/step);
      smi = 2.0*pi/ceil(2.0*pi*0.8165*sin(0.5*beta)/step);
      omi = 0.25;
      for ( thpl = opl*spl; thpl < 4.0*pi; thpl += spl ) {
	omi = 1.0 - omi;
	for ( thmi = omi*smi; thmi < 2.0*pi; thmi += smi ) {
	  // convert to alpha/gamma
	  alpha = clipper::Util::mod(0.5*(thpl+thmi),2.0*pi);
	  gamma = clipper::Util::mod(0.5*(thpl-thmi),2.0*pi);
	  rot = clipper::Rotation(clipper::Euler_ccp4(alpha,beta,gamma));
	  result.push_back( std::pair<double,clipper::Rotation>( 0.0, rot ) );
	}
      }
    }
  }

  // now compare the two maps
  Target_fn_nxmap_rotation tf( source, target, rot_step_, trn_step_, type_ );
  {
    const Search_nxmap_rotation& mr = *this;
    for ( int r = 0; r < result.size(); r++ ) 
      result[r].first = tf( result[r].second );
  }

  // make z-score params
  double sm, sd;
  {
    double sn(result.size()), s1(0.0), s2(0.0);
    for ( int r = 0; r < result.size(); r++ ) {
      s1 += result[r].first;
      s2 += result[r].first*result[r].first;
    }
    sm = s1 / sn;
    sd = sqrt( sn*s2 - s1*s1 ) / sn;
  }

  // copy the best 5%
  std::sort( result.begin(), result.end(), Compare_pair_score_rotation() );
  result.resize( result.size()/20 );

  // eliminate any points with higher neighbours (peak search)
  {
    std::vector<std::pair<double,clipper::Rotation> > rescut = result;
    result.clear();
    result.push_back( rescut[0] );
    double step = rot_step_;
    for ( int i = 1; i < rescut.size(); i++ ) {
      int j;
      double diff;
      for ( j = 0; j < i; j++ ) {
	diff = 2.0 * acos( (rescut[j].second.inverse()*rescut[i].second).w() );
	if ( diff < 3.0 * step ) break;
      }
      if ( j == result.size() ) result.push_back( rescut[i] );
    }
  }

  // refine the best results and resort
  for ( int i = 0; i < result.size(); i++ )
    if ( i < 10 ) {
      result[i].second = tf.refine( result[i].second );
      result[i].first  = tf( result[i].second );
    }
  std::sort( result.begin(), result.end(), Compare_pair_score_rotation() );

  // convert to z-score
  for ( int r = 0; r < result.size(); r++ )
    result[r].first = ( result[r].first - sm ) / sd;

  // return the result
  return result;
}


Target_fn_coord_rotation::Target_fn_coord_rotation( const Gaussian_orth_list& source, const Gaussian_orth_list& target, const double& rot_step )
{
  source_ = &source;
  target_ = &target;
  rot_step_ = rot_step;
}

double Target_fn_coord_rotation::operator() ( const clipper::Rotation& rot ) const
{
  const Gaussian_orth_list& source = *source_;
  const Gaussian_orth_list& target = *target_;
  clipper::Mat33<> mat = rot.matrix();
  double r2, result = 0.0;
  for ( int s = 0; s < source.size(); s++ ) {
    clipper::Coord_orth source_coord( mat * source[s].coord_orth() );
    for ( int t = 0; t < target.size(); t++ ) {
      r2 = ( source_coord - target[t].coord_orth() ).lengthsq();
      result += source[s].evaluate(0.5*r2) * target[t].evaluate(0.5*r2);
    }
  }
  return -result;
}

double Target_fn_coord_rotation::operator() ( const std::vector<double>& args ) const
{
  return (*this)( rotation( args ) );
}

clipper::Rotation Target_fn_coord_rotation::rotation( const std::vector<double>& args ) const
{
  return clipper::Rotation(clipper::Euler<clipper::Rotation::EulerXYZs>(args[0],args[1],args[2])) * rot_;
}


clipper::Rotation Target_fn_coord_rotation::refine( const clipper::Rotation& rot )
{
  // store initial rotation
  rot_ = rot;

  // calculate initial params
  std::vector<std::vector<double> > args;
  std::vector<double> arg(3);
  double step = 0.5 * rot_step_;
  clipper::Euler<clipper::Rotation::EulerXYZs> euler( 0.0, 0.0, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( step, 0.0, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( 0.0, step, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( 0.0, 0.0, step );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );

  // simple refinement
  Optimiser_simplex os( 0.0001, 50 );
  return rotation( os( *this, args ) );
}


Search_coord_rotation::Search_coord_rotation( const double& rot_step, const int& max_match )
{
  rot_step_ = rot_step;
  max_match_ = max_match;
}

std::vector<std::pair<double,clipper::Rotation> > Search_coord_rotation::operator() ( const Gaussian_orth_list& source, const Gaussian_orth_list& target ) const
{
  std::vector<std::pair<double,clipper::Rotation> > result;

  // limit maximum number of atoms
  int nmatch  = clipper::Util::min( source.size(), target.size() );
  nmatch = clipper::Util::min( nmatch, max_match_ );

  // first make a list of rotation candidates
  std::vector<clipper::Coord_orth> v1(3), v2(3);
  // loop over two more source atoms
  for ( int s1 = 1; s1 < nmatch-1; s1++ ) {
    double ls1 = sqrt(source[s1].coord_orth().lengthsq());
    double ws1 = source[s1].half_width();
    for ( int s2 = s1+1; s2 < nmatch; s2++ ) {
      double ls2 = sqrt(source[s2].coord_orth().lengthsq());
      double ws2 = source[s2].half_width();
      // find two appropriate target atoms
      for ( int t1 = 1; t1 < nmatch; t1++ ) {
	double lt1 = sqrt(target[t1].coord_orth().lengthsq());
	double wt1 = target[t1].half_width();
	if ( fabs(lt1-ls1) < 1.0*(ws1+wt1) ) {
	  for ( int t2 = 1; t2 < nmatch; t2++ ) {
	    double lt2 = sqrt(target[t2].coord_orth().lengthsq());
	    double wt2 = target[t2].half_width();
	    if ( fabs(lt2-ls2) < 1.0*(ws2+wt2) ) {
	      double ls12 = sqrt( (source[s2].coord_orth()-source[s1].coord_orth()).lengthsq() );
	      double lt12 = sqrt( (target[t2].coord_orth()-target[t1].coord_orth()).lengthsq() );
	      if ( fabs(ls12-lt12) < 0.5*(ws1+ws2+wt1+wt2) ) {
		v1[0] = source[0].coord_orth();
		v1[1] = source[s1].coord_orth();
		v1[2] = source[s2].coord_orth();
		v2[0] = target[0].coord_orth();
		v2[1] = target[t1].coord_orth();
		v2[2] = target[t2].coord_orth();
		result.push_back( std::pair<double,clipper::Rotation>( 0.0,
		  clipper::Rotation( clipper::RTop_orth( v1, v2 ).rot() ) ) );
	      }
	    }
	  }
	}
      }
    }
  }

  // now compare the two coord sets
  Target_fn_coord_rotation tf( source, target, rot_step_ );
  {
    const Search_coord_rotation& mr = *this;
    for ( int r = 0; r < result.size(); r++ ) 
      result[r].first = tf( result[r].second );
  }

  // make z-score params
  double sm, sd;
  {
    double sn(result.size()), s1(0.0), s2(0.0);
    for ( int r = 0; r < result.size(); r++ ) {
      s1 += result[r].first;
      s2 += result[r].first*result[r].first;
    }
    sm = s1 / sn;
    sd = sqrt( sn*s2 - s1*s1 ) / sn;
  }

  // copy the best 5%
  std::sort( result.begin(), result.end(), Compare_pair_score_rotation() );
  result.resize( result.size()/20 );

  // refine the best results and resort
  for ( int i = 0; i < result.size(); i++ )
    if ( i < 10 ) {
      result[i].second = tf.refine( result[i].second );
      result[i].first  = tf( result[i].second );
    }
  std::sort( result.begin(), result.end(), Compare_pair_score_rotation() );

  // convert to z-score
  for ( int r = 0; r < result.size(); r++ )
    result[r].first = ( result[r].first - sm ) / sd;

  // return the result
  return result;
}


// target function for NCS operator based on sphere of density

Target_fn_xmap_rtop::Target_fn_xmap_rtop( const clipper::Xmap<float>& source, const clipper::Xmap<float>& target, const double& rad, const double& rot_step, const double& trn_step, const TYPE type )
{
  source_ = &source;
  target_ = &target;
  rad_ = rad;
  rot_step_ = rot_step;
  trn_step_ = trn_step;
  type_ = type;
}


double Target_fn_xmap_rtop::operator() ( const Local_rtop& rot ) const
{
  typedef clipper::Interp_cubic InterpType;
  const clipper::Xmap<float>& source( *source_ );
  const clipper::Xmap<float>& target( *target_ );

  clipper::Mat33<> mat = rot.rot().matrix();
  double step = trn_step_;
  double radius = rad_;
  double radlim = step * floor( radius/step );
  double x, y, z;
  double score, r1, r2;
  double sn, sd, s1, s2, s11, s22, s12;
  clipper::Vec3<> c0;
  clipper::Coord_orth c1, c2;
  sn = sd = s1 = s2 = s11 = s22 = s12 = 0.0;
  for ( z = -radlim; z <= radius; z += step )
    for ( y = -radlim; y <= radius; y += step )
      for ( x = -radlim; x <= radius; x += step )
        if ( x*x+y*y+z*z <= radius*radius ) {
	  c0 = clipper::Vec3<>( x, y, z );
          c1 = clipper::Coord_orth( c0 + rot.src() );
          c2 = clipper::Coord_orth( mat*c0 + rot.tgt() );
          r1 = source.interp<InterpType>( c1.coord_frac( source.cell() ) );
          r2 = target.interp<InterpType>( c2.coord_frac( target.cell() ) );
          sn += 1.0;
          sd += clipper::Util::sqr( r1 - r2 );
          s1 += r1;
          s2 += r2;
          s11 += r1*r1;
          s22 += r2*r2;
          s12 += r1*r2;
        }
  if ( type_ == Target_fn_xmap_rtop::CORREL ) {
    score = -(sn*s12-s1*s2) / sqrt((sn*s11-s1*s1)*(sn*s22-s2*s2));
  } else if ( type_ == Target_fn_xmap_rtop::RMSD ) {
    score = sqrt( sd / sn );
  } else {
    score = -1.0;
  }
  return score; 
}


double Target_fn_xmap_rtop::operator() ( const std::vector<double>& args ) const
{
  return (*this)( local_rtop( args ) );
}


Local_rtop Target_fn_xmap_rtop::local_rtop( const std::vector<double>& args ) const
{
  return Local_rtop( clipper::Rotation(clipper::Euler<clipper::Rotation::EulerXYZs>(args[0],args[1],args[2])) * rot_.rot(), rot_.src(), clipper::Coord_orth(args[3],args[4],args[5]) + rot_.tgt() );
}


Local_rtop Target_fn_xmap_rtop::refine( const Local_rtop& rot )
{
  // store initial rotation
  rot_ = rot;

  // calculate initial params
  std::vector<std::vector<double> > args;
  std::vector<double> arg(6,0.0);
  // identity
  clipper::Euler<clipper::Rotation::EulerXYZs> euler( 0.0, 0.0, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );
  // rotation steps
  double step = 0.5 * rot_step_;
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( step, 0.0, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( 0.0, step, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( 0.0, 0.0, step );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );
  // translation steps
  step = 0.5 * trn_step_;
  arg = args[0];
  arg[3] = step;
  arg[4] = 0.0;
  arg[5] = 0.0;
  args.push_back( arg );
  arg[3] = 0.0;
  arg[4] = step;
  arg[5] = 0.0;
  args.push_back( arg );
  arg[3] = 0.0;
  arg[4] = 0.0;
  arg[5] = step;
  args.push_back( arg );
  // simple refinement
  Optimiser_simplex os( 0.0001, 50 );
  return local_rtop( os( *this, args ) );
}


// target function for NCS operator based on sphere of density

Target_fn_xmap_mask_rtop::Target_fn_xmap_mask_rtop( const clipper::Xmap<float>& source, const clipper::Xmap<float>& target, const clipper::NXmap<float>& mask, const double& rot_step, const double& trn_step )
{
  source_ = &source;
  target_ = &target;
  mask_ = &mask;
  rot_step_ = rot_step;
  trn_step_ = trn_step;
  // precalc some mask stats
  typedef clipper::NXmap<float>::Map_reference_index MRI;
  mmin = mmax = mask[mask.first()];
  for ( MRI inx = mask.first(); !inx.last(); inx.next() ) {
    if ( mask[inx] < mmin ) mmin = mask[inx];
    if ( mask[inx] > mmax ) mmax = mask[inx];
  }
}


double Target_fn_xmap_mask_rtop::operator() ( const Local_rtop& rot ) const
{
  typedef clipper::NXmap<float>::Map_reference_index MRI;
  typedef float                 TYPE;
  typedef clipper::Xmap<TYPE>   XMAP;
  typedef clipper::Interp_cubic INTERP;
  const clipper::Xmap<float>& source( *source_ );
  const clipper::Xmap<float>& target( *target_ );
  const clipper::NXmap<float>& mask( *mask_ );

  // construct operators
  clipper::RTop_orth rtop0( clipper::RTop<>::identity() );
  clipper::RTop_orth rtop1( rot.rtop_orth() );
  // and NX_operators, which provide optimised access
  clipper::NX_operator nxop0( source, mask, rtop0 );
  clipper::NX_operator nxop1( target, mask, rtop1 );
  // calculate map agreement
  double s0(0.0), s1(0.0);
  double rho0, rho1;
  double mcut = 0.1*mmin+0.9*mmax;
  for ( MRI inx = mask.first(); !inx.last(); inx.next() )
    if ( inx.coord().u()%2==0&&inx.coord().v()%2==0&&inx.coord().w()%2==0 )
     if ( mask[inx] < mcut ) {
      rho0 = nxop0.xmap_data<INTERP, TYPE, XMAP>( source, inx.coord() );
      rho1 = nxop1.xmap_data<INTERP, TYPE, XMAP>( target, inx.coord() );
      s0 += 1.0;
      s1 += clipper::Util::sqr((rho1-rho0)/mask[inx]);
    }
  s1 /= s0;
  return s1;
}


double Target_fn_xmap_mask_rtop::operator() ( const std::vector<double>& args ) const
{
  return (*this)( local_rtop( args ) );
}


Local_rtop Target_fn_xmap_mask_rtop::local_rtop( const std::vector<double>& args ) const
{
  return Local_rtop( clipper::Rotation(clipper::Euler<clipper::Rotation::EulerXYZs>(args[0],args[1],args[2])) * rot_.rot(), rot_.src(), clipper::Coord_orth(args[3],args[4],args[5]) + rot_.tgt() );
}


Local_rtop Target_fn_xmap_mask_rtop::refine( const Local_rtop& rot )
{
  // store initial rotation
  rot_ = rot;

  // calculate initial params
  std::vector<std::vector<double> > args;
  std::vector<double> arg(6,0.0);
  // identity
  clipper::Euler<clipper::Rotation::EulerXYZs> euler( 0.0, 0.0, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );
  // rotation steps
  double step = 0.5 * rot_step_;
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( step, 0.0, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( 0.0, step, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( 0.0, 0.0, step );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );
  // translation steps
  step = 0.5 * trn_step_;
  arg = args[0];
  arg[3] = step;
  arg[4] = 0.0;
  arg[5] = 0.0;
  args.push_back( arg );
  arg[3] = 0.0;
  arg[4] = step;
  arg[5] = 0.0;
  args.push_back( arg );
  arg[3] = 0.0;
  arg[4] = 0.0;
  arg[5] = step;
  args.push_back( arg );
  // simple refinement
  Optimiser_simplex os( 0.00001, 25 );
  return local_rtop( os( *this, args ) );
}


// Automatic NCS determination

std::vector<Local_rtop> Search_NCS_from_atom_map::operator() ( const clipper::Atom_list& atoms, const clipper::Xmap<float>& xmap ) const {
  // collect sites
  std::vector<clipper::Coord_orth> coords;
  for ( int i = 0; i < atoms.size(); i++ )
    coords.push_back( atoms[i].coord_orth() );

  // useful info
  clipper::Cell       cell = xmap.cell();
  clipper::Spacegroup spgr = xmap.spacegroup();

  // results
  std::vector<double> results_atmscr;
  std::vector<Local_rtop> results_uniqop;

  // find and evaluate density rotations candidate
  for ( int c1 = 0; c1 < coords.size()-1; c1++ )
    for ( int c2 = c1+1 ; c2 < coords.size(); c2++ ) {
      // calculate density rotations
      Gaussian_orth_list src, tgt, srccut, tgtcut;
      clipper::Coord_orth co, co1, co2;
      clipper::Coord_frac cf;
      // generate list of source and target atoms
      std::vector<std::pair<double,int> > src_index, tgt_index;
      for ( int i = 0; i < coords.size(); i++ ) if ( i != c1 && i != c2 )
	for ( int sym = 0; sym < spgr.num_symops(); sym++ ) {
	  co = spgr.symop(sym).rtop_orth(cell) * coords[i];
	  co1 = co - coords[c1];
	  co2 = co - coords[c2];
	  //cf  = co1.coord_frac(cell).lattice_copy_zero();
	  cf=co1.coord_frac(cell);
	  cf=cf-clipper::Coord_frac(rint(cf.u()),rint(cf.v()),rint(cf.w()));
	  src.push_back( Gaussian_orth( cf.coord_orth(cell), 3.0, 1.0 ) );
	  //cf  = co2.coord_frac(cell).lattice_copy_zero();
	  cf=co2.coord_frac(cell);
	  cf=cf-clipper::Coord_frac(rint(cf.u()),rint(cf.v()),rint(cf.w()));
	  tgt.push_back( Gaussian_orth( cf.coord_orth(cell), 3.0, 1.0 ) );
	}
      // sort by distance from centre
      for ( int i = 0; i < src.size(); i++ ) {
	src_index.push_back( std::pair<double,int>( src[i].coord_orth().lengthsq(), i ) );
	tgt_index.push_back( std::pair<double,int>( tgt[i].coord_orth().lengthsq(), i ) );
      }
      std::sort( src_index.begin(), src_index.end() );
      std::sort( tgt_index.begin(), tgt_index.end() );
      // truncate list
      int ncut = clipper::Util::min( int(src.size()), 25 );
      for ( int i = 0; i < ncut; i++ ) {
	srccut.push_back( src[ src_index[i].second ] );
	tgtcut.push_back( tgt[ tgt_index[i].second ] );
      }

      // create list of possible rotations matching these atoms
      Search_coord_rotation mr( 3.14159/18.0, 12 );
      std::vector<std::pair<double,clipper::Rotation> > rots = mr( srccut, tgtcut );

      // convert rotations to local_rtops
      std::vector<std::pair<double,Local_rtop> > result;
      for ( int r = 0; r < rots.size(); r++ )
	result.push_back( std::pair<double,Local_rtop>( rots[r].first,
	  Local_rtop( rots[r].second, coords[c1], coords[c2] ) ) );

      int nrot = clipper::Util::min( int(result.size()), 5 );

      // examine all candidates for this density rotation
      for ( int rot = 0; rot < nrot; rot++ ) {

	// calculate atom relationship score
	std::pair<int,double> atmscr = result[rot].second.atom_match( spgr, cell, coords, tol_dst_ );

	// check if this matches any previous operator
	std::pair<double,Local_rtop> rcbest( 1.0e9, Local_rtop() );
	int prv_match = 0;
	for ( int prv = 0; prv < results_uniqop.size(); prv++ ) {
	  std::pair<double,Local_rtop> rcmtch = result[rot].second.symm_match( results_uniqop[prv], spgr, cell, tol_dst_, tol_ang_ );
	  if ( rcmtch.first < rcbest.first ) {
	    rcbest = rcmtch;
	    prv_match = prv;
	  }
	}

	// Deal with new and duplicate matches:
	if ( rcbest.first < 2.0 ) {

	  // this is a duplicate - update previous if needed
	  if ( atmscr.second > results_atmscr[prv_match] ) {
	    results_atmscr[prv_match] = atmscr.second;
	    results_uniqop[prv_match] = rcbest.second;
	  }

	} else {

	  // reduce density rotations to proper-like form
	  Local_rtop rcprop = result[rot].second.proper( spgr, cell );

	  // calculate loop factor
	  std::pair<int,int> atmloo = rcprop.atom_loop( spgr, cell, coords, tol_dst_ );

	  // now look for proper operators
	  int nncs = atmloo.first;
	  int nmis = atmloo.first - atmloo.second;

	  results_atmscr.push_back( atmscr.second );
	  results_uniqop.push_back( rcprop );

	}
      }
    }

  // now do some comparisons on the density matching those ops
  Target_fn_xmap_rtop mr( xmap, xmap, 6.0, 0.1, 1.0, Target_fn_xmap_rtop::CORREL );
  // replace atom scores with density scores
  double s0(0.0), s1(0.0), s2(0.0);
  std::vector<double> results_denscr( results_uniqop.size() );
  for ( int i = 0; i < results_uniqop.size(); i++ ) {
    results_denscr[i] = mr( results_uniqop[i] );
    s0 += 1.0;
    s1 += results_denscr[i];
    s2 += results_denscr[i]*results_denscr[i];
  }
  s1 = s1/s0;
  s2 = sqrt( s2/s0 - s1*s1 );
  // create pruned list
  std::vector<Local_rtop> results_final;
  double cutoff = s1 - 3.0 * s2;  // 3 sigma cutoff
  for ( int i = 0; i < results_uniqop.size(); i++ )
    if ( results_denscr[i] < cutoff )
      results_final.push_back( /*mr.refine*/( results_uniqop[i] ) );

  // include inverses and return
  return include_inverse( results_final, spgr, cell );
}


std::vector<Local_rtop> Search_NCS_from_atom_map::include_inverse( const std::vector<Local_rtop>& ncsops, const clipper::Spacegroup& spgr, const clipper::Cell& cell ) const
{
  std::vector<Local_rtop> ncsopsi;
  for ( int i = 0; i < ncsops.size(); i++ ) {
    Local_rtop rtinv = ncsops[i].inverse();
    int j;
    for ( j = 0; j < ncsops.size(); j++ )
      if ( rtinv.symm_match( ncsops[j], spgr, cell, tol_dst_, tol_ang_ ).first
	   < 2.0 ) break;
    ncsopsi.push_back( ncsops[i] );
    if ( j == ncsops.size() ) ncsopsi.push_back( rtinv );
  }
  return ncsopsi;
}


std::vector<Local_rtop> Search_NCS_from_atom_map::exclude_inverse( const std::vector<Local_rtop>& ncsops, const clipper::Spacegroup& spgr, const clipper::Cell& cell ) const
{
  std::vector<Local_rtop> ncsopsi;
  for ( int i = 0; i < ncsops.size(); i++ ) {
    Local_rtop rtinv = ncsops[i].inverse();
    int j;
    for ( j = 0; j < ncsopsi.size(); j++ )
      if ( rtinv.symm_match( ncsopsi[j], spgr, cell, tol_dst_, tol_ang_ ).first
	   < 2.0 ) break;
    if ( j == ncsopsi.size() ) ncsopsi.push_back( ncsops[i] );
  }
  return ncsopsi;
}


Local_rtop Search_NCS_from_atom_map::masked_refine_from_map( const clipper::Xmap<float>& xmap, const Local_rtop& nxop, const double& map_radius, const double& local_radius )
{
  typedef clipper::NXmap<float>::Map_reference_index MRI;
  typedef float                 TYPE;
  typedef clipper::Xmap<TYPE>   XMAP;
  typedef clipper::Interp_cubic INTERP;
  clipper::Coord_orth orth0, orth1;
  clipper::Coord_grid grid0;

  // fetch operator components
  orth0 = nxop.src();
  // get grid offset for source coordinate
  grid0 = orth0.coord_frac( xmap.cell() ).coord_grid( xmap.grid_sampling() );
  // calc grid containing the desired volume
  clipper::Grid_range gr0( xmap.cell(), xmap.grid_sampling(), map_radius );
  // and offset by the base coordinate
  clipper::Grid_range gr1( gr0.min() + grid0, gr0.max() + grid0 );

  // construct new operators
  clipper::RTop_orth rtop0( clipper::RTop<>::identity() );
  clipper::RTop_orth rtop1( nxop.rtop_orth() );

  // make 2 nxmaps, one each for unrotated and rotated density
  clipper::NXmap<float> nxmap0( xmap.cell(), xmap.grid_sampling(), gr1 );
  clipper::NXmap<float> nxmap1( xmap.cell(), xmap.grid_sampling(), gr1 );
  // and NX_operators, which provide optimised access
  clipper::NX_operator nxop0( xmap, nxmap0, rtop0 );
  clipper::NX_operator nxop1( xmap, nxmap1, rtop1 );

  // populate the unrotated and rotated nxmap
  clipper::Xmap<float>::Map_reference_coord ix( xmap );
  clipper::Coord_frac cf;
  for ( MRI inx = nxmap0.first(); !inx.last(); inx.next() ) {
    nxmap0[inx] = nxop0.xmap_data<INTERP, TYPE, XMAP>( xmap, inx.coord() );
    nxmap1[inx] = nxop1.xmap_data<INTERP, TYPE, XMAP>( xmap, inx.coord() );
  }

  // now calculate the local correlation contributions
  // prepare local averages
  clipper::MapFilterFn_step step( local_radius );
  clipper::MapFilter_fft<float>
    filter( step, 1.0, clipper::MapFilter_fft<float>::Relative );
  // make the maps
  clipper::NXmap<float> x0( xmap.cell(), xmap.grid_sampling(), gr1 );
  clipper::NXmap<float> x1( xmap.cell(), xmap.grid_sampling(), gr1 );
  clipper::NXmap<float> x00( xmap.cell(), xmap.grid_sampling(), gr1 );
  clipper::NXmap<float> x11( xmap.cell(), xmap.grid_sampling(), gr1 );
  clipper::NXmap<float>& prd = x0;
  clipper::NXmap<float>& sub = x1;
  clipper::NXmap<float>& cov = x00;
  clipper::NXmap<float>& var = x11;
  clipper::NXmap<float>& correl = x0;
  clipper::NXmap<float>& rmsdev = x1;
  // do local averages
  // variance terms
  for ( MRI inx = nxmap0.first(); !inx.last(); inx.next() ) {
    x0[inx] = nxmap0[inx] * nxmap0[inx];
    x1[inx] = nxmap1[inx] * nxmap1[inx];
  }
  filter( x00, x0 );     // second moment
  filter( x11, x1 );     // second moment
  filter( x0, nxmap0 );  // first moment
  filter( x1, nxmap1 );  // first moment
  for ( MRI inx = nxmap0.first(); !inx.last(); inx.next() ) {
    var[inx] = sqrt( (x00[inx]-x0[inx]*x0[inx]) * (x11[inx]-x1[inx]*x1[inx]) );
    sub[inx] = x0[inx] * x1[inx];
    prd[inx] = nxmap0[inx] * nxmap1[inx];
  }
  filter( cov, prd );    // covariance
  for ( MRI inx = nxmap0.first(); !inx.last(); inx.next() ) {
    cov[inx]    = cov[inx] - sub[inx];
    correl[inx] = cov[inx] / var[inx];
    sub[inx]    = clipper::Util::sqr( nxmap0[inx] - nxmap1[inx] );
  }
  filter( rmsdev, sub ); // rms difference

  // calculate 'remote' correlation statistics
  double rad2_lo = 0.7 * map_radius * map_radius;
  double rad2_hi = 0.8 * map_radius * map_radius;
  double c0(0.0), c1(0.0), c2(0.0), cncs(0.0), dncs(0.0), dmap(0.0);
  for ( MRI inx = nxmap0.first(); !inx.last(); inx.next() ) {
    double rad2 = (inx.coord_orth()-orth0).lengthsq();
    if ( correl[inx] > cncs ) {
      cncs = correl[inx];                            // max correl
      dncs = rmsdev[inx];                            // rmsd at max correl
    }
    if ( rad2 > rad2_lo && rad2 < rad2_hi ) {
      c0 += 1.0;
      c1 += correl[inx];
      c2 += correl[inx]*correl[inx];
      if ( rmsdev[inx] > dmap ) dmap = rmsdev[inx];  // max rmsd
    }
  }
  c1 /= c0; c2 /= c0;
  c2 = sqrt( c2 - c1*c1 );  // mean and std dev of remote correl

  // calculate map variance and centre of mass
  clipper::Coord_orth cmass( 0.0, 0.0, 0.0 );
  double mass(0.0), weight(0.0);
  double pncs0 = 0.01;
  for ( MRI inx = nxmap0.first(); !inx.last(); inx.next() ) {
    double dc = clipper::Util::max( (correl[inx]-c1)/c2, 0.0 );
    double pncs = 1.0/(1.0+exp(-0.5*dc*dc)/pncs0);
    double pmap = 1.0 - pncs;
    double vmap = (dmap*(pmap*dmap+pncs*dncs)) / (pncs*(dmap-dncs));
    var[inx] = sqrt(vmap);             // map of standard deviations
    weight = 1.0/vmap;
    mass  += weight;
    cmass += weight*inx.coord_orth();  // inverse variance weighted centre
  }
  cmass = (1.0/mass) * cmass;

  // refine operator using mask
  Target_fn_xmap_mask_rtop mr( xmap, xmap, var, 0.01, 0.1 );
  Local_rtop refop = mr.refine( nxop );
  Local_rtop newop( refop.rot(), cmass, refop.rtop_orth()*cmass );

  return newop;
}




void Xmap_ncs::restrain_ncs( const clipper::Xmap<float>& xmap, const Local_rtop& nxop, const double& map_radius, const double& local_radius )
{
  typedef clipper::NXmap<float>::Map_reference_index MRI;
  typedef float                 TYPE;
  typedef clipper::Xmap<TYPE>   XMAP;
  typedef clipper::Interp_cubic INTERP;
  clipper::Coord_orth orth0, orth1;
  clipper::Coord_grid grid0;

  // fetch operator components
  orth0 = nxop.src();
  // get grid offset for source coordinate
  grid0 = orth0.coord_frac( xmap.cell() ).coord_grid( xmap.grid_sampling() );
  // reset source coordinate
  orth0 = grid0.coord_frac( xmap.grid_sampling() ).coord_orth( xmap.cell() );
  // and target coordinate
  orth1 = nxop.rtop_orth() * orth0;
  // construct new operators
  clipper::RTop_orth rtop0( clipper::Mat33<>::identity(), orth0 );
  clipper::RTop_orth rtop1( nxop.rtop_orth().rot(), orth1 );

  // make 2 nxmaps, one each for unrotated and rotated density
  clipper::Grid_range grid_range( xmap.cell(), xmap.grid_sampling(), map_radius );
  clipper::NXmap<float> nxmap0( xmap.cell(), xmap.grid_sampling(), grid_range );
  clipper::NXmap<float> nxmap1( xmap.cell(), xmap.grid_sampling(), grid_range );
  // and NX_operators, which provide optimised access
  clipper::NX_operator nxop0( xmap, nxmap0, rtop0 );
  clipper::NX_operator nxop1( xmap, nxmap1, rtop1 );

  // populate the unrotated and rotated nxmap
  clipper::Xmap<float>::Map_reference_coord ix( xmap );
  clipper::Coord_frac cf;
  for ( MRI inx = nxmap0.first(); !inx.last(); inx.next() ) {
    nxmap0[inx] = nxop0.xmap_data<INTERP, TYPE, XMAP>( xmap, inx.coord() );
    nxmap1[inx] = nxop1.xmap_data<INTERP, TYPE, XMAP>( xmap, inx.coord() );
  }

  // now calculate the local correlation contributions
  // prepare local averages
  clipper::MapFilterFn_step step( local_radius );
  clipper::MapFilter_fft<float>
    filter( step, 1.0, clipper::MapFilter_fft<float>::Relative );
  // make the maps
  clipper::NXmap<float> x0( xmap.cell(), xmap.grid_sampling(), grid_range );
  clipper::NXmap<float> x1( xmap.cell(), xmap.grid_sampling(), grid_range );
  clipper::NXmap<float> x00( xmap.cell(), xmap.grid_sampling(), grid_range );
  clipper::NXmap<float> x11( xmap.cell(), xmap.grid_sampling(), grid_range );
  clipper::NXmap<float>& prd = x0;
  clipper::NXmap<float>& sub = x1;
  clipper::NXmap<float>& cov = x00;
  clipper::NXmap<float>& var = x11;
  clipper::NXmap<float>& correl = x0;
  clipper::NXmap<float>& rmsdev = x1;
  // do local averages
  // variance terms
  for ( MRI inx = nxmap0.first(); !inx.last(); inx.next() ) {
    x0[inx] = nxmap0[inx] * nxmap0[inx];
    x1[inx] = nxmap1[inx] * nxmap1[inx];
  }
  filter( x00, x0 );     // second moment
  filter( x11, x1 );     // second moment
  filter( x0, nxmap0 );  // first moment
  filter( x1, nxmap1 );  // first moment
  for ( MRI inx = nxmap0.first(); !inx.last(); inx.next() ) {
    var[inx] = sqrt( (x00[inx]-x0[inx]*x0[inx]) * (x11[inx]-x1[inx]*x1[inx]) );
    sub[inx] = x0[inx] * x1[inx];
    prd[inx] = nxmap0[inx] * nxmap1[inx];
  }
  filter( cov, prd );    // covariance
  for ( MRI inx = nxmap0.first(); !inx.last(); inx.next() ) {
    cov[inx]    = cov[inx] - sub[inx];
    correl[inx] = cov[inx] / var[inx];
    sub[inx]    = clipper::Util::sqr( nxmap0[inx] - nxmap1[inx] );
  }
  filter( rmsdev, sub ); // rms difference

  // calculate 'remote' correlation statistics
  double rad2_lo = 0.7 * map_radius * map_radius;
  double rad2_hi = 0.8 * map_radius * map_radius;
  double c0(0.0), c1(0.0), c2(0.0), cncs(0.0), dncs(0.0), dmap(0.0);
  for ( MRI inx = nxmap0.first(); !inx.last(); inx.next() ) {
    double rad2 = inx.coord_orth().lengthsq();
    if ( correl[inx] > cncs ) {
      cncs = correl[inx];                            // max correl
      dncs = rmsdev[inx];                            // rmsd at max correl
    }
    if ( rad2 > rad2_lo && rad2 < rad2_hi ) {
      c0 += 1.0;
      c1 += correl[inx];
      c2 += correl[inx]*correl[inx];
      if ( rmsdev[inx] > dmap ) dmap = rmsdev[inx];  // max rmsd
    }
  }
  c1 /= c0; c2 /= c0;
  c2 = sqrt( c2 - c1*c1 );  // mean and std dev of remote correl

  // calculate map variance
  double pncs0 = 0.01;
  for ( MRI inx = nxmap0.first(); !inx.last(); inx.next() ) {
    double dc = clipper::Util::max( (correl[inx]-c1)/c2, 0.0 );
    double pncs = 1.0/(1.0+exp(-0.5*dc*dc)/pncs0);
    double pmap = 1.0 - pncs;
    double vmap = (dmap*(pmap*dmap+pncs*dncs)) / (pncs*(dmap-dncs));
    clipper::Coord_grid xcoord = inx.coord() + grid_range.min() + grid0;
    Gaussian_probability_1d g1 = get_data(xcoord);
    Gaussian_probability_1d g2( nxmap1[inx], sqrt(vmap) );
    set_data( xcoord, g1*g2 );
  }
}
