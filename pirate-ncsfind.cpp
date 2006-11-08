/*! \file pirate-ncsfind.cpp pirate library */
/* Copyright 2003-2006 Kevin Cowtan & University of York all rights reserved */

#include "pirate-ncsfind.h"


struct Compare_pair_score_rotation{ bool operator() ( const std::pair<double,clipper::Rotation>& r1, const std::pair<double,clipper::Rotation>& r2 ) const { return ( r1.first < r2.first ); } };


/*
 Transform an operator by pre- and post- RToperators.
*/
Local_rtop Local_rtop::transform( const clipper::RTop_orth& r1, const clipper::RTop_orth& r2 ) const
{
  Local_rtop res;
  res.src() = r1 * src();
  res.tgt() = r2 * tgt();
  res.rot() = clipper::Rotation(r2.rot() * rot().matrix() * r1.rot().inverse());
  return res;
}


/*
 Calculate whether two operators match to within a given tolerance.
*/
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
      double a = 2.0*acos(clipper::Util::min(fabs(rot.w()),1.0));
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


/*
 Calculate closest operator to proper form.
*/
Local_rtop Local_rtop::proper( const clipper::Spacegroup& spgr, const clipper::Cell& cell ) const
{
  // prepare orthogonal symops
  std::vector<clipper::RTop_orth> symop_orth( spgr.num_symops() );
  for ( int i = 0; i < symop_orth.size(); i++ )
    symop_orth[i] = spgr.symop(i).rtop_orth(cell);

  // most proper version
  int maxncs = 20;
  double wloop( 1.0 ), wscrw( 1.0 ), wtran( 0.001 );
  // make list of candidate operators
  std::vector<std::pair<double,Local_rtop> > resultsym;
  clipper::RTop_orth rtid = clipper::RTop_orth::identity();
  for ( int j = 0; j < spgr.num_symops(); j++ ) {
    Local_rtop rsym, rcel;
    rsym = transform( rtid, symop_orth[j] );
    clipper::Coord_frac cf = (rsym.tgt()-rsym.src()).coord_frac(cell);
    clipper::Coord_frac df( floor(cf.u()), floor(cf.v()), floor(cf.w()) );
    clipper::Coord_frac d;
    for ( int u = 0; u <= 1; u++ )
      for ( int v = 0; v <= 1; v++ )
	for ( int w = 0; w <= 1; w++ ) {
	  d = df + clipper::Coord_frac(double(u),double(v),double(w));
	  rcel.src() = rsym.src();
	  rcel.tgt() = rsym.tgt() - d.coord_orth(cell);
	  rcel.rot() = rsym.rot();
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
    clipper::RTop_orth rtop = resultsym[k].second.rtop_orth();
    double scrwshift = clipper::Util::sqr( v0 * rtop.trn() );
    double transhift = ( resultsym[k].second.tgt() - resultsym[k].second.src() ).lengthsq();
    double loopshift = 1.0e9;
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


/*
 Calculate a proper closed group based on this one operator, if possible.
 (Operator must already be in proper form)
*/
std::vector<Local_rtop> Local_rtop::proper_closed_group( const clipper::Spacegroup& spgr, const clipper::Cell& cell, const std::vector<int>& orders, double tol_dst ) const
{
  std::vector<Local_rtop> result;

  // first find the best loop factor
  clipper::RTop_orth rtop = rtop_orth();
  int maxncs = 20;
  int rmin = 0;
  double d, dmin( tol_dst*tol_dst );
  clipper::Coord_orth init, next;
  init = src();
  next = rtop * init;
  for ( int r = 1; r < maxncs; r++ ) {
    d = (next-init).lengthsq();
    if ( d < dmin ) {
      dmin = d;
      rmin = r;
    }
    next = rtop * next;
  }

  // check if that factor is allowed
  bool allowed = false;
  for ( int i = 0; i < orders.size(); i++ )
    if ( rmin == orders[i] ) allowed = true;

  // if so, form the exact operator
  if ( allowed ) {
    // find centre of mass of group
    clipper::Coord_orth cent( 0.0, 0.0, 0.0 );
    next = src();
    for ( int r = 0; r < rmin; r++ ) {
      cent += next;
      next = rtop * next;
    }
    cent = (1.0/double(rmin)) * cent;
    // and exact rotation
    clipper::Polar_ccp4 p = rot().polar_ccp4();
    for ( int r = 1; r < rmin; r++ )
      result.push_back( Local_rtop( clipper::Rotation( clipper::Polar_ccp4( p.psi(), p.phi(), clipper::Util::twopi()*double(r)/double(rmin) ) ), cent, cent ) );
  }

  // and return the result
  return result;
}


/*
 Ensure all distinct inverse operators are present
*/
std::vector<Local_rtop> Local_rtop::include_inverse( const std::vector<Local_rtop>& ncsops, const clipper::Spacegroup& spgr, const clipper::Cell& cell, double tol_dst, double tol_ang )
{
  std::vector<Local_rtop> ncsopsi;
  for ( int i = 0; i < ncsops.size(); i++ ) {
    Local_rtop rtinv = ncsops[i].inverse();
    int j;
    for ( j = 0; j < ncsops.size(); j++ )
      if ( rtinv.symm_match( ncsops[j], spgr, cell, tol_dst, tol_ang ).first
	   < 2.0 ) break;
    ncsopsi.push_back( ncsops[i] );
    if ( j == ncsops.size() ) ncsopsi.push_back( rtinv );
  }
  return ncsopsi;
}


/*
 Ensure all distinct inverse operators are absent
*/
std::vector<Local_rtop> Local_rtop::exclude_inverse( const std::vector<Local_rtop>& ncsops, const clipper::Spacegroup& spgr, const clipper::Cell& cell, double tol_dst, double tol_ang )
{
  std::vector<Local_rtop> ncsopsi;
  for ( int i = 0; i < ncsops.size(); i++ ) {
    Local_rtop rtinv = ncsops[i].inverse();
    int j;
    for ( j = 0; j < ncsopsi.size(); j++ )
      if ( rtinv.symm_match( ncsopsi[j], spgr, cell, tol_dst, tol_ang ).first
	   < 2.0 ) break;
    if ( j == ncsopsi.size() ) ncsopsi.push_back( ncsops[i] );
  }
  return ncsopsi;
}


/*
 Remove duplicates and identity operators
*/
std::vector<Local_rtop> Local_rtop::tidy( const std::vector<Local_rtop>& ncsops, const clipper::Spacegroup& spgr, const clipper::Cell& cell, double tol_dst, double tol_ang )
{
  std::vector<Local_rtop> ncsopsi;
  Local_rtop rtid( clipper::Rotation::zero(),
		   clipper::Coord_orth(0.0,0.0,0.0),
		   clipper::Coord_orth(0.0,0.0,0.0) );
  for ( int i = 0; i < ncsops.size(); i++ ) {
    Local_rtop rtop = ncsops[i];
    bool keep = true;
    for ( int j = 0; j < ncsopsi.size(); j++ )
      if ( rtop.symm_match( ncsopsi[j], spgr, cell, tol_dst, tol_ang ).first
	   < 2.0 ) { keep = false; break; }
    if ( rtop.symm_match( rtid, spgr, cell, tol_dst, tol_ang ).first
	 < 2.0 ) { keep = false; }
    if ( keep ) ncsopsi.push_back( ncsops[i] );
  }
  return ncsopsi;
}


// Target function for superimposing coordinates by rotation

Target_fn_coord_rotation::Target_fn_coord_rotation( const Atom_overlap_list& source, const Atom_overlap_list& target, const double& rot_step )
{
  source_ = &source;
  target_ = &target;
  rot_step_ = rot_step;
}

double Target_fn_coord_rotation::operator() ( const clipper::Rotation& rot ) const
{
  const Atom_overlap_list& source = *source_;
  const Atom_overlap_list& target = *target_;
  clipper::Mat33<> mat = rot.matrix();
  double r2, result = 0.0;
  for ( int s = 0; s < source.size(); s++ ) {
    clipper::Coord_orth source_coord( mat * source[s].coord_orth() );
    for ( int t = 0; t < target.size(); t++ ) {
      r2 = ( source_coord - target[t].coord_orth() ).lengthsq();
      result += source[s].evaluate(r2) * target[t].evaluate(r2);
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
  Optimiser_simplex os( 0.0001, 25 );
  clipper::Rotation op = rotation( os( *this, args ) );
  return op.norm();
}


Search_coord_rotation::Search_coord_rotation( const double& rot_step, const int& max_match )
{
  rot_step_ = rot_step;
  max_match_ = max_match;
}

std::vector<std::pair<double,clipper::Rotation> > Search_coord_rotation::operator() ( const Atom_overlap_list& source, const Atom_overlap_list& target ) const
{
  std::vector<std::pair<double,clipper::Rotation> > result, result_tmp;
  clipper::Rotation r0, r1;

  // limit maximum number of atoms
  int nmatch  = clipper::Util::min( source.size(), target.size() );
  nmatch = clipper::Util::min( nmatch, max_match_ );

  // first make a list of rotation candidates
  std::vector<clipper::Coord_orth> v1(3), v2(3);
  // loop over first source and target atom
  for ( int s1 = 1; s1 < nmatch-1; s1++ ) {
    double ls1 = sqrt(source[s1].coord_orth().lengthsq());
    double ws1 = source[s1].half_width();
    for ( int t1 = 1; t1 < nmatch; t1++ ) {
      double lt1 = sqrt(target[t1].coord_orth().lengthsq());
      double wt1 = target[t1].half_width();
      if ( fabs(lt1-ls1) < 0.50*(ws1+wt1) ) {
	// loop over second source and target atom
	for ( int s2 = s1+1; s2 < nmatch; s2++ ) {
	  double ls2 = sqrt(source[s2].coord_orth().lengthsq());
	  double ws2 = source[s2].half_width();
	  for ( int t2 = 1; t2 < nmatch; t2++ ) {
	    double lt2 = sqrt(target[t2].coord_orth().lengthsq());
	    double wt2 = target[t2].half_width();
	    if ( fabs(lt2-ls2) < 0.50*(ws2+wt2) ) {
	      // check third edge of triangle
	      double ls12 = sqrt( (source[s2].coord_orth()-source[s1].coord_orth()).lengthsq() );
	      double lt12 = sqrt( (target[t2].coord_orth()-target[t1].coord_orth()).lengthsq() );
	      if ( fabs(ls12-lt12) < 0.25*(ws1+ws2+wt1+wt2) ) {
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
  for ( int r = 0; r < result.size(); r++ ) 
    result[r].first = tf( result[r].second );
  std::sort( result.begin(), result.end(), Compare_pair_score_rotation() );

  // remove duplicates
  for ( int r = 0; r < result.size(); r++ ) {
    r0 = result[r].second.inverse();
    bool keep = true;
    for ( int prv = 0; prv < result_tmp.size(); prv++ ) {
      r1 = r0 * result_tmp[prv].second;
      double a = 2.0*acos(clipper::Util::min(fabs(r1.w()),1.0));
      if ( a < 0.5 ) { keep = false; break; }       // must be 30 deg apart
    }
    if ( keep ) result_tmp.push_back( result[r] );  // only add unique ops
    if ( result_tmp.size() == 5 ) break;            // stop after 5
  }
  result = result_tmp;

  // refine and resort
  for ( int i = 0; i < result.size(); i++ ) {
      result[i].second = tf.refine( result[i].second );
      result[i].first  = tf( result[i].second );
  }
  std::sort( result.begin(), result.end(), Compare_pair_score_rotation() );

  // remove duplicates
  for ( int r = 0; r < result.size(); r++ ) {
    r0 = result[r].second.inverse();
    bool keep = true;
    for ( int prv = 0; prv < result_tmp.size(); prv++ ) {
      r1 = r0 * result_tmp[prv].second;
      double a = 2.0*acos(clipper::Util::min(fabs(r1.w()),1.0));
      if ( a < 2.0*rot_step_ ) { keep = false; break; }
    }
    if ( keep ) result_tmp.push_back( result[r] );  // only add unique ops
  }
  result = result_tmp;

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
  Local_rtop op = local_rtop( os( *this, args ) );
  return Local_rtop( op.rot().norm(), op.src(), op.tgt() );
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
  Optimiser_simplex os( 0.00001, 50 );
  Local_rtop op = local_rtop( os( *this, args ) );
  return Local_rtop( op.rot().norm(), op.src(), op.tgt() );
}


// Automatic NCS determination

/*
 Find NCS operators from atoms and density.
*/
std::vector<Local_rtop> NCSfind::operator() ( const clipper::Atom_list& atoms, const clipper::Xmap<float>& xmap ) const {
  // useful info
  clipper::Spacegroup spgr = xmap.spacegroup();
  clipper::Cell       cell = xmap.cell();

  // get results from atom matching
  std::vector<Local_rtop> results_atom =
    find_ncs_candidates( atoms, spgr, cell );

  // filter results from xmap correlation
  std::vector<Local_rtop> results_xmap =
    filter_ncs_candidates( results_atom, atoms, xmap, 2.0 );

  const int orders_init[] = {2,3,4,6};
  std::vector<int> orders( orders_init, orders_init+4 );
  std::vector<int> order2( orders_init, orders_init+1 );
  // now get the best proper operator
  std::vector<Local_rtop> locked_best, locked_temp;
  for ( int i = 0; i < results_xmap.size(); i++ ) {
    locked_temp =
      results_xmap[i].proper_closed_group( spgr, cell, orders, 1.5*tol_dst_ );
    if ( locked_temp.size() > locked_best.size() )
      locked_best = locked_temp;
  }

  // did we find a proper operator?
  if ( locked_best.size() > 0 ) {
    // remove the operator from the candidate list
    std::vector<Local_rtop> results_temp;
    for ( int i = 0; i < results_xmap.size(); i++ ) {
      double match = false;
      for ( int j = 0; j < locked_best.size(); j++ )
	if ( results_xmap[i].symm_match( locked_best[j], spgr, cell, tol_dst_, tol_ang_ ).first < 8.0 ) match = true;
      if ( !match ) results_temp.push_back( results_xmap[i] );
    }

    // if there is a proper operator, store it and look for a second operator
    results_xmap = locked_best;

    // look for a perpendicular operator
    Local_rtop rts, rtb, rto( results_xmap[0] );
    for ( int i = 0; i < results_temp.size(); i++ ) {
      // check for pure 2-fold
      locked_temp =
	results_temp[i].proper_closed_group( spgr, cell, order2, 1.0*tol_dst_ );
      if ( locked_temp.size() == 1 ) {
	// find best symmetry copy
	double dmin = 1.0e9;
	for ( int s = 0; s < spgr.num_symops(); s++ ) {
	  clipper::RTop_orth so = spgr.symop(s).rtop_orth( cell );
	  rts = locked_temp[0].transform( so, so );
	  // check for orthogonality
	  clipper::Coord_orth v1( rto.rot().x(), rto.rot().y(), rto.rot().z() );
	  clipper::Coord_orth v2( rts.rot().x(), rts.rot().y(), rts.rot().z() );
	  if ( v1.unit() * v2.unit() < tol_ang_ ) {
	    // find closest lattice equivalent
	    //rts.src() = rts.tgt() = rts.src().coord_frac(cell).lattice_copy_near(rto.src().coord_frac(cell)).coord_orth(cell);
	    clipper::Coord_frac cf = ( rts.src() - rto.src() ).coord_frac(cell);
	    cf = clipper::Coord_frac(rint(cf.u()),rint(cf.v()),rint(cf.w()));
	    rts.src() = rts.tgt() = rts.src() - cf.coord_orth(cell);
	    double d = minimum_distance_between_lines( rto.src(), v1,
						       rts.src(), v2 );
	    if ( d < dmin ) {
	      dmin = d;
	      rtb = rts;
	    }
	  }
	}
	// if the centre of the symm is nearby, then keep
	if ( dmin < 0.5*tol_dst_ ) results_xmap.push_back( rtb );
	if ( debug_ ) std::cout << "Perpendicular op: " << dmin << "\n";
      }
    }
  } else {
    // if no proper operators, try for improper ones
    results_xmap = filter_ncs_candidates( results_atom, atoms, xmap, 4.0 );
    results_xmap = Local_rtop::include_inverse( results_xmap, spgr, cell );
  }

  // include inverses and return
  return Local_rtop::tidy( results_xmap, spgr, cell );
}


/*
 Find NCS candidate operators from atoms
*/
std::vector<Local_rtop> NCSfind::find_ncs_candidates( const clipper::Atom_list& atoms, const clipper::Spacegroup& spgr, const clipper::Cell& cell ) const {
  // collect sites
  std::vector<clipper::Coord_orth> coords;
  for ( int i = 0; i < atoms.size(); i++ )
    coords.push_back( atoms[i].coord_orth() );

  // results
  std::vector<Local_rtop> results_final;
  std::vector<std::pair<double,Local_rtop> > results_all;
  std::vector<std::pair<double,Local_rtop> > results_unq;

  // find and evaluate density rotations candidate
  for ( int c1 = 0; c1 < coords.size()-1; c1++ )
    for ( int c2 = c1+1 ; c2 < coords.size(); c2++ ) {
      // calculate density rotations
      Atom_overlap_list src, tgt, srccut, tgtcut;
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
	  src.push_back( Atom_overlap( cf.coord_orth(cell), 4.0, 1.0 ) );
	  //cf  = co2.coord_frac(cell).lattice_copy_zero();
	  cf=co2.coord_frac(cell);
	  cf=cf-clipper::Coord_frac(rint(cf.u()),rint(cf.v()),rint(cf.w()));
	  tgt.push_back( Atom_overlap( cf.coord_orth(cell), 4.0, 1.0 ) );
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
      Search_coord_rotation mr( 2.0*tol_ang_, 12 );
      std::vector<std::pair<double,clipper::Rotation> > rots = mr( srccut, tgtcut );

      // convert rotations to local_rtops
      int nrot = clipper::Util::min( int(rots.size()), 5 );
      for ( int r = 0; r < nrot; r++ )
	results_all.push_back( std::pair<double,Local_rtop>( rots[r].first,
	  Local_rtop( rots[r].second, coords[c1], coords[c2] ) ) );
    }

  // check for empty list
  if ( results_all.size() == 0 ) return results_final;

  // Calculate cutoff
  std::sort( results_all.begin(), results_all.end() );
  double cutoff = 0.5 * ( results_all[0].first + 3.0 ) - 3.0;

  // count the useful points
  int nrot;
  for ( nrot = 0; nrot < results_all.size(); nrot++ )
    if ( results_all[nrot].first > cutoff ) break;
  if ( nrot > 12*atoms.size() ) nrot = 12*atoms.size();  // <24 matches/atom

  // debug info
  if ( debug_ ) {
    std::cout << "Atom list filter: " << results_all.size() << "\t" << nrot << "\t" << cutoff << std::endl;
    for ( int i = 0; i < clipper::Util::min(2*nrot,int(results_all.size())); i++ ) std::cout << i << "\t" << results_all[i].first << "\t" << results_all[i].second.rot().polar_ccp4().format() << std::endl;
  }  // debug info

  // return the result
  for ( int i = 0; i < nrot; i++ )
    results_final.push_back( results_all[i].second );

  return results_final;
}


/*
 Filter NCS candidate operators against map
*/
std::vector<Local_rtop> NCSfind::filter_ncs_candidates( const std::vector<Local_rtop>& rtops, const clipper::Atom_list& atoms, const clipper::Xmap<float>& xmap, double sigma ) const
{
  // useful info
  clipper::Spacegroup spgr = xmap.spacegroup();
  clipper::Cell       cell = xmap.cell();
  // results
  std::vector<Local_rtop> results_final;
  if ( rtops.size() == 0 || atoms.size() < 2 ) return results_final;

  // calculate density scores
  Target_fn_xmap_rtop mr( xmap, xmap, 6.0, 0.1, 1.0, Target_fn_xmap_rtop::CORREL );
  std::vector<std::pair<double,Local_rtop> > results_all(rtops.size());
  for ( int i = 0; i < rtops.size(); i++ ) {
    results_all[i].first  = mr( rtops[i] );
    results_all[i].second = rtops[i];
  }
  std::sort( results_all.begin(), results_all.end() );

  // calculate density cutoff from arbitrary RTops
  double s0(0.0), s1(0.0), s2(0.0), cutoff;
  for ( int r = 0; r < 5; r++ )
    for ( int a0 = 0; a0 < atoms.size()-1; a0++ )
      for ( int a1 = a0+1; a1 < atoms.size(); a1++ ) {
	if ( s0 > 1000.0 ) break;
	double r1 = double(r+a0+a1);
	clipper::Euler_ccp4 rot( 2.0*r1, 3.0*r1, 5.0*r1 );
        clipper::Rotation rtmp( rot );
        Local_rtop rtrand( rtmp,
                          atoms[a0].coord_orth(), atoms[a1].coord_orth() );
	double result = mr( rtrand );
	s0 += 1.0;
	s1 += result;
	s2 += result*result;
      }
  s1 = s1/s0;
  s2 = sqrt( s2/s0 - s1*s1 );
  cutoff = s1 - sigma * s2;  // n sigma cutoff

  // count the useful points
  int nrot;
  for ( nrot = 0; nrot < results_all.size(); nrot++ )
    if ( results_all[nrot].first > cutoff ) break;

  // debug info
  if ( debug_ ) {
    std::cout << "Density filter: " << results_all.size() << "\t" << nrot << "\t" << cutoff << std::endl;
    for ( int i = 0; i < results_all.size(); i++ )
      std::cout << i << "\t" << results_all[i].first << "\t" << results_all[i].second.src().format() << "\t" << results_all[i].second.tgt().format() << "\t" << results_all[i].second.rot().polar_ccp4().format() << std::endl;
  }  // debug info

  // and extract them
  for ( int i = 0; i < nrot; i++ )
    results_final.push_back( results_all[i].second );

  // remove identities, duplicates and inverses
  results_final = Local_rtop::tidy( results_final, spgr, cell,
				    tol_dst_, tol_ang_ );
  results_final = Local_rtop::exclude_inverse( results_final, spgr, cell,
					       tol_dst_, tol_ang_ );

  // refine
  //for ( int i = 0; i < nrot; i++ )
  //results_final[i] = mr.refine( results_final[i] );

  // properise
  for ( int i = 0; i < results_final.size(); i++ )
    results_final[i] = results_final[i].proper( spgr, cell );

  // tidy and return
  return results_final;
}


/*
 Calculate rotated densities, local rmsd and correl
*/
std::vector<double> NCSfind::local_variance_from_map( clipper::NXmap<float>& nxmap0, clipper::NXmap<float>& nxmap1, clipper::NXmap<float>& correl, clipper::NXmap<float>& rmsdev, clipper::NXmap<float>& mskdev, const clipper::Xmap<float>& xmap, const Local_rtop& nxop, const double& map_radius, const double& local_radius, const double& pncs0 )
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
  // calc grid containing the desired volume
  clipper::Grid_range gr0( xmap.cell(), xmap.grid_sampling(), map_radius );
  // and offset by the base coordinate
  clipper::Grid_range gr1( gr0.min() + grid0, gr0.max() + grid0 );

  // construct new operators
  clipper::RTop_orth rtop0( clipper::RTop<>::identity() );
  clipper::RTop_orth rtop1( nxop.rtop_orth() );

  // init 2 nxmaps, one each for unrotated and rotated density
  nxmap0.init( xmap.cell(), xmap.grid_sampling(), gr1 );
  nxmap1.init( xmap.cell(), xmap.grid_sampling(), gr1 );
  // init 3 nxmaps for correl, rmsdev, mask dev
  correl.init( xmap.cell(), xmap.grid_sampling(), gr1 );
  rmsdev.init( xmap.cell(), xmap.grid_sampling(), gr1 );
  mskdev.init( xmap.cell(), xmap.grid_sampling(), gr1 );
  // make 3 aliases and 1 temporary nxmap
  clipper::NXmap<float>& x0 = correl;
  clipper::NXmap<float>& x1 = rmsdev;
  clipper::NXmap<float>& x11 = mskdev;
  clipper::NXmap<float> x00( xmap.cell(), xmap.grid_sampling(), gr1 );

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
  // alias the maps
  clipper::NXmap<float>& prd = x0;
  clipper::NXmap<float>& sub = x1;
  clipper::NXmap<float>& cov = x00;
  clipper::NXmap<float>& var = x11;
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
    double rad2 = (inx.coord_orth()-nxop.src()).lengthsq();
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
  for ( MRI inx = nxmap0.first(); !inx.last(); inx.next() ) {
    double dc = clipper::Util::max( (correl[inx]-c1)/c2, 0.0 );
    double pncs = 1.0/(1.0+exp(-0.5*dc*dc)/pncs0);
    double pmap = 1.0 - pncs;
    double vmap = (dmap*(pmap*dmap+pncs*dncs)) / (pncs*(dmap-dncs));
    mskdev[inx] = sqrt(vmap);             // map of standard deviations
  }

  // return stats
  std::vector<double> info(5);
  info[0] = dncs; info[1] = dmap;
  info[2] = cncs; info[3] = c1; info[4] = c2;
  return info;
}


/*
 Refine operators using masked density
*/
Local_rtop NCSfind::masked_refine_from_map( const clipper::Xmap<float>& xmap, const Local_rtop& nxop, const double& map_radius, const double& local_radius )
{
  typedef clipper::NXmap<float>::Map_reference_index MRI;

  // determine NCS correlation
  clipper::NXmap<float> nxmap0, nxmap1, correl, rmsdev, mskdev;
  std::vector<double> info = NCSfind::local_variance_from_map( nxmap0, nxmap1, correl, rmsdev, mskdev, xmap, nxop, map_radius, local_radius, 0.01 );

  // check that the mask doesn't exceed the sample sphere
  clipper::Coord_orth null( clipper::Coord_orth::null() );
  Local_rtop newop( clipper::Rotation::null(), null, null );
  if ( info[3] > 0.15 ) return newop;

  // calculate map centre of mass
  clipper::Coord_orth cmass( 0.0, 0.0, 0.0 );
  double mass(0.0), weight(0.0);
  for ( MRI inx = nxmap0.first(); !inx.last(); inx.next() ) {
    weight = 1.0/clipper::Util::sqr(mskdev[inx]);
    mass  += weight;
    cmass += weight*inx.coord_orth();  // inverse variance weighted centre
  }
  cmass = (1.0/mass) * cmass;

  // refine operator using mask
  Target_fn_xmap_mask_rtop mr( xmap, xmap, mskdev, 0.01, 0.1 );
  Local_rtop refop = mr.refine( nxop );
  newop = Local_rtop( refop.rot(), cmass, refop.rtop_orth()*cmass );

  return newop;
}


/*
  Geometrical subroutine to find the distance of closest approach of two lines.
*/
double NCSfind::minimum_distance_between_lines( const clipper::Coord_orth& a0, const clipper::Coord_orth& b0, const clipper::Coord_orth& a1, const clipper::Coord_orth& b1 )
{
  double r00 = 2.0*(b0*b0);
  double r11 = 2.0*(b1*b1);
  double r01 = -2.0*(b0*b1);
  double r10 = -2.0*(b0*b1);
  double q0 = 2.0*(b0*(a0-a1));
  double q1 = 2.0*(b1*(a1-a0));
  double c0 = (r01*q1-r11*q0)/(r00*r11-r01*r10);
  double c1 = (r10*q0-r00*q1)/(r00*r11-r01*r10);
  clipper::Coord_orth p0 = a0 + c0*b0;
  clipper::Coord_orth p1 = a1 + c1*b1;
  return sqrt( (p0-p1).lengthsq() );
}
