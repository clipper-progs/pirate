// Clipper pirate
/* Copyright 2003-2005 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-contrib.h>
#include "simulate-lib.h"
#include "pirate-lib.h"
#include <time.h>

extern "C" {
#include <stdlib.h>
}


int main( int argc, char** argv )
{
  CCP4Program prog( "cpirate", "0.4.9", "$Date: 2005/08/05" );

  std::cout << "\nCopyright 2003-2005 Kevin Cowtan and University of York\n";
  std::cout << "All rights reserved. Please reference:\n";
  std::cout << " Cowtan K. (2000) Acta Cryst. D56, 1612-1621.\n";

  int tm = time( NULL );
  //if ( tm > 1136073600 ) clipper::Message::message( clipper::Message_warn( "This program has expired: download a new version" ) );

  // defaults
  clipper::String title;
  clipper::String ipfile_ref = "csimulate.mtz";
  clipper::String ipfile_wrk = "NONE";
  clipper::String ipcol_ref_fo = "/*/*/[FP,SIGFP]";
  clipper::String ipcol_ref_hl = "/*/*/FC";
  clipper::String ipcol_wrk_fo = "NONE";
  clipper::String ipcol_wrk_hl = "NONE";
  clipper::String ipcol_wrk_fr = "NONE";
  clipper::String ipcol_wrk_xx = "NONE";
  clipper::String ipfile_ha = "NONE";
  clipper::String opfile = "pirate.mtz";
  clipper::String opcol = "pirate";
  clipper::String opcol_hl = "NONE";
  clipper::String opcol_fc = "NONE";
  std::vector<Local_rtop> ncsops, ncsnul;
  bool autocontent = false;
  bool automapllk  = false;
  bool autoncsllk  = false;
  bool unbias   = false;
  bool evaluate = false;
  bool strictfr = false;
  double filter_rad1 = 9.0;    // Initial radius of local mean in Angstroms
  double filter_rad2 = 3.0;    // Final   radius of local mean in Angstroms
  double ncs_radius = 6.0;     // Radius of NCS local correlation
  double ncs_volume = 1.0;     // Number of ASUs to allow for NCS
  double wt_expllk = 1.00;     // Scale factor for exp llk
  double wt_mapllk = 0.10;     // Scale factor for map llk
  double wt_ncsllk = 0.05;     // Scale factor for ncs llk
  double wt_ramp = 0.0;        // Ramp factor for weights
  double skew_mean = 0.0;      // skew factor for 1st moment
  double skew_sigm = 0.0;      // skew factor for 2nd moment
  double res_in = 1.0;         // Resolution limit
  double oversampling = 1.5;
  int ncycles = 1;
  int seed = 54321;
  int verbose = 0;

  // fixed values
  double centre_radius = 2.0;
  const int n_bins_mean = 9;
  const int n_bins_sigm = 9;
  const int bin_min = 1000;

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if        ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-mtzin-ref" ) {
      if ( ++arg < args.size() ) ipfile_ref = args[arg];
    } else if ( args[arg] == "-mtzin-wrk" ) {
      if ( ++arg < args.size() ) ipfile_wrk = args[arg];
    } else if ( args[arg] == "-pdbin-ha" ) {
      if ( ++arg < args.size() ) ipfile_ha = args[arg];
    } else if ( args[arg] == "-mtzout" ) {
      if ( ++arg < args.size() ) opfile = args[arg];
    } else if ( args[arg] == "-colin-ref-fo" ) {
      if ( ++arg < args.size() ) ipcol_ref_fo = args[arg];
    } else if ( args[arg] == "-colin-ref-hl" ) {
      if ( ++arg < args.size() ) ipcol_ref_hl = args[arg];
    } else if ( args[arg] == "-colin-wrk-fo" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fo = args[arg];
    } else if ( args[arg] == "-colin-wrk-hl" ) {
      if ( ++arg < args.size() ) ipcol_wrk_hl = args[arg];
    } else if ( args[arg] == "-colin-wrk-free" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fr = args[arg];
    } else if ( args[arg] == "-colin-wrk-solved-fc" ) {
      if ( ++arg < args.size() ) ipcol_wrk_xx = args[arg];
    } else if ( args[arg] == "-colout" ) {
      if ( ++arg < args.size() ) opcol = args[arg];
    } else if ( args[arg] == "-colout-hl" ) {
      if ( ++arg < args.size() ) opcol_hl = args[arg];
    } else if ( args[arg] == "-colout-fc" ) {
      if ( ++arg < args.size() ) opcol_fc = args[arg];
    } else if ( args[arg] == "-ncycles" ) {
      if ( ++arg < args.size() ) ncycles = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-weight-expllk" ) {
      if ( ++arg < args.size() ) wt_expllk = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-weight-mapllk" ) {
      if ( ++arg < args.size() ) wt_mapllk = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-weight-ncsllk" ) {
      if ( ++arg < args.size() ) wt_ncsllk = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-weight-ramp" ) {
      if ( ++arg < args.size() ) wt_ramp = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-stats-radius" ) {
      if ( ++arg < args.size() ) {
	filter_rad1 = clipper::String(args[arg]).split(",").front().f();
	filter_rad2 = clipper::String(args[arg]).split(",").back().f();
      }
    } else if ( args[arg] == "-skew-content" ) {
      if ( ++arg < args.size() ) {
	skew_mean = clipper::String(args[arg]).split(",").front().f();
	skew_sigm = clipper::String(args[arg]).split(",").back().f();
      }
    } else if ( args[arg] == "-auto-content" ) {
      autocontent = true;
    } else if ( args[arg] == "-auto-mapllk" ) {
      automapllk  = true;
    } else if ( args[arg] == "-auto-ncsllk" ) {
      autoncsllk  = true;
    } else if ( args[arg] == "-unbias" ) {
      unbias   = true;
    } else if ( args[arg] == "-evaluate" ) {
      evaluate = true;
    } else if ( args[arg] == "-strict-free" ) {
      strictfr = true;
    } else if ( args[arg] == "-resolution" ) {
      if ( ++arg < args.size() ) res_in = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-ncs-radius" ) {
      if ( ++arg < args.size() ) ncs_radius = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-ncs-volume" ) {
      if ( ++arg < args.size() ) ncs_volume = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-ncs-operator" ) {
      if ( ++arg < args.size() ) {
	std::vector<clipper::String> op = clipper::String(args[arg]).split(",");
	if ( op.size() == 9 ) {
	  clipper::Euler_ccp4 eul( clipper::Util::d2rad(op[0].f64()),
				   clipper::Util::d2rad(op[1].f64()),
				   clipper::Util::d2rad(op[2].f64()) );
	  clipper::Coord_orth src( op[3].f64(), op[4].f64(), op[5].f64() );
	  clipper::Coord_orth tgt( op[6].f64(), op[7].f64(), op[8].f64() );
	  clipper::Rotation rot( eul );
	  ncsops.push_back( Local_rtop( rot, src, tgt ) );
	} else {
	  std::cout << "\nInvalid ncs operator:\t" << args[arg] << "\n";
	  args.clear();
	}
      }
    } else if ( args[arg] == "-seed" ) {
      if ( ++arg < args.size() ) seed = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "\nUnrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "\nUsage: cpirate\n\t-mtzin-ref <filename>\t\tCOMPULSORY\n\t-mtzin-wrk <filename>\t\tCOMPULSORY\n\t-pdbin-ha <filename>\n\t-mtzout <filename>\n\t-colin-ref-fo <colpath>\n\t-colin-ref-hl <colpath>\n\t-colin-wrk-fo <colpath>\t\tCOMPULSORY\n\t-colin-wrk-hl <colpath>\t\tCOMPULSORY\n\t-colin-wrk-free <colpath>\n\t-colout <colpath>\n\t-colout-hl <colpath>\n\t-colout-fc <colpath>\n\t-ncycles <cycles>\n\t-weight-expllk <weight>\n\t-weight-mapllk <weight>\n\t-weight-ncsllk <weight>\n\t-weight-ramp <weight>\n\t-resolution <resolution/A>\n\t-stats-radius <radius/A>,<radius/A>\n\t-skew-content <factor>,<factor>\n\t-auto-content\n\t-auto-mapllk\n\t-auto-ncsllk\n\t-unbias\n\t-evaluate\n\t-strict-free\n\t-ncs-volume <factor>\n\t-ncs-radius <radius>\n\t-ncs-operator <alpha>,<beta>,<gamma>,<x>,<y>,<z>,<x>,<y>,<z>\n\t-seed <seed>\nAn input mtz are specified for both reference and work structures.\nFor the reference structure, F's and calc HL coefficients are required.\nFor the work structure, F's and HL coefficients are required.\nAll output data and the reference structure inputs default correctly.\n";
    exit(1);
  }

  if ( opcol_hl == "NONE" ) opcol_hl = opcol;
  if ( opcol_fc == "NONE" ) opcol_fc = opcol;
  if ( evaluate ) oversampling = 1.25;
  if ( wt_ramp == 0.0 ) wt_ramp = unbias ? 1.0 : 2.0;

  // numbers to output
  clipper::Resolution resol;
  clipper::CCP4MTZfile mtzfile;
 
  // Get resolution for calculation
  mtzfile.open_read( ipfile_ref );
  double res_ref = clipper::Util::max( mtzfile.resolution().limit(), res_in );
  mtzfile.close_read();
  mtzfile.open_read( ipfile_wrk );
  double res_wrk = clipper::Util::max( mtzfile.resolution().limit(), res_in );
  mtzfile.close_read();
  resol = clipper::Resolution( clipper::Util::max( res_ref, res_wrk ) );
  if ( res_ref > res_wrk ) std::cout << "\nWARNING: resolution of work structure truncated to reference:\n Ref: " << res_ref << " Wrk: " << res_wrk << std::endl;

  // Get reference reflection data
  clipper::HKL_info hkls_ref;
  mtzfile.open_read( ipfile_ref );
  hkls_ref.init( mtzfile.spacegroup(), mtzfile.cell(), resol, true );
  clipper::HKL_data<clipper::data32::F_sigF> ref_f( hkls_ref );
  clipper::HKL_data<clipper::data32::ABCD> ref_hl( hkls_ref );
  mtzfile.import_hkl_data( ref_f, ipcol_ref_fo );
  mtzfile.import_hkl_data( ref_hl, ipcol_ref_hl );
  mtzfile.close_read();

  // Get work reflection data
  clipper::HKL_info hkls;
  mtzfile.open_read( ipfile_wrk );
  hkls.init( mtzfile.spacegroup(), mtzfile.cell(), resol, true );
  clipper::HKL_data<clipper::data32::F_sigF> fobs( hkls );
  clipper::HKL_data<clipper::data32::ABCD>   abcd( hkls );
  clipper::HKL_data<clipper::data32::Flag>   flag( hkls );
  clipper::HKL_data<clipper::data32::F_phi>  fphi( hkls );
  mtzfile.import_hkl_data( fobs, ipcol_wrk_fo );
  mtzfile.import_hkl_data( abcd, ipcol_wrk_hl );
  if ( ipcol_wrk_fr != "NONE" ) mtzfile.import_hkl_data( flag, ipcol_wrk_fr );
  if ( ipcol_wrk_xx != "NONE" ) mtzfile.import_hkl_data( fphi, ipcol_wrk_xx );
  clipper::String oppath = mtzfile.assigned_paths()[0].notail() + "/";
  mtzfile.close_read();

  // apply free flag
  clipper::HKL_data<clipper::data32::F_sigF> fsig = fobs;
  { using namespace clipper::datatypes;  // 'using' for irix6.5
    fsig.mask( flag != 0 );
  }

  // apply llk scale
  clipper::HKL_info::HKL_reference_index ih;
  for ( ih = abcd.first(); !ih.last(); ih.next() ) {
    if ( abcd[ih].missing() ) {
      abcd[ih].a() = abcd[ih].b() = abcd[ih].c() = abcd[ih].d() = 0.0;
    } else {
      abcd[ih].a() *= wt_expllk;
      abcd[ih].b() *= wt_expllk;
      abcd[ih].c() *= wt_expllk;
      abcd[ih].d() *= wt_expllk;
    }
  }

  // get NCS operators for heavy atom (HA) file
  if ( ipfile_ha != "NONE" ) {
    // atomic model
    clipper::MMDBManager mmdb;
    mmdb.SetFlag( ::mmdb::MMDBF_AutoSerials | ::mmdb::MMDBF_IgnoreDuplSeqNum );
    mmdb.ReadPDBASCII( (char*)ipfile_ha.c_str() );

    // get a list of all the atoms
    clipper::mmdb::PPCAtom psel;
    int hndl, nsel;
    hndl = mmdb.NewSelection();
    mmdb.SelectAtoms( hndl, 0, 0, ::mmdb::SKEY_NEW );
    mmdb.GetSelIndex( hndl, psel, nsel );
    clipper::MMDBAtom_list atoms( psel, nsel );
    mmdb.DeleteSelection( hndl );

    // make a map
    clipper::HKL_data<clipper::data32::F_phi> fphi( hkls );
    clipper::HKL_data<clipper::data32::Phi_fom> phiw( hkls );
    phiw.compute( abcd, clipper::data32::Compute_phifom_from_abcd() );
    fphi.compute( fobs, phiw, clipper::data32::Compute_fphi_from_fsigf_phifom() );
    clipper::Grid_sampling grid( hkls.spacegroup(), hkls.cell(), hkls.resolution(), oversampling );
    clipper::Xmap<float> xmap( hkls.spacegroup(), hkls.cell(), grid );
    xmap.fft_from( fphi );

    // search
    NCSfind ncsatom( 3.0, 0.1 );
    ncsops = ncsatom( atoms, xmap );

    // refine
    NCSfind ncsref( 3.0, 0.1 );
    const clipper::Spacegroup& spgr = hkls.spacegroup();
    const clipper::Cell&       cell = hkls.cell();
    double mrad = ncs_volume*pow(cell.volume()/double(spgr.num_symops()),0.333);
    ncsops = Local_rtop::exclude_inverse( ncsops, spgr, cell );
    Local_rtop rtop;
    std::vector<Local_rtop> ncsops_tmp;
    for ( int i = 0; i < ncsops.size(); i++ ) {
      rtop = ncsref.masked_refine_from_map( xmap, ncsops[i], mrad, ncs_radius );
      if ( !rtop.is_null() ) ncsops_tmp.push_back( rtop );
    }
    ncsops = Local_rtop::include_inverse( ncsops_tmp, spgr, cell );
    ncsops = Local_rtop::tidy( ncsops, spgr, cell );

    // output NCS operators
    for ( int i = 0; i < ncsops.size(); i++ )
      std::cout << "\nNCS operator: " << i << "\n Rotation: " << ncsops[i].rot().polar_ccp4().format() << "\n Source:   " << ncsops[i].src().format() << "\n Target:   " << ncsops[i].tgt().format() << std::endl;
  }

  // other params
  srand( seed );

  // result objects
  clipper::HKL_data<clipper::data32::F_sigF> sim_f( hkls_ref );
  clipper::HKL_data<clipper::data32::ABCD> sim_hl( hkls_ref );
  clipper::HKL_data<clipper::data32::ABCD> abcd_new( hkls );    

  // DO INITIAL MAP SIMULATION

  MapSimulate mapsim( 100, 20 );
  mapsim( sim_f, sim_hl, ref_f, ref_hl, fobs, abcd );
  double filter_radius = filter_rad1;

  // DO AUTOMATIC PARAMETER DETERMINATION

  // content skewing
  if ( autocontent ) {
    double y[] = { -1.0, -1.0, -1.0 };
    double x[] = { -0.5,  0.0,  0.5 };
    for ( int ccyc = 0; ccyc < 5; ccyc++ ) {
      for ( int cpar = 0; cpar < 3; cpar++ ) if ( y[cpar] < 0.0 ) {
	skew_mean = x[cpar];
	Refine_HL_simulate phaseref( unbias, centre_radius, filter_radius,
				     ncs_radius, ncs_volume,
				     wt_mapllk, wt_ramp, skew_mean, skew_sigm,
				     n_bins_mean, n_bins_sigm, 12, 1.25 );
	y[cpar] = phaseref.correl_moment_1( fphi, abcd_new, fsig, fobs, abcd, sim_f, ref_hl, sim_hl, ncsnul );
      }
      if ( verbose >= 3 ) std::cout << "Skew " << ccyc << ":\t" << x[0] << " " << x[1] << " " << x[2] << " =\t" << y[0] << " " << y[1] << " " << y[2] << std::endl; 
      if      ( y[0] > y[1] && y[0] > y[2] ) { x[1] = x[0]; y[1] = y[0]; }
      else if ( y[2] > y[1] && y[2] > y[0] ) { x[1] = x[2]; y[1] = y[2]; }
      double dx = 0.25 * ( x[2] - x[0] );
      x[0] = x[1] - dx; x[2] = x[1] + dx;
      y[0] = y[2] = -1.0;
    }
    skew_mean = x[1];
    x[0] = -0.5; x[1] = 0.0; x[2] = 0.5;
    for ( int ccyc = 0; ccyc < 5; ccyc++ ) {
      for ( int cpar = 0; cpar < 3; cpar++ ) if ( y[cpar] < 0.0 ) {
	skew_sigm = x[cpar];
	Refine_HL_simulate phaseref( unbias, centre_radius, filter_radius,
				     ncs_radius, ncs_volume,
				     wt_mapllk, wt_ramp, skew_mean, skew_sigm,
				     n_bins_mean, n_bins_sigm, 12, 1.25 );
	y[cpar] = phaseref.correl_moment_2( fphi, abcd_new, fsig, fobs, abcd, sim_f, ref_hl, sim_hl, ncsnul );
      }
      if ( verbose >= 3 ) std::cout << "Skew " << ccyc << ":\t" << x[0] << " " << x[1] << " " << x[2] << " =\t" << y[0] << " " << y[1] << " " << y[2] << std::endl; 
      if      ( y[0] > y[1] && y[0] > y[2] ) { x[1] = x[0]; y[1] = y[0]; }
      else if ( y[2] > y[1] && y[2] > y[0] ) { x[1] = x[2]; y[1] = y[2]; }
      double dx = 0.25 * ( x[2] - x[0] );
      x[0] = x[1] - dx; x[2] = x[1] + dx;
      y[0] = y[2] = -1.0;
    }
    skew_sigm = x[1];
    std::cout << "\nAutomatic content fitting:\n Skew by: " << skew_mean << "," << skew_sigm << " (dense, ordered)\n";
  }

  // mapllk weight
  if ( automapllk ) {
    double y[] = { -1.0, -1.0, -1.0 };
    double x[] = { -3.0, -2.0, -1.0 };
    for ( int ccyc = 0; ccyc < 3; ccyc++ ) {
      for ( int cpar = 0; cpar < 3; cpar++ ) if ( y[cpar] < 0.0 ) {
	double wt_tmp = exp( x[cpar] );
	Refine_HL_simulate phaseref( unbias, centre_radius, filter_radius,
				     ncs_radius, ncs_volume,
				     wt_tmp, wt_ramp, skew_mean, skew_sigm,
				     n_bins_mean, n_bins_sigm, 12, 1.25 );
	phaseref( fphi, abcd_new, fsig, fobs, abcd, sim_f, ref_hl, sim_hl,
		  ncsnul );
	y[cpar] = phaseref.llk_gain_free();
      }
      if ( verbose >= 3 ) std::cout << "MapW " << ccyc << ":\t" << x[0] << " " << x[1] << " " << x[2] << " =\t" << y[0] << " " << y[1] << " " << y[2] << std::endl; 
      if      ( y[0] > y[1] && y[0] > y[2] ) { x[1] = x[0]; y[1] = y[0]; }
      else if ( y[2] > y[1] && y[2] > y[0] ) { x[1] = x[2]; y[1] = y[2]; }
      double dx = 0.25 * ( x[2] - x[0] );
      x[0] = x[1] - dx; x[2] = x[1] + dx;
      y[0] = y[2] = -1.0;
    }
    double wt_total = exp( x[1] );
    const double n = double( ncycles );
    wt_mapllk = wt_total / n;
    // upweight slightly when using multiple cycles
    const double dn = 0.5;
    wt_mapllk = wt_mapllk * (n+dn*n)/(n+dn);
    std::cout << "\nAutomatic likelihood weighting:\n mapllk-wt: " << wt_total << "\t" << wt_mapllk << " (total, per-cycle)\n";
  }

  // ncsllk weight
  if ( autoncsllk && ncsops.size() > 0 ) {
    double y[] = { -1.0, -1.0, -1.0 };
    double x[] = { -4.0, -3.0, -2.0 };
    for ( int ccyc = 0; ccyc < 3; ccyc++ ) {
      for ( int cpar = 0; cpar < 3; cpar++ ) if ( y[cpar] < 0.0 ) {
	double wt_tmp = exp( x[cpar] );
	Refine_HL_simulate phaseref( unbias, centre_radius, filter_radius,
				     ncs_radius, ncs_volume,
				     wt_tmp, wt_ramp, skew_mean, skew_sigm,
				     n_bins_mean, n_bins_sigm, 12, 1.25 );
	phaseref( fphi, abcd_new, fsig, fobs, abcd, ncsops );
	y[cpar] = phaseref.llk_gain_free();
      }
      if ( verbose >= 3 ) std::cout << "MapW " << ccyc << ":\t" << x[0] << " " << x[1] << " " << x[2] << " =\t" << y[0] << " " << y[1] << " " << y[2] << std::endl; 
      if      ( y[0] > y[1] && y[0] > y[2] ) { x[1] = x[0]; y[1] = y[0]; }
      else if ( y[2] > y[1] && y[2] > y[0] ) { x[1] = x[2]; y[1] = y[2]; }
      double dx = 0.25 * ( x[2] - x[0] );
      x[0] = x[1] - dx; x[2] = x[1] + dx;
      y[0] = y[2] = -1.0;
    }
    double wt_total = exp( x[1] );
    wt_ncsllk = 0.4 * wt_total;
    std::cout << "\nAutomatic likelihood weighting:\n ncsllk-wt: " << wt_total << "\t" << wt_ncsllk << " (base, adjusted)\n";
  }

  // DO FIRST PHASE IMPROVEMENT CYCLE

  Refine_HL_simulate phaseref;

  // do ncs refinement
  if ( ncsops.size() > 0 ) {
    phaseref=Refine_HL_simulate( unbias, centre_radius, filter_radius,
				 ncs_radius, ncs_volume,
				 wt_ncsllk, wt_ramp, skew_mean, skew_sigm,
				 n_bins_mean, n_bins_sigm, 12, oversampling );
    phaseref( fphi, abcd_new, fsig, fobs, abcd, ncsops );
    for ( ih = abcd.first(); !ih.last(); ih.next() ) abcd[ih] = abcd_new[ih];
    std::cout << "\nUnbiased results from NCS cycle:"
	      << "\n R-factor     : " << phaseref.r_factor_work()
	      << "\n Free R-factor: " << phaseref.r_factor_free()
	      << "\n E-correl     : " << phaseref.e_correl_work()
	      << "\n Free E-correl: " << phaseref.e_correl_free()
	      << "\n LgLkGain     : " << phaseref.llk_gain_work()
	      << "\n Free LgLkGain: " << phaseref.llk_gain_free()
	      << std::endl;
  }
  // do phase refinement calc
  if ( ncycles > 0 ) {
    phaseref=Refine_HL_simulate( unbias, centre_radius, filter_radius,
				 ncs_radius, ncs_volume,
				 wt_mapllk, wt_ramp, skew_mean, skew_sigm,
				 n_bins_mean, n_bins_sigm, 12, oversampling );
    if ( verbose >= 10 ) phaseref.debug();
    phaseref( fphi, abcd_new,
	      fsig, fobs, abcd, sim_f, ref_hl, sim_hl, ncsnul );
    // output stats
    std::cout << "\nUnbiased results from initial cycle:"
	      << "\n R-factor     : " << phaseref.r_factor_work()
	      << "\n Free R-factor: " << phaseref.r_factor_free()
	      << "\n E-correl     : " << phaseref.e_correl_work()
	      << "\n Free E-correl: " << phaseref.e_correl_free()
	      << "\n LgLkGain     : " << phaseref.llk_gain_work()
	      << "\n Free LgLkGain: " << phaseref.llk_gain_free()
	      << std::endl;
  }

  if ( evaluate ) goto done;

  // DO ADDITIONAL PHASE IMPROVEMENT CYCLES

  for ( int cyc = 1; cyc < ncycles; cyc++ ) {
    // update starting point
    for ( ih = abcd.first(); !ih.last(); ih.next() ) abcd[ih] = abcd_new[ih];
    filter_radius += (filter_rad2-filter_rad1)/double(ncycles-1);

    // update simulation (REMOVING THIS STEP DOESN'T MAKE MUCH DIFFERENCE)
    mapsim( sim_f, sim_hl, ref_f, ref_hl, fobs, abcd );
    // do phase refinement calc
    phaseref=Refine_HL_simulate( false, centre_radius, filter_radius,
				 ncs_radius, ncs_volume,
				 wt_mapllk, wt_ramp, skew_mean, skew_sigm,
				 n_bins_mean, n_bins_sigm, 12, oversampling );
    if ( strictfr )
      phaseref( fphi, abcd_new,
		fsig, fobs, abcd, sim_f, ref_hl, sim_hl, ncsnul );
    else
      phaseref( fphi, abcd_new,
		fobs, fobs, abcd, sim_f, ref_hl, sim_hl, ncsnul );
    // output stats
    std::cout << "\nBiased results from cycle " << cyc+1 << ":"
	      << "\n Biased R-factor     : " << phaseref.r_factor_work()
	      << "\n Biased Free R-factor: " << phaseref.r_factor_free()
	      << "\n Biased E-correl     : " << phaseref.e_correl_work()
	      << "\n Biased Free E-correl: " << phaseref.e_correl_free()
	      << "\n Biased LgLkGain     : " << phaseref.llk_gain_work()
	      << "\n Biased Free LgLkGain: " << phaseref.llk_gain_free()
	      << std::endl;
  }

  // DONE PHASE IMPROVEMENT

  // output new results
  if ( opcol_hl[0] != '/' ) opcol_hl = oppath + opcol_hl;
  if ( opcol_fc[0] != '/' ) opcol_fc = oppath + opcol_fc;
  mtzfile.open_append( ipfile_wrk, opfile );
  mtzfile.export_hkl_data( abcd_new, opcol_hl );
  mtzfile.export_hkl_data( fphi    , opcol_fc );
  mtzfile.close_append();

  // output new NCS operators
  std::cout << std::endl;
  for ( int i = 0; i < ncsops.size(); i++ )
    std::cout << "ncs-operator " << clipper::Util::rad2d(ncsops[i].rot().euler_ccp4().alpha()) << "," << clipper::Util::rad2d(ncsops[i].rot().euler_ccp4().beta()) << "," << clipper::Util::rad2d(ncsops[i].rot().euler_ccp4().gamma()) << "," << ncsops[i].src().x() << "," << ncsops[i].src().y() << "," << ncsops[i].src().z() << "," << ncsops[i].tgt().x() << "," << ncsops[i].tgt().y() << "," << ncsops[i].tgt().z() << std::endl;

 done:
  // debug output
  if ( verbose >= 1 )
    std::cout << "\nMap RMSDs:\n Ref calc: " << phaseref.rmsd_calc()
	      << "\n Ref simu: " << phaseref.rmsd_simu()
	      << "\n Wrk obs : " << phaseref.rmsd_work() << std::endl;
  if ( verbose >= 2 ) {
    std::cout << "\nHistograms:\n";
    for ( int k = 0; k < 50; k++ ) {
      clipper::Range<double> xrng = phaseref.hist_mod()(0,0);
      double x = ((double(k)+0.5)/50.0)*xrng.range() + xrng.min();
      std::cout << clipper::String(x).substr(0,6) << " ";
      for ( int i = 0; i < n_bins_mean; i+= n_bins_mean/2 )
	for ( int j = 0; j < n_bins_sigm; j+= n_bins_sigm/2 )
	  std::cout << clipper::String(int(phaseref.hist_mod()(i,j).y(x)),6) << " ";
      std::cout << std::endl;
    }
  }
  if ( verbose >= 3 ) {
    std::cout << "\nE/FOM bins:\n";
    for ( int i = 0; i < n_bins_mean; i++ ) {
      for ( int j = 0; j < n_bins_sigm; j++ )
	std::cout << "  " << clipper::String(int(phaseref.hist_raw()(i,j).sum()),4) << "/" << clipper::String(int(phaseref.hist_mod()(i,j).sum()),4);
      std::cout << std::endl;
    }
  }
}
