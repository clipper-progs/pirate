// Clipper pirate
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <ccp4_program.h>
#include "piralib.h"
#include "ccp4-extras.h"
#include <time.h>


int main( int argc, char** argv )
{
  CCP4program prog( "cpirate", "0.2.1", "$Date: 2004/09/01" );

  std::cout << "\nCopyright 2003-2004 Kevin Cowtan and University of York\n";
  std::cout << "All rights reserved. Please reference:\n";
  std::cout << " Cowtan K. (2000) Acta Cryst. D56, 1612-1621.\n";

  int tm = time( NULL );
  if ( tm < 1088636400 || tm > 1104537600 ) clipper::Message::message( clipper::Message_fatal( "This program has expired: download a new version" ) );

  // defaults
  clipper::String ipfile_ref = "csimulate.mtz";
  clipper::String ipfile_wrk = "NONE";
  clipper::String ipcol_ref_fo = "/*/*/[FP,SIGFP]";
  clipper::String ipcol_ref_hl = "/*/*/FC";
  clipper::String ipcol_wrk_fo = "NONE";
  clipper::String ipcol_wrk_hl = "NONE";
  clipper::String ipcol_wrk_fr = "NONE";
  clipper::String opfile = "pirate.mtz";
  clipper::String opcol = "pirate";
  clipper::String opcol_hl = "NONE";
  clipper::String opcol_fc = "NONE";
  clipper::String biasmode = "bias";  // DEPRECATED
  bool acontent = false;
  bool unbias   = false;
  bool evaluate = false;
  bool strictfr = false;
  double filter_rad1 = 9.0;    // Initial radius of local mean in Angstroms
  double filter_rad2 = 3.0;    // Final   radius of local mean in Angstroms
  double wt_expllk = 1.0;      // Scale factor for exp llk
  double wt_mapllk = 0.1;      // Scale factor for map llk
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
  CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-mtzin-ref" ) {
      if ( ++arg < args.size() ) ipfile_ref = args[arg];
    } else if ( args[arg] == "-mtzin-wrk" ) {
      if ( ++arg < args.size() ) ipfile_wrk = args[arg];
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
      acontent = true;
    } else if ( args[arg] == "-unbias" ) {
      unbias   = true;
    } else if ( args[arg] == "-evaluate" ) {
      evaluate = true;
    } else if ( args[arg] == "-strict-free" ) {
      strictfr = true;
    } else if ( args[arg] == "-resolution" ) {
      if ( ++arg < args.size() ) res_in = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-seed" ) {
      if ( ++arg < args.size() ) seed = clipper::String(args[arg]).i();
    } else {
      std::cout << "\nUnrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "\nUsage: cpirate\n\t-mtzin-ref <filename>\t\tCOMPULSORY\n\t-mtzin-wrk <filename>\t\tCOMPULSORY\n\t-mtzout <filename>\n\t-colin-ref-fo <colpath>\n\t-colin-ref-hl <colpath>\n\t-colin-wrk-fo <colpath>\t\tCOMPULSORY\n\t-colin-wrk-hl <colpath>\t\tCOMPULSORY\n\t-colin-wrk-free <colpath>\n\t-colout <colpath>\n\t-colout-hl <colpath>\n\t-colout-fc <colpath>\n\t-ncycles <cycles>\n\t-weight-expllk <weight>\n\t-weight-mapllk <weight>\n\t-weight-ramp <weight>\n\t-resolution <resolution/A>\n\t-stats-radius <radius/A>,<radius/A>\n\t-skew-content <factor>,<factor>\n\t-auto-content\n\t-unbias\n\t-evaluate\n\t-strict-free\n\t-seed <seed>\nAn input mtz are specified for both reference and work structures.\nFor the reference structure, F's and calc HL coefficients are required.\nFor the work structure, F's and HL coefficients are required.\nAll output data and the reference structure inputs default correctly.\n";
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
  if ( res_ref > res_wrk ) std::cout << "\nWARNING: resolution of work structure truncated to reference:\n Ref: " << res_ref << " Wrk: " << res_wrk << "\n";

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
  mtzfile.import_hkl_data( fobs, ipcol_wrk_fo );
  mtzfile.import_hkl_data( abcd, ipcol_wrk_hl );
  if ( ipcol_wrk_fr != "NONE" ) mtzfile.import_hkl_data( flag, ipcol_wrk_fr );
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

  // other params
  srand( seed );

  // result objects
  clipper::HKL_data<clipper::data32::F_sigF> sim_f( hkls_ref );
  clipper::HKL_data<clipper::data32::ABCD> sim_hl( hkls_ref );
  clipper::HKL_data<clipper::data32::F_phi> fphi( hkls );
  clipper::HKL_data<clipper::data32::ABCD> abcd_new( hkls );    

  // DO INITIAL MAP SIMULATION

  MapSimulate mapsim( 100, 20 );
  mapsim( sim_f, sim_hl, ref_f, ref_hl, fobs, abcd );
  double filter_radius = filter_rad1;

  // DO AUTOMATIC PARAMETER DETERMINATION

  if ( acontent ) {
    double y[] = { -1.0, -1.0, -1.0 };
    double x[] = { -0.5,  0.0,  0.5 };
    for ( int ccyc = 0; ccyc < 3; ccyc++ ) {
      for ( int cpar = 0; cpar < 3; cpar++ ) if ( y[cpar] < 0.0 ) {
	skew_mean = skew_sigm = x[cpar];
	Refine_HL_simulate phaseref( unbias, centre_radius, filter_radius,
				     wt_mapllk, wt_ramp, skew_mean, skew_sigm,
				     n_bins_mean, n_bins_sigm, 12, 1.25 );
	phaseref( fphi, abcd_new, fsig, fobs, abcd, sim_f, ref_hl, sim_hl );
	y[cpar] = phaseref.e_correl_free();
      }
      if      ( y[0] > y[1] && y[0] > y[2] ) { x[1] = x[0]; y[1] = y[0]; }
      else if ( y[2] > y[1] && y[2] > y[0] ) { x[1] = x[2]; y[1] = y[2]; }
      double dx = 0.25 * ( x[2] - x[0] );
      x[0] = x[1] - dx; x[2] = x[1] + dx;
      y[0] = y[2] = -1.0;
    }
    double skew_base = x[1];
    x[0] = -0.1; x[1] = 0.0; x[2] = 0.1;
    for ( int ccyc = 0; ccyc < 2; ccyc++ ) {
      for ( int cpar = 0; cpar < 3; cpar++ ) if ( y[cpar] < 0.0 ) {
	skew_mean = skew_base - x[cpar];
	skew_sigm = skew_base + x[cpar];
	Refine_HL_simulate phaseref( unbias, centre_radius, filter_radius,
				     wt_mapllk, wt_ramp, skew_mean, skew_sigm,
				     n_bins_mean, n_bins_sigm, 12, 1.25 );
	phaseref( fphi, abcd_new, fsig, fobs, abcd, sim_f, ref_hl, sim_hl );
	y[cpar] = phaseref.e_correl_free();
      }
      if      ( y[0] > y[1] && y[0] > y[2] ) { x[1] = x[0]; y[1] = y[0]; }
      else if ( y[2] > y[1] && y[2] > y[0] ) { x[1] = x[2]; y[1] = y[2]; }
      double dx = 0.25 * ( x[2] - x[0] );
      x[0] = x[1] - dx; x[2] = x[1] + dx;
      y[0] = y[2] = -1.0;
    }
    skew_mean = skew_base - x[1];
    skew_sigm = skew_base + x[1];
    std::cout << "\nAutomatic content fitting:\n Skew by: " << skew_mean << "," << skew_sigm << " (dense, ordered)\n";
  }


  // DO FIRST PHASE IMPROVEMENT CYCLE

  // do phase refinement calc
  Refine_HL_simulate phaseref( unbias, centre_radius, filter_radius,
			       wt_mapllk, wt_ramp, skew_mean, skew_sigm,
			       n_bins_mean, n_bins_sigm, 12, oversampling );
  if ( verbose >= 10 ) phaseref.debug();
  phaseref( fphi, abcd_new, fsig, fobs, abcd, sim_f, ref_hl, sim_hl );
  // output stats
  std::cout << "\nUnbiased results from initial cycle:"
	    << "\n R-factor     : " << phaseref.r_factor_work()
	    << "\n Free R-factor: " << phaseref.r_factor_free()
	    << "\n E-correl     : " << phaseref.e_correl_work()
	    << "\n Free E-correl: " << phaseref.e_correl_free() <<"\n";

  if ( evaluate ) goto done;

  // DO ADDITIONAL PHASE IMPROVEMENT CYCLES

  for ( int cyc = 1; cyc < ncycles; cyc++ ) {
    // update starting point
    for ( ih = abcd.first(); !ih.last(); ih.next() ) abcd[ih] = abcd_new[ih];
    filter_radius += (filter_rad2-filter_rad1)/double(ncycles-1);

    // update simulation (REMOVING THIS STEP DOESN'T MAKE MUCH DIFFERENCE)
    mapsim( sim_f, sim_hl, ref_f, ref_hl, fobs, abcd );
    // do phase refinement calc
    Refine_HL_simulate phaseref( false, centre_radius, filter_radius,
				 wt_mapllk, wt_ramp, skew_mean, skew_sigm,
				 n_bins_mean, n_bins_sigm, 12, oversampling );
    if ( strictfr )
      phaseref( fphi, abcd_new, fsig, fobs, abcd, sim_f, ref_hl, sim_hl );
    else
      phaseref( fphi, abcd_new, fobs, fobs, abcd, sim_f, ref_hl, sim_hl );
    // output stats
    std::cout << "\nBiased results from cycle " << cyc+1 << ":"
	      << "\n Biased R-factor     : " << phaseref.r_factor_work()
	      << "\n Biased Free R-factor: " << phaseref.r_factor_free()
	      << "\n Biased E-correl     : " << phaseref.e_correl_work()
	      << "\n Biased Free E-correl: " << phaseref.e_correl_free() <<"\n";
  }

  // DONE PHASE IMPROVEMENT

  // output new results
  if ( opcol_hl[0] != '/' ) opcol_hl = oppath + opcol_hl;
  if ( opcol_fc[0] != '/' ) opcol_fc = oppath + opcol_fc;
  mtzfile.open_append( ipfile_wrk, opfile );
  mtzfile.export_hkl_data( abcd_new, opcol_hl );
  mtzfile.export_hkl_data( fphi    , opcol_fc );
  mtzfile.close_append();

 done:
  // debug output
  if ( verbose >= 1 )
    std::cout << "\nMap RMSDs:\n Ref calc: " << phaseref.rmsd_calc()
	      << "\n Ref simu: " << phaseref.rmsd_simu()
	      << "\n Wrk obs : " << phaseref.rmsd_work() << "\n";
  if ( verbose >= 2 ) {
    std::cout << "\nHistograms:\n";
    for ( int k = 0; k < 50; k++ ) {
      clipper::Range<double> xrng = phaseref.hist_mod()(0,0);
      double x = ((double(k)+0.5)/50.0)*xrng.range() + xrng.min();
      std::cout << clipper::String(x).substr(0,6) << " ";
      for ( int i = 0; i < n_bins_mean; i+= n_bins_mean/2 )
	for ( int j = 0; j < n_bins_sigm; j+= n_bins_sigm/2 )
	  std::cout << clipper::String(int(phaseref.hist_mod()(i,j).y(x)),6) << " ";
      std::cout << "\n";
    }
  }
  if ( verbose >= 3 ) {
    std::cout << "\nE/FOM bins:\n";
    for ( int i = 0; i < n_bins_mean; i++ ) {
      for ( int j = 0; j < n_bins_sigm; j++ )
	std::cout << "  " << clipper::String(int(phaseref.hist_raw()(i,j).sum()),4) << "/" << clipper::String(int(phaseref.hist_mod()(i,j).sum()),4);
      std::cout << "\n";
    }
  }
}
