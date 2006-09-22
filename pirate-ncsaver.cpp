/*! \file pirate-ncsaver.cpp pirate library */
/* Copyright 2003-2006 Kevin Cowtan & University of York all rights reserved */

#include "pirate-ncsaver.h"


void Xmap_ncs::restrain_ncs( const clipper::Xmap<float>& xmap, const Local_rtop& nxop, const double& map_radius, const double& local_radius )
{
  typedef clipper::NXmap<float>::Map_reference_index MRI;

  // determine NCS correlation
  clipper::NXmap<float> nxmap0, nxmap1, correl, rmsdev, mskdev;
  std::vector<double> info = NCSfind::local_variance_from_map( nxmap0, nxmap1, correl, rmsdev, mskdev, xmap, nxop, map_radius, local_radius, 0.01 );

  //std::cout << "ncsc " << info[2] << " " << info[3] << " " << info[4] << " " << info[1] << " " << info[0] << "\n";

  // calculate map constraint
  clipper::Coord_grid offset = xmap.coord_map( nxmap0.coord_orth( clipper::Coord_map(0.0,0.0,0.0) ) ).coord_grid();
  for ( MRI inx = nxmap0.first(); !inx.last(); inx.next() ) {
    clipper::Coord_grid xcoord = inx.coord() + offset;
    Gaussian_probability_1d g1 = get_data(xcoord);
    Gaussian_probability_1d g2( nxmap1[inx], mskdev[inx] );
    set_data( xcoord, g1*g2 );
  }
}
