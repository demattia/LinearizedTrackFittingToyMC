//
//  Geometry.h
//  LinearizedTrackFitting
//
//  Created by Luciano Ristori on 5/1/14.
//  Copyright (c) 2014 Luciano Ristori. All rights reserved.
//

#ifndef __LinearizedTrackFitting__Geometry__
#define __LinearizedTrackFitting__Geometry__

#include <vector>
#include <cmath>
#include <string>

class Geometry {
    
public:
    
    const double Pi = 3.14159265;
    
    // track generation parameters
    
    double t_eta = 0.5;
    double t_phi = 0.;
    
    double t_invPt = 1./2.;// Gev/c^(-1)
    
    double t_x0 = 0.;
    double t_y0 = 0.;
    double t_z0 = 0.;
    
    double t_deltaX0 = 0.;
    double t_deltaY0 = 0.;
    double t_deltaZ0 = 0.1;
    
    double t_deltaEta = 0.1;
    double t_deltaPhi = 0.1;
    
    // detector geometry
 
    
    double ps_x = 30.e-6;
    double ps_z = 1.5e-3/sqrt(12.);
    double ss_x = 30.e-6;
    double ss_z = 0.05/sqrt(12.);
    
    
    static const int nBarrels = 6;
    const double rBarrel[nBarrels] = {0.2, 0.3, 0.5, 0.65, 0.9, 1.1};
    
    const double dx_b[nBarrels]= {ps_x,ps_x,ps_x,ss_x,ss_x,ss_x}; // error of x coordinate measurement
    const double dz_b[nBarrels]= {ps_z,ps_z,ps_z,ss_z,ss_z,ss_z} ; // error of z coordinate measurement
    
    static const int nDiscs = 5;
    const double zDisc[nDiscs] = {1.4, 1.6, 1.9, 2.2, 2.6};
    const double rBoundary = 0.5; // boundary between ps and ss modules
    
    // layer configuration
    
    static const int nLayers = 6;
    
    char detType[nLayers] = {'B','B','B','B','B','B'}; // barrel (B) or disc (D)
    int detInd[nLayers] = {0, 1, 2, 3, 4, 5};

    Geometry(){};
    
};

#endif /* defined(__LinearizedTrackFitting__Geometry__) */
