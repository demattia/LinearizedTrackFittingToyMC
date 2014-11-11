//
//  Track.h
//  LinearizedTrackFitting
//
//  Created by Luciano Ristori on 4/25/14.
//

#ifndef __LinearizedTrackFitting__Track__
#define __LinearizedTrackFitting__Track__

class Track {
    
private:
    
    void init(double x0_, double y0_, double z0_, double invPt_, double eta_, double phi_);
    
public:
    
    // primary parameters
    
    double x0, y0, z0; // origin
    double invPt, eta, phi; // pt in GeV/c, phi [-Pi, +Pi]
    
    // derived parameters
    
    //double px, py, pz, pt;
    //double R; // radius of curvature
    double c; // curvature
    double charge; // sign of invPt
    double k; // constant for z evolution = R*pz/pt
    double cotTheta; // cotangent of polar angle
    
    const double cMin = 0.001; // minimum curvature - use 1st order approx if less than that
    
    
    // constructors
    
    Track(); // default - generates random track
    
    Track(Geometry &g); // specific geometry
    
    Track(int special); // special track
    
    Track(double x0_, double y0_, double z0_, double invPt_, double eta_, double phi_); // specific track constructor
    
    void print(); // dumps text to cout
    
    // calculate coordinates of intersections with detector planes
    
    bool xzBarrel(double yDet, double &x, double &z); // yDet is the position of the barrel plane (parallel to xz)
    bool xyDisc(Geometry &g, int iDisc, double &X, double &R); // disc plane (parallel to xy)
    bool phizBarrel(Geometry &g, int iBarrel, double &hitPhi, double &hitZ);// hitPhi is a linear coordinate measured along phi
    bool getHit(Geometry &g, int iLayer, double &x, double &z);// find hit in generic plane
    
    // rigidly rotate track around z
    
    void rotate(double phiRot);
};

#endif /* defined(__LinearizedTrackFitting__Track__) */
