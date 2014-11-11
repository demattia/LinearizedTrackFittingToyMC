//
//  testRegion.cpp
//  LinearizedTrackFitting
//
//  Created by Luciano Ristori on 5/23/14.
//  Copyright (c) 2014 Luciano Ristori. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <string>
#include <fstream>

#include "Eigen/Eigenvalues"


#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

#include "Geometry.h"
#include "Track.h"
#include "Statistics.h"


using namespace std;
using namespace Eigen;

extern     std::default_random_engine generator;
extern     std::uniform_real_distribution<double> distribution;
extern     std::normal_distribution<double> gauss;

extern string pathToMyDirectory;
bool verbose = false;

const int nTrackParameters = 4; // number of track parameters
const int nVars = 12; // number of coordinates
const int nZeroConstr = 0; // number of constraint values to be ignored


//#include "testRegion.h"

int testRegion(long int nTracks){
    
     
     
    // open matrix file and read V and D arrays
    cout<< "opening matrixVD.txt for reading" << endl;
    
    std::ifstream infile;
    infile.open(pathToMyDirectory + "matrixVD.txt");
    if(!infile) {
        cout << "error opening matrixVD.txt" << endl;
        return -1;
    }
    
    // read transformation matrix V from file
    
    double x;
    
    MatrixXd V(nVars,nVars);
    
    for(int i = 0; i != nVars; ++i)
        for(int j = 0; j != nVars; ++j){
            infile >> x;
            V(i,j) = x;
        }
    cout << "V:"<< endl;
    cout << setprecision(4) << V << endl;
    
    // read transformation matrix D from file
    
    MatrixXd D(nTrackParameters,nVars);
    
    for(int i = 0; i != nTrackParameters; ++i)
        for(int j = 0; j != nVars; ++j){
            infile >> x;
            if(j < nZeroConstr)D(i,j) = 0.;// zero out constraint values to be ignored
            else D(i,j) = x;
        }
    cout << "D:" << endl;
    cout << D << endl;
    
    // DV is direct transformation from coordinates to track parameters
    
    MatrixXd DV(nTrackParameters,nVars);
    DV = D*V;
    
    cout << "DV:" << endl;
    cout << DV << endl;
    
    
    Statistics statPar[nTrackParameters];
    Statistics statConstr[nVars];
    
    // define region in parameter space
    
    Geometry g;
    
    g.t_phi = 0.;
    g.t_eta = 0.5;
    g.t_z0 = 0.;
    g.t_invPt = 1./2.;
    
    g.t_deltaEta = 0.5;
    g.t_deltaPhi = 2.*g.Pi/8.;
    g.t_deltaZ0 = 0.1;
    
    
    // begin first loop on tracks
    
    cout << endl << "begin first loop on tracks" << endl;
    
    for(long int iT = 1; iT != nTracks+1; ++ iT){// begin loop on tracks
        
        if(iT%100000 == 0) cout << iT << endl;
        
        Track t(g); // construct random track
        
        // find intersections with errors
        
        double x_, z_ ;
        VectorXd vars(nVars);
        vars = VectorXd::Constant(nVars,0.);
        int iV = 0;
        
        for(int iL = 0; iL != g.nLayers; ++ iL){
            if(t.getHit(g, iL, x_, z_)) {
                vars(iV++) = x_;
                vars(iV++) = z_;
            }
        }
        
        // skip track if not all layers are hit        
        if(iV != nVars) continue;
        
        
        // get principal components
        
        VectorXd principal(nVars);
        principal = V*vars;        
               
        VectorXd par(nTrackParameters);
        par = DV*vars;
        
        // parameter errors
        
        VectorXd errPar(nTrackParameters);
        errPar(0) = par[0] - t.phi;
        errPar(1) = par[1] - t.cotTheta;
        errPar(2) = par[2] - t.z0;
        errPar(3) = par[3] - t.invPt;

        // accumulate statistics
        
        for(int i = 0; i != nTrackParameters; ++i) statPar[i].fill(errPar(i));
        for(int i = 0; i != nVars; ++i) statConstr[i].fill(principal(i));
        
        
    } // end first loop on tracks
    
    // print statistics
    
    string parName[nTrackParameters] = {"phi     ", "cotTheta", "z0      ", "invPt   "};

    cout << endl;
    
    for(int i = 0; i != nTrackParameters; ++i) {
        cout << "par "<< i << ": " << parName[i] << ": " << statPar[i].getEntries() << " " << setw(10) << statPar[i].getMean() << " " << setw(10) << statPar[i].getSigma() << endl;
    }
    
    cout << endl;
    
    for(int i = 0; i != nVars; ++i) {
        cout << "constr "<< setw(2) << i  << ": " << statConstr[i].getEntries() << " " << setw(10)<< statConstr[i].getMean() << " " << setw(10) << statConstr[i].getSigma() << endl;
    }
    
    
    //********************************************************************
    //
    // book histograms
    
    TFile* histFile = new TFile((pathToMyDirectory + "regionHists.root").c_str(),"RECREATE");
    
    // coordinates
    
    TH1D hVar[nVars], hPC[nVars];
    for(int iV = 0; iV != nVars; ++ iV){
        ostringstream s;
        s << "Var_" << iV;
        hVar[iV] = TH1D((s.str()).c_str(), (s.str()).c_str(), 100, -2., +2.);
    }
    
    // principal components
    
    double nSigmas = 5.0;    
    for(int iV = 0; iV != nVars; ++ iV){
        ostringstream t;
        t << "PC_" << iV;
        hPC[iV] = TH1D((t.str()).c_str(), (t.str()).c_str(), 100,
                       statConstr[iV].getMean() - nSigmas*statConstr[iV].getSigma(),
                       statConstr[iV].getMean() + nSigmas*statConstr[iV].getSigma());
    }
    
    // track parameters
 
    TH1D hPhi("Phi","Phi",100, -1., +1.);
    TH1D hCotTheta("CotTheta","CotTheta",100, -1., +4.);
    TH1D hZ0("Z0","Z0",100, -0.2, +0.2);
    TH1D hInvPt("InvPt","InvPt",100, -0.6, +0.6);
    
    // errors on parameters
    
    TH1D hErrPhi("ErrPhi","ErrPhi",100,
                 statPar[0].getMean() - nSigmas*statPar[0].getSigma(),
                 statPar[0].getMean() + nSigmas*statPar[0].getSigma());
    TH1D hErrEta("ErrCotTheta","ErrCotTheta",100,
                 statPar[1].getMean() - nSigmas*statPar[1].getSigma(),
                 statPar[1].getMean() + nSigmas*statPar[1].getSigma());
    TH1D hErrZ0("ErrZ0","ErrZ0",100,
                statPar[2].getMean() - nSigmas*statPar[2].getSigma(),
                statPar[2].getMean() + nSigmas*statPar[2].getSigma());
    TH1D hErrInvPt("ErrInvPt","ErrInvPt",100,
                   statPar[3].getMean() - nSigmas*statPar[3].getSigma(),
                   statPar[3].getMean() + nSigmas*statPar[3].getSigma());
    
    
    //********************************************************************
    


    
    // begin second loop on tracks
    
    cout << endl << "begin second loop on tracks" << endl;
    
    for(long int iT = 1; iT != nTracks+1; ++ iT){// begin loop on tracks
        
        if(iT%100000 == 0) cout << iT << endl;
        
        Track t(g); // construct random track
        
        // find intersections with errors
        
        double x_, z_ ;
        VectorXd vars(nVars);
        vars = VectorXd::Constant(nVars,0.);
        int iV = 0;
        
        for(int iL = 0; iL != g.nLayers; ++ iL){
            if(t.getHit(g, iL, x_, z_)) {
                vars(iV++) = x_;
                vars(iV++) = z_;
            }
        }
        
        // skip track if not all layers are hit
        if(iV != nVars) continue;
        
        
        // get principal components
        
        VectorXd principal(nVars);
        principal = V*vars;
        
        VectorXd par(nTrackParameters);
        par = DV*vars;
        
        // parameter errors
        
        VectorXd errPar(nTrackParameters);
        errPar(0) = par[0] - t.phi;
        errPar(1) = par[1] - t.cotTheta;
        errPar(2) = par[2] - t.z0;
        errPar(3) = par[3] - t.invPt;
         
        // fill histograms
        
        for(int iV = 0; iV != nVars; ++ iV){
            hVar[iV].Fill(vars(iV));
            hPC[iV].Fill(principal(iV));
        }
        
        hPhi.Fill(par(0));
        hCotTheta.Fill(par(1));
        hZ0.Fill(par(2));
        hInvPt.Fill(par(3));
        
        hErrPhi.Fill(errPar(0));
        hErrEta.Fill(errPar(1));
        hErrZ0.Fill(errPar(2));
        hErrInvPt.Fill(errPar(3));
        
        
        
    } // end second loop on tracks

    
    histFile->Write();
    
    cout << endl << "************* Done *************" << endl;
    
    return 0;
};

