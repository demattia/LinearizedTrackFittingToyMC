//
//  makeMatrix.cpp
//  LinearizedTrackFitting
//
//  Created by Luciano Ristori on 5/23/14.
//  Copyright (c) 2014 Luciano Ristori. All rights reserved.
//

#include <iostream>
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

int makeMatrix(long int maxTracks1, long int maxTracks2, long int maxTracks3) {
    
    // maxTracks1 = number of tracks for first loop
    // maxTracks2 = number of tracks for second loop
    // maxTracks3 = number of tracks for third loop

    
 
    
    bool verbose = false;
    
    // instantiate geometry
    
    Geometry g;

    // number of hit coordinates for a track
    
    const int nVars = 12;
    
    // number of track parameters
    // here we consider 4 degrees of freedom for an
    // helix intersecting the beam line (x=y=0)
    
    const int nTrackParameters = 4;
    
    
    
    //********************************************************************
    //
    // book histograms
    
    TFile* histFile = new TFile((pathToMyDirectory + "matrixHists.root").c_str(),"RECREATE");
    
    TH1D hVar[nVars], hPC[nVars];
    for(int iV = 0; iV != nVars; ++ iV){
        ostringstream s;
        s << "Var_" << iV;
        hVar[iV] = TH1D((s.str()).c_str(), (s.str()).c_str(), 100, -1., +1.);
    }
    
    for(int iV = 0; iV != nVars; ++ iV){
        ostringstream t;
        t << "PC_" << iV;
        hPC[iV] = TH1D((t.str()).c_str(), (t.str()).c_str(), 100, -0.001, +0.001);
    }
    
    TH1D hCotTheta("CotTheta","CotTheta",100, -1., +1.);
    TH1D hPhi("Phi","Phi",100, -1., +1.);
    TH1D hZ0("Z0","Z0",100, -0.2, +0.2);
    TH1D hInvPt("InvPt","InvPt",100, -0.6, +0.6);
    
    TH1D hErrEta("ErrCotTheta","ErrCotTheta",100, -0.01, +0.01);
    TH1D hErrPhi("ErrPhi","ErrPhi",100, -0.01, +0.01);
    TH1D hErrZ0("ErrZ0","ErrZ0",100, -0.005, +0.005);
    TH1D hErrInvPt("ErrInvPt","ErrInvPt",100, -0.005, +0.005);
    
    
    //********************************************************************
    
    
    
    
    
    // initialize covariance matrix and mean vector
    
    MatrixXd cov(nVars, nVars);
    cov = MatrixXd::Constant(nVars, nVars, 0.);
    VectorXd meanValues(nVars);
    meanValues = VectorXd::Constant(nVars, 0.);
    
    
    // begin first loop on tracks
    
    cout << "begin first loop on tracks" << endl;
    
    long int nTracks1 = 0;
    
    for(long int iT = 1; iT != maxTracks1+1; ++ iT){
        
        if(iT%100000 == 0) cout << iT << endl;
        
        Track t; // construct random track
        
        std::vector<double> x,z, vars;
        
        double x_, z_ ;
        
        for(int iL = 0; iL != g.nLayers; ++ iL){
                if(t.getHit(g, iL, x_, z_)) {
                    x.push_back(x_);
                    vars.push_back(x_);
                    z.push_back(z_);
                    vars.push_back(z_);
                }
            
        }
        
        // skip track if not all layers are hit
        
        if(x.size() != g.nBarrels){
            cout << "*** track " << iT << " x: " << x.size() << " layers instead of " << g.nBarrels << endl;
            continue;
        }
        
        if(z.size() != g.nBarrels){
            cout << "*** track " << iT << " z: " << z.size() << " layers instead of " << g.nBarrels << endl;
            continue;
        }
        
        
        if(verbose){ // dumps track data
            cout << endl;
            cout << "1st round track " << iT << " primary vertex: " << t.x0 << " " << t.y0 << " " << t.z0 << endl;
            cout << "invPt: " << t.invPt << " eta: " << t.eta << " phi: " << t.phi << endl;
            cout << "x: ";
            for(int iL = 0; iL != g.nBarrels; ++ iL)cout << " "<< x[iL];
            cout << endl;
            
            cout << "z: ";
            for(int iL = 0; iL != g.nBarrels; ++ iL)cout << " "<< z[iL];
            cout << endl;
        }
        
        // update means
        
        ++nTracks1;
        
        for(int iVar = 0; iVar != nVars; ++iVar){
            meanValues(iVar) += (vars[iVar] - meanValues(iVar))/nTracks1;
        };
        
        // update covariance matrix
        
        if(iT == 1) continue; // skip first track
        
        for(int i = 0; i != nVars; ++i)
            for(int j = 0; j != nVars; ++j){
                cov(i, j) += (vars[i] - meanValues(i))*(vars[j] - meanValues(j))/(nTracks1-1) - cov(i, j)/nTracks1;
            }
        
        
    } // end first loop on tracks
    
    if(verbose){ // dumps covariance matrix and mean values
        cout << meanValues << endl << endl;
        cout << cov << endl;
    }
    
    
    
    // find eigenvectors of covariance matrix
    
    SelfAdjointEigenSolver<MatrixXd> es(cov);
    cout << "Sqrt(eigenvalues) of cov:" << endl;
    for(int i = 0; i != nVars; ++i) cout << " " << sqrt(es.eigenvalues()[i]);
    cout << endl;
    
    // V in the ortogonal transformation from variable space to parameter space
    // Parameters are constraints + rotated track parameters
    
    MatrixXd V = (es.eigenvectors()).transpose();
    
    
    // prepare for second loop on tracks
    
    
    // correlation between diagonalized coordinates (V) and parameters (P)
    
    MatrixXd corrPV(nTrackParameters, nVars);
    corrPV = MatrixXd::Constant(nTrackParameters, nVars, 0.);
    VectorXd meanP(nTrackParameters);
    meanP = VectorXd::Constant(nTrackParameters, 0.);
    VectorXd meanV(nVars);
    meanV = VectorXd::Constant(nVars, 0.);
    
    long int nTracks2 = 0; // number of tracks actually processed
    
    
    cout << "begin second loop on tracks" << endl;
    
    for(long int iT = 1; iT != maxTracks2+1; ++ iT){// begin second loop on tracks
        
        if(iT%100000 == 0) cout << iT << endl;
        
        Track t; // construct random track
        
        std::vector<double> x,z;
        
        double x_, z_ ;
        VectorXd vars(nVars);
        vars = VectorXd::Constant(nVars,0.);
        
        // find intersections with errors
        
        int iV =0;
        for(int iL = 0; iL != g.nLayers; ++ iL){
            if(t.getHit(g, iL, x_, z_)) {
                x.push_back(x_);
                vars(iV++) = x_;
                z.push_back(z_);
                vars(iV++) = z_;
            }
        }
        if(x.size() != g.nBarrels){
            cout << "*** track " << iT << " x: " << x.size() << " layers instead of " << g.nBarrels << endl;
            continue;
        }
        
        if(z.size() != g.nBarrels){
            cout << "*** track " << iT << " z: " << z.size() << " layers instead of " << g.nBarrels << endl;
            continue;
        }
        
        if(verbose){ // dumps track data
            cout << endl;
            cout << "2st round track " << iT << " primary vertex: " << t.x0 << " " << t.y0 << " " << t.z0 << endl;
            cout << "invPt: " << t.invPt << " eta: " << t.eta << " phi: " << t.phi << endl;
            cout << "x: ";
            for(int iL = 0; iL != g.nBarrels; ++ iL)cout << " "<< x[iL];
            cout << endl;
            
            cout << "z: ";
            for(int iL = 0; iL != g.nBarrels; ++ iL)cout << " "<< z[iL];
            cout << endl;
        }
        
        // transform coordinates to principal components
        
        VectorXd principal(nVars);
        principal = V*(vars - meanValues);
        
        // create vector of track parameters
        
        VectorXd par(nTrackParameters);
        par(0) = t.phi;
        par(1) = t.cotTheta;
        par(2) = t.z0;
        par(3) = t.invPt;
        
        
        // update correlation matrix between parameters and principal components
        
        ++nTracks2;
        
        // update means
        
        for(int iVar = 0; iVar != nVars; ++iVar){
            meanV(iVar) += (principal[iVar] - meanV(iVar))/nTracks2;
        };
        
        for(int iPar = 0; iPar != nTrackParameters; ++iPar){
            meanP(iPar) += (par[iPar] - meanP(iPar))/nTracks2;
        };
        
        // update covariance matrix
        
        if(iT == 1) continue; // skip first track
        
        for(int i = 0; i != nTrackParameters; ++i)
            for(int j = 0; j != nVars; ++j){
                corrPV(i, j) += (par[i] - meanP(i))*(principal[j] - meanV(j))/(nTracks2-1) - corrPV(i, j)/nTracks2;
            }
        
    } // end second loop on tracks
    
    
    // invert (diagonal)  correlation matrix dividing by eigenvalues
    // *** check the math ***
    
    MatrixXd D(nTrackParameters,nVars); // transformation from coordinates to track parameters
    
    for(int iP = 0; iP != nTrackParameters; ++iP)
        for(int iV = 0; iV != nVars; ++iV)
            D(iP, iV) = corrPV(iP, iV)/es.eigenvalues()[iV];
    
    
    cout << endl;
    cout << "V:" << endl;
    cout << V << endl;
    cout << "D:" << endl;
    cout << D << endl;
    
    
    // open matrix file and write V and D arrays
    
    cout<< "opening matrixVD.txt for writing" << endl;
    
    std::ofstream outfile;
    outfile.open(pathToMyDirectory + "matrixVD.txt");
    if(!outfile) {
        cout << "error opening matrixVD.txt" << endl;
        return -1;
    }
    outfile << V;
    outfile << endl << endl;
    outfile << D;
    outfile << endl;
    
    
    // begin third loop on tracks
    
    //Statistics s[nVars];
    
    cout << "begin third loop on tracks" << endl;
    
    for(long int iT = 1; iT != maxTracks3+1; ++ iT){// begin third loop on tracks
        
        if(iT%100000 == 0) cout << iT << endl;
        
        Track t;
        
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
        
        // skip track if not all layers are hit !!!! ***************************
        
        
        // get principal components
        
        VectorXd principal(nVars);
        principal = V*(vars - meanValues);
        
        // fill histograms
        
        for(int iV = 0; iV != nVars; ++ iV){
            hVar[iV].Fill(vars(iV));
            hPC[iV].Fill(principal(iV));
        }
        
        // zero a number of components
        
        for(int i = 0; i != 6; ++i) principal(i) = 0.;
        
        // get parameters
        
        VectorXd par(nTrackParameters);
        par = D*principal + meanP;
        
        // parameter errors
        
        VectorXd errPar(nTrackParameters);
        errPar(0) = par[0] - t.phi;
        errPar(1) = par[1] - t.cotTheta;
        errPar(2) = par[2] - t.z0;
        errPar(3) = par[3] - t.invPt;
        
        
        if(verbose){
            cout << "track " << iT << endl;
            cout << "errPar: " << errPar.transpose() << endl;
        }
        
        
        hPhi.Fill(par(0));
        hCotTheta.Fill(par(1));
        hZ0.Fill(par(2));
        hInvPt.Fill(par(3));
        
        hErrPhi.Fill(errPar(0));
        hErrEta.Fill(errPar(1));
        hErrZ0.Fill(errPar(2));
        hErrInvPt.Fill(errPar(3));
        
        
    } // end third loop on tracks
    
    //for(int i = 0; i != nSignificantEigenvectors; ++i)
    //cout << "v(" << i << ") mean: " << s[i].getMean() << " sigma: " << s[i].getSigma() << endl;
    
    
    
    
    
    histFile->Write();
    
    return 0;
}

