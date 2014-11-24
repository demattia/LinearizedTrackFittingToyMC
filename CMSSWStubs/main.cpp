#include "L1TrackTriggerTree.h"
#include "Eigen/Eigenvalues"
#include <unordered_set>
#include <sstream>
#include <math.h>
#include "TFile.h"
#include "TH1D.h"

using namespace Eigen;

int fillVars(const L1TrackTriggerTree * tree, VectorXd & vars)
{
  std::unordered_set<int> layersFound;
  int iV = 0;
  for (int k=0; k<tree->m_stub; ++k) {
    // Use only stubs from muons (skip also antimuons)
    if (tree->m_stub_pdg->at(k) == 13) {
      // Need to skip the opposite side otherwise it will cause a discontinuity
      if (tree->m_stub_etaGEN->at(k) < 0) continue;
      // std::cout << "eta gen = " << tree->m_stub_etaGEN->at(k) << std::endl;
      // float phiGen = atan2(tree->m_stub_pyGEN->at(k), tree->m_stub_pxGEN->at(k));
      // std::cout << "phi gen = " << phiGen << std::endl;
      // Avoid duplicates
      // std::cout << "all layers("<<k<<") = " << tree->m_stub_layer->at(k) << std::endl;
      if (layersFound.insert(tree->m_stub_layer->at(k)).second) {
        // std::cout << "layer("<<k<<") = " << tree->m_stub_layer->at(k) << std::endl;
        float phi = atan2(tree->m_stub_y->at(k), tree->m_stub_x->at(k));
        // phi = (phi>=0) ? phi : phi+2*M_PI;
        // std::cout << "before phi = " << phi << std::endl;
        // phi < 0. ? phi = M_PI/2. + (M_PI + phi) : phi = M_PI/2. - phi;
        // std::cout << "after phi = " << phi << std::endl;
        vars(iV++) = phi;
        vars(iV++) = tree->m_stub_z->at(k);
      }
    }
  }
  return layersFound.size();
}

int main()
{
  int nLayers = 6;
  int nVars = nLayers*2;

  L1TrackTriggerTree * tree = new L1TrackTriggerTree("extracted.root");


  TFile* histFile = new TFile("matrixHists.root", "RECREATE");

  TH1D * hVar[nVars], * hPC[nVars], * hPCNorm[nVars];
  float xRange = 0.001;
  float varRange = 1.;
  for(int iV = 0; iV != nVars; ++iV){
    std::ostringstream s;
    std::ostringstream t;
    std::ostringstream tNorm;
    s << "Var_" << iV;
    iV%2 != 0 ? varRange = 100. : varRange = 1.;
    hVar[iV] = new TH1D((s.str()).c_str(), (s.str()).c_str(), 100, -varRange, +varRange);
    t << "PC_" << iV;
    tNorm << "PCNorm_" << iV;
    if (iV > 8) xRange = 0.1;
    else if (iV > 5) xRange = 0.01;
    hPC[iV] = new TH1D((t.str()).c_str(), (t.str()).c_str(), 100, -xRange, +xRange);
    hPCNorm[iV] = new TH1D((tNorm.str()).c_str(), (tNorm.str()).c_str(), 100, -10., +10.);
  }

//  TH1D hCotTheta("CotTheta","CotTheta",100, -1., +1.);
//  TH1D hPhi("Phi","Phi",100, -1., +1.);
//  TH1D hZ0("Z0","Z0",100, -0.2, +0.2);
//  TH1D hInvPt("InvPt","InvPt",100, -0.6, +0.6);

//  TH1D hErrEta("ErrCotTheta","ErrCotTheta",100, -0.1, +0.1);
//  TH1D hErrPhi("ErrPhi","ErrPhi",100, -0.01, +0.01);
//  TH1D hErrZ0("ErrZ0","ErrZ0",100, -0.01, +0.01);
//  TH1D hErrInvPt("ErrInvPt","ErrInvPt",100, -0.01, +0.01);

  TH1D hNormChi2("NormChi2", "NormChi2", 100, 0, 10);
//  TH1D hNormChi2Params("NormChi2Params", "NormChi2Params", 100, 0, 10);
//  TH1D hNormChi2Diff("NormChi2Diff", "NormChi2Diff", 100, -10, 10);


  MatrixXd cov(nVars, nVars);
  cov = MatrixXd::Constant(nVars, nVars, 0.);
  VectorXd meanValues(nVars);
  meanValues = VectorXd::Constant(nVars, 0.);

  int nTracks = 0;

  for (int i=0; i<tree->n_entries/2; ++i) {
  // for (int i=0; i<10; ++i) {
    tree->getEntry(i);
    // tree->printInfo();

    if (tree->m_stub < 6) continue;

    VectorXd vars(nVars);
    vars = VectorXd::Constant(nVars,0.);
    int layersFound = fillVars(tree, vars);

    if (layersFound == nLayers) {
      ++nTracks;
      for (int iVar=0; iVar<nVars; ++iVar) {
        // update mean
        meanValues(iVar) += (vars[iVar] - meanValues(iVar))/nTracks;

        // update covariance matrix
        if(nTracks == 1) continue; // skip first track
        for (int jVar=0; jVar<nVars; ++jVar) {
          cov(iVar, jVar) += (vars[iVar] - meanValues(iVar))*(vars[jVar] - meanValues(jVar))/(nTracks-1) - cov(iVar, jVar)/nTracks;
        }
      }
    }
  }

  // Diagonalize matrix to find principal components
  SelfAdjointEigenSolver<MatrixXd> es(cov);
  std::cout << "Sqrt(eigenvalues) of cov:" << std::endl;
  std::vector<float> sqrtEigenvalues;
  for(int i = 0; i != nVars; ++i) {
    sqrtEigenvalues.push_back(sqrt(es.eigenvalues()[i]));
    std::cout << " " << sqrt(es.eigenvalues()[i]);
  }
  std::cout << std::endl;

  // V in the ortogonal transformation from variable space to parameter space
  // Parameters are constraints + rotated track parameters

  MatrixXd V = (es.eigenvectors()).transpose();


  // Second loop on tracks

  for (int i=tree->n_entries/2+1; i<tree->n_entries; ++i) {
    tree->getEntry(i);

    VectorXd vars(nVars);
    vars = VectorXd::Constant(nVars, 0.);
    int layersFound = fillVars(tree, vars);
    if (layersFound != nLayers) continue;

    // get principal components
    VectorXd principal(nVars);
    principal = V*(vars - meanValues);

    float chi2 = 0.;
    int nDof = 0;
    // fill histograms
    for(int iV = 0; iV != nVars; ++iV){
      hVar[iV]->Fill(vars(iV));
      // std::cout << "vars("<<iV<<") = " << vars(iV) << std::endl;
      hPC[iV]->Fill(principal(iV));
      hPCNorm[iV]->Fill(principal(iV)/sqrtEigenvalues[iV]);
      if (iV > 8) continue;
      chi2 += (principal(iV)/sqrtEigenvalues[iV])*(principal(iV)/sqrtEigenvalues[iV]);
      ++nDof;
    }
    hNormChi2.Fill(chi2/nDof);

  }
  histFile->Write();
}
