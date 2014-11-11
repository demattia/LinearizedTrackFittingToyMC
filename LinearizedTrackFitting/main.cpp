//
//  main.cpp
//  LinearizedTrackFitting
//
//  Created by Luciano Ristori on 4/25/14.
//
//  Revised on 5/19/14
//  Revised on 5/23/14

#include <iostream>
#include <string>
#include <random>
#include "TFile.h"

#include "makeMatrix.h"
#include "testRegion.h"

using namespace std;


// random generator

 std::default_random_engine generator;
 std::uniform_real_distribution<double> distribution(0.,1.);
 std::normal_distribution<double> gauss(0.0,1.0);


string pathToMyDirectory = "/Users/demattia/Desktop/LinearizedTrackFitting/";
// string pathToMyDirectory = "/Users/Luciano/Documents/Ccode/LinearizedTrackFitting/";

TFile* histFile ;

int main(int argc, const char * argv[]) {
    
    // generate tracks, make principal components, (write matrix to file - TBI), make some histograms    
    
    // first parameter = number of tracks for first loop
    // second parameter = number of tracks for second loop
    // third parameter = number of tracks for third loop
    
    makeMatrix(100e6,100e6,100e3);
    
    // test region of linearity
    
    testRegion(100e3);
    
    
    
    return 0;
}

