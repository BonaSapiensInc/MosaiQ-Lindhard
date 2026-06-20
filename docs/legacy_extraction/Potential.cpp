//
//  Potential.cpp
//  rpa-plasma
//
//  Created by In-Gee Kim on 7/9/14.
//  Copyright (c) 2014 In-Gee Kim. All rights reserved.
//

// The input q value is rationalized in Hatree atomic unit (a_0^-1).
// The output potential is given in units of E_H / e.

#include "Potential.h"


Potential::Potential(void)
: _Ehe(-0.90767486e-1)
{
    // default constructor
}

Potential::~Potential(void)
{
    // default destructor
}

double Potential::constant(double q, vector<double> &parameterList)
{
    const double constant = parameterList[0];
    const double chargeMultiple = parameterList[1];
    
    double potential = constant * chargeMultiple;
    
    return potential;
}

double Potential::yukawa(double q, vector<double> &parameterList)
{
    double cutoff = parameterList[0];
    double chargeMultiple = parameterList[1];
    
    const double fourPi = 4.0 * ( 4.0 * atan( 1.0 ) );
    
    double qSquare = q * q;
    double cutoffSquare = cutoff * cutoff;
    double fourier = fourPi;
    double potential = ( ( fourier * chargeMultiple ) / ( qSquare + cutoffSquare ) );

    return potential;
}

double Potential::lennard_jones(double q, vector<double> &parameterList)
{
    // be sure epsilon and sigma are given in atomic units
    double epsilon = parameterList[0];
    double sigma = parameterList[1];

    double sgn = 1.0;
    if ( q < 0 ) {
        sgn = -1.0;
    }
    const double pi = 4.0 * atan(1.0);
    double rootHalfPi = sqrt( pi / 2.0 );
    
    // Fourier transformed Lennard-Jones potential
    double potential = rootHalfPi*sigma * sgn;
    potential *= ( 332640.0 + pow(q, 6.0) );
    potential *= pow(q, 5.0);
    potential *= epsilon;
    potential /= 39916800.0;
    
    return potential;
}

double Potential::modified_coulomb(double q, vector<double> &parameterList)
{
    double cutoff = parameterList[0];
    double chargeMultiple = parameterList[1];
    
    const double fourPi = 4.0 * ( 4.0 * atan( 1.0 ) );
    
    double qSquare = q * q;
    double qQuad = qSquare * qSquare;
    double cutoffSquare = cutoff * cutoff;
    double fourier = fourPi;
    double potential = ( ( fourier * chargeMultiple ) / ( qSquare + cutoffSquare * qQuad ) );
    
    return potential;
}