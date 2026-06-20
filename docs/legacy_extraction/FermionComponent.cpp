//
//  FermionComponent.cpp
//  rpa-plasma
//
//  Created by In-Gee Kim on 6/12/14.
//  Copyright (c) 2014 In-Gee Kim. All rights reserved.
//

#include "FermionComponent.h"

FermionComponent::FermionComponent()
{
    // default constructor is left as an empty-body
}

FermionComponent::~FermionComponent()
{
    // default destructor is left as an empty-body
}

FermionComponent::FermionComponent(const char name, const double mass, const double charge, const double radiusWS)
: Component(name, mass, charge, radiusWS)
{
    const double pi = _pi;
    const double oneThird = 1.0 / 3.0;
    const double ninePiFour = pow(9.0 * pi / 4.0, oneThird);
    
    double r_s = _radiusWignerSeitz;
    double mu = _mass; // mass ratio w.r.t. m_e
    
    _FermiWavevector = ( ninePiFour / r_s );
    _FermiEnergy = ( 1.0 / ( 2.0 * mu ) ) * (ninePiFour * ninePiFour ) * ( 1.0 / (r_s * r_s) );
    _DOSatFermiEnergy = ( 3.0 * _density ) / ( 2.0 * _FermiEnergy );
}

const complex<double> FermionComponent::complexLindhard(const double q, const double omega) const
{
    double real = get_realLindhard(q, omega);
    double imag = get_imagLindhard(q, omega);
    // const double DeF = _DOSatFermiEnergy;
    // real *= DeF;
    // imag *= DeF;
    
    complex<double> Lindhard(real, imag);
    
    return Lindhard;
}

const double FermionComponent::get_realLindhard(const double q, const double w) const
{
    // zero-temperature Lindhard function / DeF
    // J. Lindhard, K. Dan. Vidensk. Selsk. Mat.-Fys. Medd. 28 (1954) 8.
    // implemented based on Eq. (4.24) of Giuliani and Vignale (2005).
    
    double nu_p = ( w/q + q )/2.0;
    double nu2_p = nu_p * nu_p;
    double argNum_p = nu_p + 1.0;
    double argDen_p = nu_p - 1.0;
    
    double nu_n = ( w/q - q )/2.0;
    double nu2_n = nu_n * nu_n;
    double argNum_n = nu_n + 1.0;
    double argDen_n = nu_n - 1.0;
    
    double arg_p = sqrt( (argNum_p / argDen_p) * (argNum_p / argDen_p) );
    double term_p = log(arg_p);
    term_p *= ((1.0 - nu2_p) / (2.0 * q));
    
    double arg_n = sqrt( (argNum_n / argDen_n) * (argNum_n / argDen_n) );
    double term_n = log(arg_n);
    term_n *= ((1.0 - nu2_n) / (2.0 * q));
    
    double real = -(1.0 - term_n + term_p);
    
    return real;
}

const double FermionComponent::get_imagLindhard(const double q, const double w) const
{
    const double pi = _pi;
    // zero temperature Lindhard function / DeF
    // J. Lindhard, K. Dan. Vidensk, Selsk. Mat.-Fys. Medd. 28 (1954) 8.
    // implementation based on Eq. (4.25) of Giuliani and Vignale (2005).
    
    double nu_p = ( w/q + q ) / 2.0;
    double nu2_p = nu_p * nu_p;
    
    double nu_n = ( w/q - q ) / 2.0;
    double nu2_n = nu_n * nu_n;
    
    double xn = 1.0 - nu2_n;
    double xp = 1.0 - nu2_p;
    
    double imag = step(xn) * xn - step(xp) * xp;
    imag *= -( pi / (2.0 * q) );
    
    return imag;
}

