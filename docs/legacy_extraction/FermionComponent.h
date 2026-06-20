//
//  FermionComponent.h
//  rpa-plasma
//
//  Created by In-Gee Kim on 6/12/14.
//  Copyright (c) 2014 In-Gee Kim. All rights reserved.
//

#ifndef __rpa_plasma__FermionComponent__
#define __rpa_plasma__FermionComponent__

// system header(s)
#include <cmath>
using std::exp;
using std::log;
using std::pow;
using std::sqrt;

#include <complex> 
using std::complex;

// user defined function
double step(double);

// user defined header(s)
#include "Component.h"

class FermionComponent : public Component {
public:
    // default constructor
    FermionComponent();
    // default destructor
    ~FermionComponent();
    // preferred constructor
    FermionComponent(const char, const double, const double, const double);
    
    // zero-temperature Lindhard function (1954)
    const complex<double> complexLindhard(const double, const double) const; // return complex Lindhard function
    const double get_imagLindhard(const double, const double) const;
    const double get_realLindhard(const double, const double) const; 
    
    // extracts the private data
    const double get_FermiWaveNumber(void) const { return _FermiWavevector; };
    const double get_FermiEnergy(void) const { return _FermiEnergy; };
    const double get_DOSatFermiEnergy(void) const { return _DOSatFermiEnergy; };
    
protected:
    double _FermiWavevector;  // (9 pi/4)^(1/3) (1/r_s) in units of a_0^{-1}
    double _FermiEnergy; // (\hbar^2 k_F^2)/(2 m) in units of Hartree
    double _DOSatFermiEnergy;  // density of states at Fermi Energy D(E_F) in units of states/Hartree
};

#endif /* defined(__rpa_plasma__FermionComponent__) */
