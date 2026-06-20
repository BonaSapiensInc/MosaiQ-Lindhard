//
//  TemperatureFermionComponent.h
//  rpa-plasma
//
//  Created by In-Gee Kim on 6/12/14.
//  Copyright (c) 2014 In-Gee Kim. All rights reserved.
//

#ifndef __rpa_plasma__TemperatureFermionComponent__
#define __rpa_plasma__TemperatureFermionComponent__

// system header(s)
#include <iostream>
using std::cerr;
using std::endl;
using std::string;

#include <fstream>
using std::ofstream;

#include <iomanip>
using std::fixed;
using std::scientific;
using std::setprecision;
using std::setw;

#include <cmath>
using std::atan;
using std::cos;
using std::exp;
using std::log;
using std::sin;
using std::sqrt;
using std::pow;
using std::signbit;

#include <complex>
using std::arg;
using std::complex;
using std::polar;

// user defined headers and functions
#include "FermionComponent.h"

class TemperatureFermionComponent : public FermionComponent {
public:
    // default constructor
    TemperatureFermionComponent();
    // default destructor
    ~TemperatureFermionComponent();
    // preferred constructors
    TemperatureFermionComponent(const char, const double, const double, const double, const double, const bool);
    
    // extracts the private data
    void printStatus(ofstream &) const;
    const double get_temperature(void) const { return _temperature; };
    const double get_tempHartree(void) const { return _tempHartree; };
    const double get_reducedTemperature(void) const { return _reducedTemperature; };
    const double get_DebyeWaveNumber(void) const { return _DebyeWaveNumber; };
    const double get_chemicalPotential(void) const { return _chemicalPotential; };
    const double get_reducedChemicalPotential(void) const { return _reducedChemicalPotential; };
    const bool get_isMaxwell(void) const { return _isMaxwell; };
    
    // imaginary Lindhard function
    const double imaginaryLindhard(const double, const double, const double) const; // selector
    const double get_FCimagLindhard(const double, const double, const double) const; // Fried and Conte (1961) classical limit
    const double get_KGimaginaryLindhard(double const, double const, double const) const; // Khanna and Glyde (1976)
    const double get_GVimagLindhard(const double, const double, const double) const; // Giuliani and Vignale (2005)

    // several trials for the Lindhard function
    const double realLindhard(double const, double const, double const) const; // selector
    const double get_FCrealLindhard(const double, const double, const double) const; // Fried and Conte (1961) classical limit
    const double get_GDrealLindhard(double const, double const, double const) const; // Gouedard and Deutsch (1978)
    const double get_KLrealLindhard(double const, double const, double const) const; // Kim and Lee (1998)
    const double get_KGrealLindhard(double const, double const, double const) const; // Khanna and Glyde (1976)
    const double get_KKrealLindhard(double const, double const, double const) const; // Kramers-Kronig via Stenger (1981)
    const double get_GVrealLindhard(double const, double const, double const) const; // Giuliani and Vignale (2005)
    
    // complex and square of Lindhard functions using the above functions.
    const complex<double> complexLindhard(double const, double const, double const) const;
    const double squareLindhard(double const, double const, double const) const;
    
    // setting up function(s)
    void set_ChemicalPotential( void );
    
    // intrinsic function
    const double theFermiDiracIntegral(long, long, double) const;
    
private:
    const double _constantBoltzmann; // Boltzmann constant in units of Hartree/K

    double _temperature;  // temperature in units of K
    double _tempHartree; // temperature in units of Hartree
    double _reducedTemperature;  // \tau in units of k_B T / E_F
    double _chemicalPotential;  // in units of Hartree
    double _reducedChemicalPotential; // \gamma in units of Fermi Energy
    double _DebyeWaveNumber; // k_D = sqrt( 3 Z^2 (E_H / k_B T) r_s^3) in units of 1/a_0
    bool _isMaxwell; // choose the Maxwellian distribution function
    
    // plasma dispersion function Phi by Fried and Conte (1961)
    const double approxPhi(double) const;
    const double plasmaDispersion(double) const;

    // integrals for real part of the Lindhard function by Kim and Lee (1998)
    const double KLresidueSum(long const, double const, double const, double const) const;
    const double KLsumEulerMaclaurin(long const, double const, double const, double const) const;
};

#endif /* defined(__rpa_plasma__TemperatureFermionComponent__) */
