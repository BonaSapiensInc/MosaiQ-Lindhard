//
//  Plasma.h
//  RPA-Plasma
//
//  Created by In-Gee Kim on 12/2/14.
//  Copyright (c) 2014 In-Gee Kim. All rights reserved.
//

#ifndef __RPA_Plasma__Plasma__
#define __RPA_Plasma__Plasma__

// system headers
#include <cmath>
using std::atan;
using std::tanh;
using std::sqrt;

#include <complex>
using std::complex;
using std::conj;

#include <fstream>
using std::ofstream;

#include <iostream>
using std::endl;

#include <string>
using std::string;

#include <vector>
using std::vector;

// user defined headers
#include "BosonComponent.h"
#include "TemperatureFermionComponent.h"
#include "MaxwellianComponent.h"
#include "Potential.h"
#include "Plot.h"

class Plasma {
public:
    Plasma();  // default constructor
    ~Plasma(); // default destructor
    
    // set up the plasma state
    void add_component(const BosonComponent *);
    void add_component(const TemperatureFermionComponent *);
    void add_component(const MaxwellianComponent *);
    
    void set_eePotentialMode(const long, const vector<double> &);
    void set_eiPotentialMode(const long, const vector<double> &);
    void set_iiPotentialMode(const long, const vector<double> &);
    void set_PotentialvTable(void);
    
    void set_weightedTemperatures(void);

    // get the primitive information
    long get_nComponent(void) const;
    long get_componentType(const long) const;
    double get_charge(const long) const;
    double get_density(const long) const;
    double get_plasmaFrequency(const long) const;
    double get_DebyeWaveNumber(const long) const;
    double get_intrinsicFrequency(const long) const;
    double get_intrinsicWaveNumber(const long) const;
    double get_intrinsicTemperature(const long) const;
    string get_designationString(const long) const;
    complex<double> get_chi0(const long, const double, const double, const double) const;
    complex<double> get_chi0(const long, const double, const double) const;
    complex<double> get_dielectric(const long, const long, const long, const double, const double) const;
    complex<double> get_T0dielectric(const long, const long, const long, const double, const double) const;
    complex<double> get_chiRPA(const long, const long, const long, const double, const double) const;
    complex<double> get_T0chiRPA(const long, const long, const long, const double, const double) const;
    double get_Skw(const long, const long, const long, const double, const double) const;
    double get_Sq(const long, const long, const long, const double) const;
    double get_Fq(const long, const long, const long, const double) const;
    
    // object actions
    void print_Status(ofstream &) const;
    void print_chi0(const Plot &, ofstream &) const;
    void print_dielectric(const Plot &, ofstream &) const;
    void print_chiRPA(const Plot &, ofstream &) const;
    void print_Skw(const Plot &, ofstream &) const;
    void print_sumRules(const Plot &, ofstream &) const;
    
private:
    // Temperature Fermion Components
    vector<TemperatureFermionComponent> _fermions;
    
    // Thermal Maxwellican Components
    vector<MaxwellianComponent> _maxwells;
    
    // Thermal Boson Components
    vector<BosonComponent> _bosons;
    
    // Component indexing
    vector<long> _componentType;  // the type (boson -1, maxwell 0, fermion +1) of the s-th component
    vector<long> _componentIndex; // the index of the component type of _componentType[s]
    
    // Potential information
    long _eePotentialMode;
    long _eiPotentialMode;
    long _iiPotentialMode;
    
    vector<double> _eePotentialParameter;
    vector<double> _eiPotentialParameter;
    vector<double> _iiPotentialParameter;
    
    double ( Potential::*tblPotential[4] ) (double, vector<double> &); // potential functor
    
    // weighted temperature
    vector<double> _weightedTemperature;
    const double _constantBoltzmann;
};

#endif /* defined(__RPA_Plasma__Plasma__) */
