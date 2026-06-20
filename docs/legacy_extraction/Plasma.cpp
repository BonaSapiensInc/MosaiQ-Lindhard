//
//  Plasma.cpp
//  RPA-Plasma
//
//  Created by In-Gee Kim on 12/2/14.
//  Copyright (c) 2014 In-Gee Kim. All rights reserved.
//

#include "Plasma.h"

Plasma::Plasma()
: _constantBoltzmann(3.1668114e-6)
{
    // empty body for the default constructor
}

Plasma::~Plasma()
{
    // empty body for the default destructor
}

void Plasma::add_component(const BosonComponent *component)
{
    _bosons.push_back(*component);
    _componentType.push_back(-1);
    long s = _bosons.size(); // s is the current number of elements
    s--; // s becomes the index
    _componentIndex.push_back(s);
}

void Plasma::add_component(const TemperatureFermionComponent *component)
{
    _fermions.push_back(*component);
    _componentType.push_back(+1);
    long s = _fermions.size(); // s is the current number of elements
    s--; // s becomes the index
    _componentIndex.push_back(s);
}

void Plasma::add_component(const MaxwellianComponent *component)
{
    _maxwells.push_back(*component);
    _componentType.push_back(0);
    long s = _maxwells.size(); // s is the current number of elements
    s--;  // s becomes the index
    _componentIndex.push_back(s);
}

void Plasma::set_eePotentialMode(const long potentialMode, const vector<double> &potentialParameter)
{
    _eePotentialMode = potentialMode;
    _eePotentialParameter = potentialParameter;
}

void Plasma::set_eiPotentialMode(const long potentialMode, const vector<double> &potentialParameter)
{
    _eiPotentialMode = potentialMode;
    _eiPotentialParameter = potentialParameter;
}

void Plasma::set_iiPotentialMode(const long potentialMode, const vector<double> &potentialParameter)
{
    _iiPotentialMode = potentialMode;
    _iiPotentialParameter = potentialParameter;
}

void Plasma::set_PotentialvTable(void)
{
    tblPotential[0] = &Potential::constant;
    tblPotential[1] = &Potential::yukawa;
    tblPotential[2] = &Potential::lennard_jones;
    tblPotential[3] = &Potential::modified_coulomb;
}

void Plasma::set_weightedTemperatures(void)
{
    long nComponent = this->get_nComponent();
    
    vector<double> T(nComponent, 0.0);
    vector<double> m(nComponent, 1.0);
    
    long iType;
    long idx;
    for ( long s = 0; s < nComponent; s++ ) {
        iType = _componentType[s];
        idx = _componentIndex[s];
        switch (iType) {
            case (-1):  // boson
                m[s] = _bosons[idx].get_mass();  // in units of m_e
                T[s] = _bosons[idx].get_temperature();  // in units of K
                break;
            case ( 0):  // maxwellian
                m[s] = _maxwells[idx].get_mass();
                T[s] = _maxwells[idx].get_temperature();
                break;
            case (+1):  // fermion
                m[s] = _fermions[idx].get_mass();
                T[s] = _fermions[idx].get_temperature();
                break;
        }
    }
    
    double wt = 0.0;
    for ( long s = 0; s < nComponent; s++ ) {
        for ( long t = 0; t < nComponent; t++ ) {
            wt = ( m[s] * T[t] + m[t] * T[s] ) / ( m[s] + m[t] );
            _weightedTemperature.push_back(wt);
        }
    }
}

void Plasma::print_Status(ofstream &logFile) const
{
    long nComponent = this->get_nComponent();
    
    long iType;
    long iName;
    long idx;
    
    string componentName[3] = {"bosonic", "maxwellian", "fermionic" };
    
    for ( long s = 0; s < nComponent; s++ ) {
        logFile << "<<< The generated ";
        iType = _componentType[s];
        iName = iType + 1;
        idx = _componentIndex[s];
        logFile << componentName[iName] << " component >>>" << endl;
        switch (iType) {
            case (-1): // boson
                _bosons[idx].printStatus(logFile);
                break;
            case ( 0): // maxwellian
                _maxwells[idx].printStatus(logFile);
                break;
            case (+1): // fermion
                _fermions[idx].printStatus(logFile);
                break;
        }
        logFile << endl;
        logFile << "<<< The generated " << componentName[iName]
            << " component >>>" << endl;
        logFile << endl;
    }
}

long Plasma::get_nComponent(void) const
{
    long nComponents = 0;
    nComponents = _componentType.size();
    
    return nComponents;
}

long Plasma::get_componentType(const long s) const
{
    return _componentType[s];
}

double Plasma::get_charge(const long s) const
{
    long iType = _componentType[s];
    long idx = _componentIndex[s];
    double charge = 0.0;
    switch ( iType ) {
        case (-1): // boson
            charge = _bosons[idx].get_charge();
            break;
        case (0): // maxwellian
            charge = _maxwells[idx].get_charge();
            break;
        case (+1): // fermion
            charge = _fermions[idx].get_charge();
            break;
    }
    
    return charge;
}

double Plasma::get_density(const long s) const
{
    long iType = _componentType[s];
    long idx = _componentIndex[s];
    double density = 0.0;
    switch ( iType ) {
        case (-1): // boson
            density = _bosons[idx].get_density();
            break;
        case (0): // maxwellian
            density = _maxwells[idx].get_density();
            break;
        case (+1): // fermion
            density = _fermions[idx].get_density();
            break;
    }
    
    return density;
}

double Plasma::get_plasmaFrequency(const long s) const
{
    long iType = _componentType[s];
    long idx = _componentIndex[s];
    double w_p = 0.0;
    switch ( iType ) {
        case (-1): // boson
            w_p = _bosons[idx].get_plasmaFrequency();
            break;
        case ( 0): // maxwellian
            w_p = _maxwells[idx].get_plasmaFrequency();
            break;
        case (+1): // fermions
            w_p = _fermions[idx].get_plasmaFrequency();
            break;
    }
    
    return w_p; // in energy units
}

double Plasma::get_DebyeWaveNumber(const long s) const
{
    long iType = _componentType[s];
    long idx = _componentIndex[s];
    double k_D = 0.0;
    switch ( iType ) {
        case (-1): // boson
            k_D = _bosons[idx].get_DebyeWaveNumber();
            break;
        case ( 0): // maxwellian
            k_D = _maxwells[idx].get_DebyeWaveNumber();
            break;
        case (+1): // fermions
            k_D = _fermions[idx].get_DebyeWaveNumber();
            break;
    }
    
    return k_D; // in energy units
}

double Plasma::get_intrinsicFrequency(const long s) const
{
    long iType = _componentType[s];
    long idx = _componentIndex[s];
    double w_i = 0.0;
    switch ( iType ) {
        case (-1): // boson
            w_i = _bosons[idx].get_DebyeWaveNumber();
            break;
        case ( 0): // maxwellian
            w_i = _maxwells[idx].get_DebyeWaveNumber();
            break;
        case (+1): // fermions
            w_i = _fermions[idx].get_FermiEnergy();
            break;
    }
    
    return w_i;
}

double Plasma::get_intrinsicWaveNumber(const long s) const
{
    long iType = _componentType[s];
    long idx = _componentIndex[s];
    double k_i = 0.0;
    switch ( iType ) {
        case (-1): // boson
            k_i = _bosons[idx].get_DebyeWaveNumber();
            break;
        case ( 0): // maxwellian
            k_i = _maxwells[idx].get_DebyeWaveNumber();
            break;
        case (+1): // fermions
            k_i = _fermions[idx].get_FermiWaveNumber();
            break;
    }
    
    return k_i;
}

double Plasma::get_intrinsicTemperature(const long s) const
{
    long iType = _componentType[s];
    long idx = _componentIndex[s];
    double t_i = 0.0;
    switch ( iType ) {
        case (-1): // boson
            t_i = _bosons[idx].get_tempHartree();
            break;
        case ( 0): // maxwellian
            t_i = _maxwells[idx].get_tempHartree();
            break;
        case (+1): // fermions
            t_i = _fermions[idx].get_reducedTemperature();
            break;
    }
    
    return t_i;
}

string Plasma::get_designationString(const long s) const
{
    long iType = _componentType[s];
    long idx = _componentIndex[s];
    string desig;
    switch (iType) {
        case (-1): // boson
            desig = _bosons[idx].get_desig();
            break;
        case (0): // maxwellian
            desig = _maxwells[idx].get_desig();
            break;
        case (+1):
            desig = _fermions[idx].get_desig();
            break;
    }
    
    return desig;
}

complex<double> Plasma::get_chi0(const long s, const double q, const double w, const double t) const
{
    long iType = _componentType[s];
    long idx = _componentIndex[s];
    complex<double> chi0;
    const complex<double> zero(0.0, 0.0);
    switch (iType) {
        case (-1): // boson
            chi0 = zero;
            break;
        case (0): // maxwellian
            chi0 = zero;
            break;
        case (+1):
            chi0 = _fermions[idx].TemperatureFermionComponent::complexLindhard(q, w, t);
            break;
    }
    
    return chi0;
}

complex<double> Plasma::get_chi0(const long s, const double q, const double w) const
{
    long iType = _componentType[s];
    long idx = _componentIndex[s];
    complex<double> chi0;
    const complex<double> zero(0.0, 0.0);
    switch (iType) {
        case (-1): // boson
            chi0 = zero;
            break;
        case (0): // maxwellian
            chi0 = zero;
            break;
        case (+1):
            chi0 = _fermions[idx].FermionComponent::complexLindhard(q, w);
            break;
    }
    
    return chi0;
}
