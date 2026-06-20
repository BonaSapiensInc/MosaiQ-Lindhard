//
//  TemperatureFermionComponent.cpp
//  rpa-plasma
//
//  Created by In-Gee Kim on 6/12/14.
//  Copyright (c) 2014 In-Gee Kim. All rights reserved.
//

#include "TemperatureFermionComponent.h"

TemperatureFermionComponent::TemperatureFermionComponent()
: _constantBoltzmann(3.1668114e-6)
{
    // default constructor is left as an empty-body
}

TemperatureFermionComponent::~TemperatureFermionComponent()
{
    // default destructor is left as an empty-body
}

// preferred constructor
TemperatureFermionComponent::TemperatureFermionComponent(const char name, const double mass, const double charge,
                                                         const double radiusWS, const double temperature, const bool isMaxwell)
: FermionComponent(name, mass, charge, radiusWS), _constantBoltzmann(3.1668114e-6)
{
    _temperature = temperature; // in units of K
    _tempHartree = _constantBoltzmann * temperature; // in units of Hartree
    _reducedTemperature = _tempHartree / _FermiEnergy; // reduced temperature
    _isMaxwell = isMaxwell;
    
    double ZSquare = _charge * _charge;
    double rsCube = radiusWS * radiusWS * radiusWS;
    _DebyeWaveNumber = sqrt( ( 3.0 * ZSquare * 311775.0 ) / ( temperature * rsCube ) );
   
    set_ChemicalPotential();
}

void TemperatureFermionComponent::printStatus(ofstream &logFile) const
{
    // prints the current status of the Fermion componont
    
    string indent = "   ";
    
    logFile << indent << "designation = " << _desig << endl;
    
    logFile << indent << "mass = " << _mass << " m_e" << endl;
    logFile << indent << "charge = " << _charge << " |e|" << endl << endl;
    
    logFile << indent << "Wigner-Seitz radius = "  << setprecision(9)
        << _radiusWignerSeitz << " a_0" << endl;
    logFile << indent << "Plasma frequency = " << setprecision(9)
        << _plasmaFrequency << " Hartree" << endl << endl;
    
    logFile << indent << "Fermi wavevector (k_F) = " << setprecision(9)
        << _FermiWavevector << " a_0^-1" << endl;
    logFile << indent << "Fermi energy (E_F) = "  << setprecision(9)
        << _FermiEnergy << " Hartree" << endl;
    logFile << indent << "DOS at E_F = "  << setprecision(9)
        << _DOSatFermiEnergy << " states/Hartree" << endl << endl;
    
    logFile << indent << "temperature (T) = "  << setprecision(9) << _temperature << " K" << endl;
    logFile << indent << "temperature (k_B T) = " << setprecision(9) << _tempHartree << " Hartree" << endl;
    logFile << indent << "reduced temperature (\\tau = k_B T / E_F) = "  << fixed << setprecision(9)
        << _reducedTemperature << endl;
    logFile << indent << "Debye wave number (k_D) = " << setprecision(9)
    << _DebyeWaveNumber << " a_0^-1" << endl;
    logFile << indent << "chemical potential (\\mu) = "  << setprecision(9)
        << _chemicalPotential << " Hartree" << endl;
    logFile << indent << "reduced chemical potential (\\mu / E_F) = "  << setprecision(9)
        << _reducedChemicalPotential << endl;
    string dist = "Fermi-Dirac";
    if (_isMaxwell)
        dist = "Maxwell-Boltzmann";
    logFile << indent << "distribution = " << dist << endl;
}

void TemperatureFermionComponent::set_ChemicalPotential( void )
{
    // obtains the chemical potential of "this" fermion component
    // by employing the Newton-Raphson method with the Fermi-Dirac integrals.
    
    // preparation for the Newton-Raphson inputs
    double chemicalPotential = _FermiEnergy;  // initial guess
    double gamma = chemicalPotential / _FermiEnergy;
    double tau = _reducedTemperature;
    
    double density = _density;
    double densityPrefactor;
    const double twoThird = 3.0 / 2.0;
    const double kFCube = _FermiWavevector * _FermiWavevector * _FermiWavevector;
    const double twoPiSquare = 2.0 * _pi * _pi;
    densityPrefactor = pow(tau, twoThird ) * kFCube / twoPiSquare;
    
    // these two dimensionless variables are used
    double eta = gamma / tau;
    double y = density / densityPrefactor;
    
    // the Newton-Raphson body
    double etaOld;
    double etaNew;
    double delta;
    long iter = 0;

    const long iterMax = 100;
    const double epsilon = 1.0e-14;
    
    etaOld = etaNew = eta; // the initial guess
    do {
        etaOld = etaNew;
        double F = 0.0;
        F = theFermiDiracIntegral(1, 2, etaOld);
        double Fd = 0.0;
        Fd = theFermiDiracIntegral(-1, 2, etaOld);

        delta = (y / Fd) - (F / Fd) ;
        etaNew = etaOld + delta;
    } while ( ( abs(delta) > epsilon ) && ( iter++ < iterMax ) );
    eta = etaNew;
    
    // assign the results
    gamma = eta * tau;
    _reducedChemicalPotential = gamma; // update the reduced chemical potential
    chemicalPotential = gamma * _FermiEnergy;
    _chemicalPotential = chemicalPotential;  // update the chemical potential
    
    // Sommerfeld expansion for comparison
    // cerr << "Chemical potential = " << chemicalPotential << endl;
    // cerr << "Sommerfeld chemical potential = "
    //     << _fermiEnergy * ( 1.0 - (1.0/3.0) * pow( _pi * _constantBoltzmann * _temperature / (2.0 * _fermiEnergy), 2.0) )
    // << endl;
}

const double TemperatureFermionComponent::theFermiDiracIntegral(long orderNumerator, long orderDenominator, double eta) const
{
    // calculate the Fermi-Dirac integral
    // the order of the FD integral is calculated by orderNumerator/orderDenominator
    double order;
    order = ((double) orderNumerator) / ((double) orderDenominator);

    double x0; // parameter for integration
    
    bool isNegativeEta = signbit(eta);  // true if signbit(x) is negative.

    if ( isNegativeEta ) {
        x0 = 8.0;
    } else {
        x0 = sqrt( eta + 64.0 );
    }
    const long N = 128; // define the number of sinc nodes
    const double h = x0 / N; // define the abscissa step
    
    //
    // the trapezoidal part
    //
    double F = 0.0; // the integral
    
    // the first term
    double x = 0.0;
    F = h * ( pow(x, 2.0 * order + 1.0) / (exp( x * x - eta ) + 1.0) ) ;
    
    // the summation
    for ( long n = 1; n < N; n++ ) {
        double x;
        x = ((double) n ) * h;
        F += ( pow(x, 2.0 * order + 1.0) / (exp( x * x  - eta ) + 1.0) ) ;
    }
    F = 2.0 * F;
    
    //
    // the pole correction part
    //
    long kMax = 0;  // the number of poles to be included
    double kMaxSquare;
    double d = h * 7.0;  // the width of the contour strip
    
    double dSquare;
    double etaSquare;
    const double piSquare = _pi * _pi ;
    const double twoPi = 2.0 * _pi;
    double twoDSquarePlusEta;
    double twoDSquarePlusEtaSquare;
    
    dSquare = d * d;
    twoDSquarePlusEta = 2.0 * dSquare + eta;
    twoDSquarePlusEtaSquare = twoDSquarePlusEta * twoDSquarePlusEta ;
    etaSquare = eta * eta;
    
    kMaxSquare = ( ( twoDSquarePlusEtaSquare - etaSquare ) / ( 4.0 * piSquare ) );
    kMax = (long) ( sqrt( kMaxSquare ) );

    // sign(eta) dependent of the prefactor for the pole correction
    double polePrefactor;
    if ( isNegativeEta ) {
        polePrefactor = +4.0 * _pi ;
    } else {
        polePrefactor = -4.0 * _pi ;
    }
    
    // the pole correction summation
    double Pjh = 0.0;
    for ( long k = 0 ; k <= kMax ; k++) {
        double dk = (double) k;
        double oddk = 2.0 * dk + 1.0;
        double rho_k = sqrt( etaSquare + piSquare * ( oddk * oddk ) );
        double theta_k;
        if ( isNegativeEta ) {
            theta_k = _pi + atan( oddk * _pi / eta);
        } else {
            theta_k = atan( oddk * _pi / eta );
        }
        double delta_k;
        delta_k = ( twoPi / h ) * sqrt( ( rho_k + eta ) / 2.0 );
        double lambda_k;
        lambda_k = exp ( ( twoPi / h ) * sqrt( ( rho_k - eta ) / 2.0 ) );
        double numeratorPjh;
        double denominatorPjh;
        numeratorPjh = lambda_k * sin( delta_k + order * theta_k ) - sin( order * theta_k );
        denominatorPjh = 1.0 + lambda_k * lambda_k - 2.0 * lambda_k * cos( delta_k );
        Pjh += pow( rho_k, order ) * ( numeratorPjh / denominatorPjh ) ;
    }
    Pjh *= polePrefactor;
    
    return ( F + Pjh );
}