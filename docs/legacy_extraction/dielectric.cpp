//
//  dielectric.cpp
//  RPA-Plasma
//
//  Created by In-Gee Kim on 12/9/14.
//  Copyright (c) 2014 In-Gee Kim. All rights reserved.
//

#include "Plasma.h"

complex<double> Plasma::get_dielectric(const long unit, const long s, const long t, const double q, const double w) const
{
    // The required unit conversion should be done in the method
    // The method get_chi0() uses only its intrinsic units
    // The reference component is the s-th one.
    
    double qRef, wRef;
    double q_t = 0.0, w_t = 0.0;
    if ( unit == 0 ) { // plasma units of the s-th component
        q_t = get_DebyeWaveNumber(s);
        qRef = q * q_t;
        w_t = get_plasmaFrequency(s);
        wRef = w * w_t;
    } else { // intrinsic units of the s-th component
        q_t = get_intrinsicWaveNumber(s);
        qRef = q * q_t;
        w_t = get_plasmaFrequency(s);
        wRef = w * w_t;
    }

    double eRef, eCouple; // the charges of the s-th and the t-th components, respectively
    eRef = get_charge(s);
    eCouple = get_charge(t);
    
    //
    // set up the potential functions
    //
    Potential potential;
    
    // update the potential parameter list
    vector<double> parameterList = _eePotentialParameter;
    double eMultiple = eRef * eRef;
    parameterList.push_back( eMultiple );
    double eePotential = ( potential.*tblPotential[ _eePotentialMode ] ) ( q, parameterList );
    
    parameterList = _iiPotentialParameter;
    eMultiple = eCouple * eCouple;
    parameterList.push_back( eMultiple );
    double iiPotential = ( potential.*tblPotential[ _iiPotentialMode ] ) ( q, parameterList );

    parameterList = _eiPotentialParameter;
    eMultiple = eRef * eCouple;
    parameterList.push_back( eMultiple );
    double eiPotential = ( potential.*tblPotential[ _eiPotentialMode ] ) ( q, parameterList );
    double eiPotSquare = eiPotential * eiPotential;
    
    //
    // set up the mass weighted temperature
    //
    long idx = -1;
    for ( long is = 0; is < s; is++ )
        for ( long it = 0; it < t; it++ )
            idx++;
    double tau = _weightedTemperature[idx]; // in K
    double k_B = _constantBoltzmann;
    tau *= k_B;  // in Hartree
    tau *= w_t;
    
    //
    // set up chi_0 with the arguments of intrinsic units
    //
    double qE = ( qRef / get_intrinsicWaveNumber(s) );
    double wE = ( wRef / get_intrinsicFrequency(s) );
    double tE = ( tau / get_intrinsicFrequency(s) );
    complex<double> e_chi0 = get_chi0(s, qE, wE, tE);

    double qI = ( qRef / get_intrinsicWaveNumber(t) );
    double wI = ( wRef / get_intrinsicFrequency(t) );
    double tI = ( tau / get_intrinsicFrequency(t) );
    complex<double> i_chi0 = get_chi0(t, qI, wI, tI);

    //
    // calculating the dielectric function epsilon_st
    //
    complex<double> dielectric(0.0, 0.0);
    const complex<double> zone(1.0, 0.0);
    dielectric = ( zone - e_chi0 * eePotential );
    dielectric *= ( zone - i_chi0 * iiPotential );
    dielectric -= ( e_chi0 * i_chi0 * eiPotSquare );
    
    return dielectric;
}

complex<double> Plasma::get_T0dielectric(const long unit, const long s, const long t, const double q, const double w) const
{
    // This calculates the zero-temperature dielectric function
    // The required unit conversion should be done in the method
    // The method get_chi0() uses only its intrinsic units
    // The reference component is the s-th one.
    
    double qRef, wRef;
    double q_t = 0.0, w_t = 0.0;
    if ( unit == 0 ) { // plasma units of the s-th component
        q_t = get_DebyeWaveNumber(s);
        qRef = q * q_t;
        w_t = get_plasmaFrequency(s);
        wRef = w * w_t;
    } else { // intrinsic units of the s-th component
        q_t = get_intrinsicWaveNumber(s);
        qRef = q * q_t;
        w_t = get_plasmaFrequency(s);
        wRef = w * w_t;
    }
    
    double eRef, eCouple; // the charges of the s-th and the t-th components, respectively
    eRef = get_charge(s);
    eCouple = get_charge(t);
    
    //
    // set up the potential functions
    //
    Potential potential;
    
    // update the potential parameter list
    vector<double> parameterList = _eePotentialParameter;
    double eMultiple = eRef * eRef;
    parameterList.push_back( eMultiple );
    double eePotential = ( potential.*tblPotential[ _eePotentialMode ] ) ( q, parameterList );
    
    parameterList = _iiPotentialParameter;
    eMultiple = eCouple * eCouple;
    parameterList.push_back( eMultiple );
    double iiPotential = ( potential.*tblPotential[ _iiPotentialMode ] ) ( q, parameterList );
    
    parameterList = _eiPotentialParameter;
    eMultiple = eRef * eCouple;
    parameterList.push_back( eMultiple );
    double eiPotential = ( potential.*tblPotential[ _eiPotentialMode ] ) ( q, parameterList );
    double eiPotSquare = eiPotential * eiPotential;
    
    //
    // set up chi_0 with the arguments of intrinsic units
    //
    double qE = ( qRef / get_intrinsicWaveNumber(s) );
    double wE = ( wRef / get_intrinsicFrequency(s) );
    complex<double> e_chi0 = get_chi0(s, qE, wE);
    
    double qI = ( qRef / get_intrinsicWaveNumber(t) );
    double wI = ( wRef / get_intrinsicFrequency(t) );
    complex<double> i_chi0 = get_chi0(t, qI, wI);
    
    //
    // calculating the dielectric function epsilon_st
    //
    complex<double> dielectric(0.0, 0.0);
    const complex<double> zone(1.0, 0.0);
    dielectric = ( zone - e_chi0 * eePotential );
    dielectric *= ( zone - i_chi0 * iiPotential );
    dielectric -= ( e_chi0 * i_chi0 * eiPotSquare );
    
    return dielectric;
}

complex<double> Plasma::get_chiRPA(const long unit, const long s, const long t, const double q, const double w) const
{
    // The required unit conversion should be done in the method
    // The method get_chi0() uses only its intrinsic units
    // The reference component is the s-th one.
    
    double qRef, wRef;
    double q_t = 0.0, w_t = 0.0;
    if ( unit == 0 ) { // plasma units of the s-th component
        q_t = get_DebyeWaveNumber(s);
        qRef = q * q_t;
        w_t = get_plasmaFrequency(s);
        wRef = w * w_t;
    } else { // intrinsic units of the s-th component
        q_t = get_intrinsicWaveNumber(s);
        qRef = q * q_t;
        w_t = get_plasmaFrequency(s);
        wRef = w * w_t;
    }
    
    double eRef, eCouple; // the charges of the s-th and the t-th components, respectively
    eRef = get_charge(s);
    eCouple = get_charge(t);
    
    //
    // set up the potential functions
    //
    Potential potential;
    
    // update the potential parameter list
    vector<double> parameterList = _eePotentialParameter;
    double eMultiple = eRef * eRef;
    parameterList.push_back( eMultiple );
    double eePotential = ( potential.*tblPotential[ _eePotentialMode ] ) ( qRef, parameterList );
    
    parameterList = _iiPotentialParameter;
    eMultiple = eCouple * eCouple;
    parameterList.push_back( eMultiple );
    double iiPotential = ( potential.*tblPotential[ _iiPotentialMode ] ) ( qRef, parameterList );
    
    parameterList = _eiPotentialParameter;
    eMultiple = eRef * eCouple;
    parameterList.push_back( eMultiple );
    double eiPotential = ( potential.*tblPotential[ _eiPotentialMode ] ) ( qRef, parameterList );
    
    //
    // set up the mass weighted temperature
    //
    long idx = 0;
    for ( long is = 0; is < s; is++ )
        for ( long it = 0; it < t; it++ )
            idx++;
    double tau = _weightedTemperature[idx]; // in K
    double k_B = _constantBoltzmann;
    tau *= k_B;  // in Hartree
    tau *= w_t;
    
    //
    // set up chi_0 with the arguments of intrinsic units
    //
    double qE = ( qRef / get_intrinsicWaveNumber(s) );
    double wE = ( wRef / get_intrinsicFrequency(s) );
    double tE = ( tau / get_intrinsicFrequency(s) );
    complex<double> e_chi0 = get_chi0(s, qE, wE, tE);
    
    double qI = ( qRef / get_intrinsicWaveNumber(t) );
    double wI = ( wRef / get_intrinsicFrequency(t) );
    double tI = ( tau / get_intrinsicFrequency(t) );
    complex<double> i_chi0 = get_chi0(t, qI, wI, tI);
    
    //
    // calculating the dielectric function epsilon_st
    //
    complex<double> chiRPA(0.0, 0.0);
    
    long id = s + t;
    const complex<double> zone(1.0, 0.0);
    switch ( id ) {
        case (0): // ee
            chiRPA = e_chi0 * ( zone - i_chi0 * iiPotential );
            break;
        case (1): // ei or ie
            chiRPA = e_chi0 * i_chi0 * eiPotential;
            break;
        case (2): // ii
            chiRPA = i_chi0 * ( zone - e_chi0 * eePotential );
            break;
    }
    
    const complex<double> epsilon = get_dielectric(unit, s, t, q, w);
    chiRPA *= conj(epsilon);
    chiRPA /= ( epsilon * conj(epsilon) );
    
    return chiRPA;
}

complex<double> Plasma::get_T0chiRPA(const long unit, const long s, const long t, const double q, const double w) const
{
    // This calculates the zero-temperature chiRPA
    // The required unit conversion should be done in the method
    // The method get_chi0() uses only its intrinsic units
    // The reference component is the s-th one.
    
    double qRef, wRef;
    double q_t = 0.0, w_t = 0.0;
    if ( unit == 0 ) { // plasma units of the s-th component
        q_t = get_DebyeWaveNumber(s);
        qRef = q * q_t;
        w_t = get_plasmaFrequency(s);
        wRef = w * w_t;
    } else { // intrinsic units of the s-th component
        q_t = get_intrinsicWaveNumber(s);
        qRef = q * q_t;
        w_t = get_plasmaFrequency(s);
        wRef = w * w_t;
    }
    
    double eRef, eCouple; // the charges of the s-th and the t-th components, respectively
    eRef = get_charge(s);
    eCouple = get_charge(t);
    
    //
    // set up the potential functions
    //
    Potential potential;
    
    // update the potential parameter list
    vector<double> parameterList = _eePotentialParameter;
    double eMultiple = eRef * eRef;
    parameterList.push_back( eMultiple );
    double eePotential = ( potential.*tblPotential[ _eePotentialMode ] ) ( qRef, parameterList );
    
    parameterList = _iiPotentialParameter;
    eMultiple = eCouple * eCouple;
    parameterList.push_back( eMultiple );
    double iiPotential = ( potential.*tblPotential[ _iiPotentialMode ] ) ( qRef, parameterList );
    
    parameterList = _eiPotentialParameter;
    eMultiple = eRef * eCouple;
    parameterList.push_back( eMultiple );
    double eiPotential = ( potential.*tblPotential[ _eiPotentialMode ] ) ( qRef, parameterList );
    
    //
    // set up chi_0 with the arguments of intrinsic units
    //
    double qE = ( qRef / get_intrinsicWaveNumber(s) );
    double wE = ( wRef / get_intrinsicFrequency(s) );
    complex<double> e_chi0 = get_chi0(s, qE, wE);
    
    double qI = ( qRef / get_intrinsicWaveNumber(t) );
    double wI = ( wRef / get_intrinsicFrequency(t) );
    complex<double> i_chi0 = get_chi0(t, qI, wI);
    
    //
    // calculating the dielectric function epsilon_st
    //
    complex<double> chiRPA(0.0, 0.0);
    
    long id = s + t;
    const complex<double> zone(1.0, 0.0);
    switch ( id ) {
        case (0): // ee
            chiRPA = e_chi0 * ( zone - i_chi0 * iiPotential );
            break;
        case (1): // ei or ie
            chiRPA = e_chi0 * i_chi0 * eiPotential;
            break;
        case (2): // ii
            chiRPA = i_chi0 * ( zone - e_chi0 * eePotential );
            break;
    }
    
    const complex<double> epsilon = get_T0dielectric(unit, s, t, q, w);
    chiRPA *= conj(epsilon);
    chiRPA /= ( epsilon * conj(epsilon) );
    
    return chiRPA;
}
