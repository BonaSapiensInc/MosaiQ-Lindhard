//
//  imaginaryLindhard.cpp
//  rpa-plasma
//
//  Created by In-Gee Kim on 6/19/14.
//  Copyright (c) 2014 In-Gee Kim. All rights reserved.
//

#include "TemperatureFermionComponent.h"

const double TemperatureFermionComponent::squareLindhard(const double q, const double omega, const double tau) const
{
    // q, omega, and tau are in rational units
    double DeF = _DOSatFermiEnergy;
    double real = realLindhard(q, omega, tau) * DeF;
    double imaginary = imaginaryLindhard(q, omega, tau) * DeF;
    
    // result in natural units
    return ( real * real + imaginary * imaginary );
}

const complex<double> TemperatureFermionComponent::complexLindhard(const double q, const double omega, const double tau) const
{
    // q, omega, and tau are in rational units
    complex<double> theLindhard;
    double DeF = _DOSatFermiEnergy;
    
    // result in natural units
    theLindhard.real( realLindhard(q, omega, tau) * DeF );
    theLindhard.imag( imaginaryLindhard(q, omega, tau) * DeF );
    return theLindhard;
}

const double TemperatureFermionComponent::imaginaryLindhard(const double q, const double omega, const double tau) const
{
    double imag = 0.0;
    if(_isMaxwell) {
        imag = get_FCimagLindhard(q, omega, tau); // Fried and Conte (1961) classical
        
    } else {
        // canonical enemble
        // imag = get_KGimaginaryLindhard(q, omega, tau); // Khanna and Glyde (1976)
        
        // grand canonical ensemble
        imag = get_GVimagLindhard(q, omega, tau);  // Giuliani and Vignale (2005)
    }
    
    return imag;
}

const double TemperatureFermionComponent::get_GVimagLindhard(const double q, const double omega, const double tau) const
{
    // Giuliani and Vignale, Quantum Theory of the Electron Liquid (Cambridge University Press, Cambridge, 2005)
    //   Eq. (4.45) p. 173.
    // in rational units with respect to k_F and E_F
    // This is implicitly a grand canonical ensemble representation.
    
    // use the quad-precision functions
    long double nu_p = ( ( omega / q + q ) / 2.0 );
    nu_p *= nu_p;
    long double nu_n = ( ( omega / q - q ) / 2.0 );
    nu_n *= nu_n;
    
    const long double gamma = _reducedChemicalPotential;
    
    long double arg4exp_n = ( ( nu_n - gamma ) / tau );
    long double arg4exp_p = ( ( nu_p - gamma ) / tau );
    
    long double numerator = 1.0 + exp( -arg4exp_n );     // this might be wrong in the GV book
    long double denominator = 1.0 + exp( -arg4exp_p );   // this might be wrong in the GV book
    long double arg4log = numerator / denominator;
    long double theLog = log( arg4log );
    
    double imag = ( tau * ((double) theLog) );
    // imag += omega;     // this might be wrong in the GV book
    imag *= -( _pi / q );

    return imag;
}

const double TemperatureFermionComponent::get_KGimaginaryLindhard(double const q, double const omega, double const tau) const
{
    // Khanna and Glyde, Can. J. Phys. 54 (1976) 648.
    // in rational units with respect to k_F and E_F of each component
    // This is implicitly a canonical ensemble representation.
    
    // use the quad-precision functions
    long double xi_p = ( ( omega / q + q ) / 2.0 );
    long double xi_n = ( ( omega / q - q ) / 2.0 );
    
    const long double gamma = _reducedChemicalPotential;
    
    long double arg4exp_n = ( ( gamma - xi_n * xi_n ) / tau );
    long double arg4exp_p = ( ( gamma - xi_p * xi_p ) / tau );
    
    long double numerator = 1.0 + exp( arg4exp_n );
    long double denominator = 1.0 + exp( arg4exp_p );
    long double arg4log = numerator / denominator;
    long double theLog = log( arg4log );
    
    // result in rational units ; in double precision
    double prefactor = ( ( tau * _pi ) / ( 4.0 * q ) ) ;
    
    return prefactor * ( (double) theLog );
}

const double TemperatureFermionComponent::get_FCimagLindhard(const double q, const double omega, const double tau) const
{
    // The classical distribution approximation
    // by B. D. Fried and S. D. Conte, "The Plasma Dispersion Function" (Academic Press, New York, 1961).
    // implementation based on Eq. (33.8) of
    // Fetter and Walecka, "Quantum theory of many-particle systems" (McGraw-Hill, New York, 1971)
    double x = ( omega / ( 2.0 * tau ) );

    double y = ( (( omega * omega ) / ( q * q )) + ( q * q ) );
    y /= ( 4.0 * tau );
    
    double imag = 0.0;
    imag = sinh(x) / x;
    imag *= exp( -y );
    
    const double pi = _pi;
    imag *= sqrt(pi / (tau * tau * tau));
    imag *= (omega / q);
    imag *= -(2.0 / 3.0);
    
    return imag;
}

const double TemperatureFermionComponent::realLindhard(const double q, const double omega, const double tau) const
{
    double real = 0.0;
    
    if (_isMaxwell) {
        real = get_FCrealLindhard(q, omega, tau); // Fried and Conte (1961)
    } else {
        // canonical ensemble
        // real = get_KKrealLindhard(q, omega, tau); // Kramers-Kronig via Stenger (1981)
        
        // grand canonical ensemble
        real = get_GVrealLindhard(q, omega, tau); // Giuliani and Vignale (2005)
    }
    
    // belows are obsolute
    // real = get_KKrealLindhard(q, omega, tau);  // Kramers-Kronig via Stenger (1981)
    // real = get_KLrealLindhard(q, omega, tau); // Kim and Lee (1998)
    // real = get_DWrealLindhard(q, omega, tau); // Dharma-wardana (1981)
    // real = get_GDrealLindhard(q, omega, tau); // Gouedard and Deutsch (1978)
    // real = get_KGrealLindhard(q, omega, tau); // Khanna and Glyde (1976)
    
    return real;
}

const double TemperatureFermionComponent::get_GVrealLindhard(const double q, const double omega, const double tau) const
{
    // The finite temperature Lindhard function in grand canonical ensemble.
    // Giuliani and Vignale, Quantum Theory of the Electron Liquid (Cambridge University Press, Cambridge, 2005)
    //   Eq. (4.44) p. 172.
    
    // The method for quadrature is based on the sinc quadrature rule by
    // Davis and Rabinowitz, Methods of Numerical Integration Second Edition (Dover Publications, Mineola, 1984)
    //   Eq. (3.4.6.21) p. 218.

    const long N = 256; // convergence parameter
    const long double h = sqrt( (_pi * 7.0) / ((double) N)); // convergence paramter
    
    const double nu_p = ( ( (omega / q) + q ) / 2.0 );
    const double nu_n = ( ( (omega / q) - q ) / 2.0 );
    const double gamma = _reducedChemicalPotential;

    double real = 0.0;
    double kh = 0.0;
    double x = 0.0;
    double argp = 0.0;
    double argn = 0.0;
    double flog = 0.0;
    double sqExp = 0.0;
    for (long k = -N; k <= N; k++) {
        kh = (h * ( (double) k ) );
        x = log( exp(kh) + sqrt( 1.0 + exp( 2.0 * kh) ) );
        
        argn = ( x - nu_n ) / ( x + nu_n );
        argn = sqrt( argn * argn );
        argp = ( x - nu_p ) / ( x + nu_p );
        argp = sqrt( argp * argp );
        
        flog = ( log(argn) - log(argp) );
        flog /= ( exp( ( x * x - gamma ) / tau ) + 1.0 );
        flog *= x;
        sqExp = 1.0 / sqrt(1.0 + exp( - 2.0 * kh ) );
        real += (sqExp * flog);
    }
    real *= h;
    
    real /= -q;
    
    return real;
}

const double TemperatureFermionComponent::get_KKrealLindhard(const double q, const double omega, const double tau) const
{
    // Lindhard function by the Kramers-Kronig transformation
    // in terms of sinc quadrature rule by F. Stanger, J. Approx. Theory 17, 222 (1976),
    // with the formulation appeared in P. J. Davis and P. Rabinowitz,
    //   "Method of Numerical Integration" Second Edition, (Dover, Mineola, 1984) Eq. (3.4.6.16).
    // This is implicitly a canonical ensemble representaion.
    
    const double gamma = _reducedChemicalPotential;
    const double pi = _pi;
    
    const long N = 128; // convergence parameter
    const long double h = sqrt( (pi * 7.0) / ((double) N)); // convergence paramter
    const long double piH = pi / h;

    const long double xi_n = ( ((long double) omega) / ((long double) q) - ((long double) q) ) / ((long double) 2.0);
    const long double xi_p = ( ((long double) omega) / ((long double) q) + ((long double) q) ) / ((long double) 2.0);

    long double kh = 0.0;
    long double arg4exp = 0.0;

    long double f_kh = 0.0;
    
    long double cauchyPrincipalValue_n = 0.0;
    long double cauchyPrincipalValue_p = 0.0;
    for(long k = -N; k <= N; k++) {
        kh = ( ((long double) k) * h );
        arg4exp = ( ((long double) gamma) - kh * kh ) / tau;
        f_kh = log( 1.0 + exp( arg4exp ) ) ;
        
        cauchyPrincipalValue_n += f_kh * (( 1.0 - cos(piH * (xi_n - kh) ) ) / ( piH * (xi_n - kh)) );
        cauchyPrincipalValue_p += f_kh * (( 1.0 - cos(piH * (xi_p - kh) ) ) / ( piH * (xi_p - kh)) );
    }
    
    double real = ((double) ( cauchyPrincipalValue_n - cauchyPrincipalValue_p ) );
    real *= ((pi * tau) / (4.0 * q));
    
    return real;
}

const double TemperatureFermionComponent::get_FCrealLindhard(const double q, const double omega, const double tau) const
{
    double alpha_p = ( (omega / q) + q ) / ( 2.0 * sqrt( tau ) );
    double alpha_n = ( (omega / q) - q ) / ( 2.0 * sqrt( tau ) );
    
    // double real = approxPhi(alpha_p) - approxPhi(alpha_n);
    double real = plasmaDispersion(alpha_p) - plasmaDispersion(alpha_n);
    
    real *= -( 1.0 / ( 3.0 * q * sqrt(tau) ) ) ;
    
    return real;
}

const double TemperatureFermionComponent::approxPhi(double x) const
{
    // approximated plasma dispersion function in Fetter-Wallecka
    double value = 0.0;
    
    if ( x <= 1.0 ) {
        value = ( 1.0 + (1.0 / ( 2.0 * x * x)));
        value /= x;
    } else {
        value = ( 1.0 - ((2.0 * x * x) / 3.0));
        value *= (2.0 * x);
    }
    
    return value;
}

const double TemperatureFermionComponent::plasmaDispersion(double x) const
{
    // plasma dispersion function through the Cauchy principal value integration
    const long N = 128;
    const double pi = _pi;
    const double h = sqrt( (pi * 7.0) / ((double) N) );
    const double piH = pi / h;
    
    double principal = 0.0;
    double kh = 0.0;
    double ySquare = 0.0;
    double f_kh = 0.0;
    
    for (long k = -N; k <= N; k++) {
        kh = ( (double) k ) * h;
        ySquare = kh * kh;
        f_kh = exp( -ySquare );
        
        principal += f_kh * (( 1.0 - cos(piH * (x - kh) ) ) / ( piH * (x - kh)) );
    }
    
    double prefactor = sqrt(pi);
    
    principal *= prefactor;
    
    return principal;
}

const double TemperatureFermionComponent::get_GDrealLindhard(double const q, double const omega, double const tau) const
{
    // C. Gouedard and C. Deutsch, J. Math. Phys. 19 (1978) 32.
    // q, omega, and tau are in rational units
    const double pi = _pi;
    const double gamma = _reducedChemicalPotential;
    const double eta = gamma / tau;
    
    // the zero-momentum and zero-frequency contribution
    double zeroIntegral = theFermiDiracIntegral(-1, 2, eta);
    zeroIntegral *= sqrt( tau );
    zeroIntegral /= 2.0;
    
    // the j-summation
    double r2_j = 0.0;
    double s_j = 0.0;
    double twojone = 0.0;
    double tau_j = 0.0;
    double a_j = 0.0;
    double a2_j = 0.0;
    double b_j = 0.0;
    double b2_j = 0.0;
    
    double twoQ = 2.0 * q;
    double arg_pp = 0.0;
    double arg_pn = 0.0;
    double arg_np = 0.0;
    double arg_nn = 0.0;
    double arctangents = 0.0;
    double ratio = 0.0;
    
    double xi_p = ( omega / q + q ) / 2.0 ;
    double xi_n = ( omega / q - q ) / 2.0 ;
    
    const long jmax = 640; // the number of terms to be summed.
    double jSummation = 0.0;
    for (long j = 0; j < jmax; j++) {
        twojone = ((double) ( 2j + 1 ));
        tau_j = pi * twojone * tau;
        s_j = sqrt(gamma * gamma + tau_j * tau_j);
        a_j = sqrt(gamma + s_j) / sqrt(2.0);
        b_j = sqrt(-gamma + s_j) / sqrt(2.0);
        a2_j = a_j * a_j;
        b2_j = b_j * b_j;
        r2_j = a2_j + b2_j;
        arg_pp = ( xi_p + a_j ) / b_j;
        arg_pn = ( xi_p - a_j ) / b_j;
        arg_np = ( xi_n + a_j ) / b_j;
        arg_nn = ( xi_n - a_j ) / b_j;
        
        arctangents = 0.0;
        arctangents = atan(arg_pp) + atan(arg_pn) - atan(arg_np) - atan(arg_nn);
        arctangents /= twoQ;
        
        ratio = ( b_j / r2_j );
        
        jSummation += (ratio - arctangents) ;
    }
    jSummation *= (pi * tau);
    
    double result = ( zeroIntegral + jSummation );
    result /= -4.0;
    
    return result;
}

const double TemperatureFermionComponent::get_KLrealLindhard(double const q, double const omega, double const tau) const
{
    // Kim and Lee, Phys. Lett. A 245 (1998) 139.
    // in rational units
    const long m = 8;  // the number of residues
    
    double integral = 0.0;
    integral += KLresidueSum(m, q, omega, tau);
    integral += KLresidueSum(m, q, -omega, tau);
    
    integral += KLsumEulerMaclaurin(m, q, omega, tau);
    integral += KLsumEulerMaclaurin(m, q, -omega, tau);
    
    const double prefactor = -3.0 / ( 2.0 * q );
    
    return prefactor * integral;
}

const double TemperatureFermionComponent::KLresidueSum(long const m, double const q, double const omega, double const tau) const
{
    // Kim and Lee, Phys. Lett. A 245 (1998) 139.
    // in rational units
    double pi = _pi;
    double gamma = _reducedChemicalPotential;
    double qHalf = q / 2.0;

    double eta_j = 0.0;
    double phi_j = 0.0;
    double twojone = 0.0;
    double tau_j = 0.0;
    double pitau_j = 0.0;
    complex<double> z_p(0.0, 0.0);
    complex<double> z_n(0.0, 0.0);
    complex<double> kappa(0.0, 0.0);
    
    double sum = 0.0;
    double direct = 0.0;
    for( long j = 0 ; j < ( m - 1 ); j++ ) {
        twojone = (double) ( 2 * j + 1 );
        tau_j = tau * twojone;
        pitau_j = pi * tau_j;
        
        // argument
        if ( gamma > 0 ) {
            phi_j = atan( pitau_j / gamma ) / 2.0 ;
        } else {
            phi_j = ( pi - atan( abs( pitau_j / gamma ) ) ) / 2.0 ;
        }
        // modulus
        eta_j = sqrt( sqrt( gamma * gamma + pitau_j * pitau_j ) );
        
        kappa = polar( eta_j, phi_j );
        
        z_p = qHalf + kappa;
        z_n = qHalf - kappa;
        
        direct = 0.0;
        direct = ( atan2( z_p.imag(), z_p.real() ) - atan2( z_n.imag(), z_n.real() ));
        
        sum += direct;
    }
    
    double twopi = 2.0 * pi;
    sum *= (twopi * tau);
    
    return sum;
}

const double TemperatureFermionComponent::KLsumEulerMaclaurin(long const m, double const q, double const omega, double tau) const
{
    // Kim and Lee, Phys. Lett. A 245 (1998) 139.
    // in rational units
    double pi = _pi;
    double gamma = _reducedChemicalPotential;
    double qHalf = q / 2.0;
    double xi_p = (q + omega/q)/2.0 ;
    double xi2_p = xi_p * xi_p;

    double sum = 0.0;
    
    // J0
    double tau_m = tau * ( (double) ( 2 * m + 1 ) );
    double arg4sqrt = gamma * gamma + ( pi * tau_m ) * ( pi * tau_m );
    double eta_m = sqrt( sqrt( arg4sqrt ) );
    double phi_m;
    if ( gamma < 0 ) {
        phi_m = atan( pi * tau_m / gamma ) / 2.0;
    } else {
        phi_m = pi - atan( pi * tau_m / gamma ) ;
        phi_m /= 2.0;
    }
    double j01 = 2.0 * xi_p * eta_m * cos( phi_m );
    
    double num4log, den4log, sin4log;
    sin4log = ( eta_m * sin( phi_m ) ) * ( eta_m * sin( phi_m ) );
    num4log = (( xi_p + eta_m * cos( phi_m ) ) * ( xi_p + eta_m * cos( phi_m ) )) + sin4log;
    den4log = (( xi_p - eta_m * cos( phi_m ) ) * ( xi_p - eta_m * cos( phi_m ) )) + sin4log;
    double arg4log = num4log / den4log;
    double j02log = log( arg4log );
    double j02 = (( gamma - xi2_p ) / 2.0) * j02log;
    double j0 = j01 + j02;
    sum += j0;
    
    // J1
    complex<double> kappa( polar(eta_m, phi_m) );
    complex<double> z_p = ( qHalf + kappa );
    complex<double> z_n = ( qHalf - kappa );
    double pitau_m = ((double) (2 * m + 1)) * pi * tau;
    double j1 = -2.0 * pitau_m * ( atan2( z_p.imag(), z_p.real() ) + atan2( z_n.imag(), z_n.real() ) );
    sum += j1;
    
    // J2
    double j2Numerator = ( xi2_p - gamma ) * cos( phi_m ) + pi * tau * sin( phi_m );
    double j2Denominator = ( xi2_p - gamma ) * ( xi2_p  - gamma ) + ( pi * tau ) * ( pi * tau ) ;
    double j2r = ( pi * tau ) * ( pi * tau ) * ( j2Numerator / j2Denominator );
    double j2 = -xi_p / ( 3.0 * eta_m ) * j2r;
    sum += j2;
    
    return sum;
}

const double TemperatureFermionComponent::get_KGrealLindhard(double const q, double const omega, double const tau) const
{
    // Eq. (20) of F. C. Khanna and H. R. Glyde, J. Can. Phys. 54 (1976) 648.
    // q, omega, and tau are in rational units
    const double pi = _pi;
    const double gamma = _reducedChemicalPotential;
    const double eta = gamma / tau;
    double gamma2 = gamma * gamma;
    double qHalf = q / 2.0;
    double qSquare = q * q;
    double qInv = 1.0 / q;
    double x = gamma - ( qHalf * qHalf );
    double y = gamma + ( qHalf * qHalf );
    
    // the zero-momentum and zero-frequency contribution
    double zeroIntegral = theFermiDiracIntegral(-1, 2, eta);
    zeroIntegral *= sqrt( tau );
    zeroIntegral /= 2.0;
    
    // the residue contributions
    long jmax = 640;
    
    double xi_p = ( omega / q + q ) / 2.0;
    double xi_n = ( omega / q - q ) / 2.0;
    double xi2_p = xi_p * xi_p;
    double xi2_n = xi_n * xi_n;
    
    // the j-summation
    double staticSum = 0.0;
    double staticCorrect = 0.0;
    
    double staticIntg1 = 0.0;
    double staticIntg2 = 0.0;
    double staticIntg3 = 0.0;
    double staticIntg4 = 0.0;
    double staticIntg5 = 0.0;
    double staticIntg6 = 0.0;
    
    // double dynamicSum = 0.0;
    double dynamicSum1 = 0.0;
    double dynamicSum2 = 0.0;
    // double dynamicIntg = 0.0;
    double dynamicIntg1 = 0.0;
    double dynamicIntg2 = 0.0;
    
    double twojone = 0.0;
    double pitau_j = 0.0;
    double pitau2_j = 0.0;
    double eta_j = 0.0;
    double kr_j = 0.0;
    double ki_j = 0.0;
    double k2_j = 0.0;
    double ks_j = 0.0;
    
    double qnlog = 0.0;
    double qplog = 0.0;
    double nlog = 0.0;
    double plog = 0.0;
    for (long j = 0; j < jmax; j++) {
        twojone = ((double) (2j + 1));
        pitau_j = pi * twojone * tau;
        pitau2_j = pitau_j * pitau_j;
        eta_j = sqrt(gamma2 + pitau2_j);
        kr_j = sqrt(gamma + eta_j) / sqrt(2.0);
        ki_j = sqrt(-gamma + eta_j) / sqrt(2.0);
        k2_j = kr_j * kr_j + ki_j * ki_j;
        ks_j = sqrt(k2_j);
        
        // the static contributions
        // the first integral
        staticIntg1 = ( kr_j / k2_j ) * (qSquare / 3.0);
        // the second integral
        staticIntg2 = 2.0 * kr_j * (log( (x * x) + pitau2_j) );
        // the first part of the third integral
        staticIntg3 = ( (2.0 * kr_j) / k2_j ) *
        ( ( (3.0/2.0) * qHalf * qHalf - 5.0 * gamma) * log( (x * x) + pitau2_j)
         - 2.0 * x - qHalf * ( 3.0 * gamma - qHalf * qHalf) - 4.0 * atan( x / pitau_j ) );
        // the second part of the third integral
        if ( (qHalf * qHalf - (gamma - pitau_j) ) > 0.0 ) {
            qnlog = log(qHalf * qHalf - (gamma - pitau_j));
        } else {
            qnlog = 0.0;
        }
        if ( (qHalf * qHalf - (gamma + pitau_j) ) > 0.0 ) {
            qnlog = log(qHalf * qHalf - (gamma + pitau_j));
        } else {
            qplog = 0.0;
        }
        if ( (gamma - pitau_j) > 0.0 ) {
            nlog = log(gamma - pitau_j);
        } else {
            nlog = 0.0;
        }
        if ( (gamma - pitau_j) > 0.0 ) {
            plog = log(gamma + pitau_j);
        } else {
            plog = 0.0;
        }
        staticIntg4 = (qInv * (64.0 * pitau2_j) / k2_j) *
        ( sqrt(gamma - pitau_j) * (qnlog - nlog)
         - sqrt(gamma + pitau_j) * (qplog - plog) );
        // the third part of the third integral
        staticIntg5 = (4.0 * kr_j * ki_j / k2_j ) * log((x * x) + pitau2_j);
        
        // summing up
        staticSum += ((staticIntg1 - staticIntg2 - staticIntg3 - staticIntg4 - staticIntg5) / pitau_j);
        
        // the dynamic contribution
        // static correction from the frequency integral
        staticIntg6 = (kr_j / k2_j) * (log((x * x) + pitau2_j) - log((y * y) + pitau2_j));
        staticCorrect += staticIntg6;
        
        // the first dynamic integral
        if ( (xi2_n - (gamma - pitau_j)) > 0.0 ) {
            qnlog = log(xi2_n - (gamma - pitau_j));
        } else {
            qnlog = 0.0;
        }
        if ( (xi2_p - (gamma - pitau_j)) > 0.0 ) {
            qplog = log(xi2_p - (gamma - pitau_j));
        } else {
            qplog = 0.0;
        }
        if ( (xi2_n - (gamma + pitau_j)) > 0.0 ) {
            nlog = log(xi2_n - (gamma + pitau_j));
        } else {
            nlog = 0.0;
        }
        if ( (xi2_p- (gamma + pitau_j)) > 0.0 ) {
            plog = log(xi2_p - (gamma + pitau_j));
        } else {
            plog = 0.0;
        }
        dynamicIntg1 = ((pitau2_j * kr_j)/k2_j) *
        (sqrt(gamma - pitau_j) * (qnlog - qplog)
         - sqrt(gamma + pitau_j) * (nlog - plog));
        dynamicSum1 += dynamicIntg1;
        
        // the second dynamic integral
        dynamicIntg2 = (kr_j / k2_j) *
        (xi_n * log((gamma - xi2_n) * (gamma - xi2_n) + pitau2_j)
         - xi_p * log((gamma - xi2_p) * (gamma - xi2_p) + pitau2_j));
        dynamicSum2 += dynamicIntg2;
        
        // The original dynamic sum
        // dynamicIntg = atan((xi_p * ki_j)/(xi2_p - ks_j)) - atan((xi_n * ki_j)/(xi2_n - ks_j));
        // dynamicSum += dynamicIntg;
    }
    // multiplication factors
    staticSum *= (pi * tau);
    staticCorrect *= ((pi * pi * tau * tau) / 4.0);
    dynamicSum1 *= (4.0 * pi * tau * qInv);
    dynamicSum2 *= ((pi * pi * tau * tau * qInv)/ 2.0);
    // dynamicSum /= (2.0 * q);
    
    double result = zeroIntegral - staticSum - staticCorrect + dynamicSum1 - dynamicSum2;
    return result;
}

