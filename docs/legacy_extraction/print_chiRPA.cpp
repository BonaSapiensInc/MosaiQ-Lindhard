//
//  print_chiRPA.cpp
//  RPA-Plasma
//
//  Created by In-Gee Kim on 12/10/14.
//  Copyright (c) 2014 In-Gee Kim. All rights reserved.
//

#include "print_chiRPA.h"

void Plasma::print_chiRPA(const Plot &params, ofstream &logFile) const
{
    long nComponent = get_nComponent();
    long unit = params.get_unitSystem();
    
    const string unitName[2] = { "plasma", "intrinsic" };
    complex<double> chiRPA(0.0, 0.0);
    
    long qDiv = params.get_qDivision();
    long wDiv = params.get_omegaDivision();
    
    for ( long s = 0; s < nComponent; s++ ) {
        
        long iType = _componentType[s];
        
        for ( long t = 0; t < nComponent; t++ ) {
            string filename = "chiRPA_";
            filename += get_designationString(s);
            filename += get_designationString(t);
            filename += ".dat";
            
            // open a text file for (s, t)-dielectric function
            ofstream textFile(filename);
            
            //
            // main body for the plotting is below
            //
            
            // write headers
            textFile << "#" << endl;
            textFile << "# RPA susceptibility function of the (" << get_designationString(s)
            << ", " << get_designationString(t) << ") components." << endl;
            textFile << "#" << endl;
            textFile << "# The system of units for plotting is " << unitName[ unit ] << " units." << endl;
            if ( unit == 0 ) {  // plasma units
                textFile << "#     q (k_D)    "  << "    w (w_p)     "   <<  "        real        "  <<  "         imaginary" << endl;
            } else if ( (unit == 1) && (iType == +1) ) {  // quantum fermion units
                textFile << "#     q (k_F)    "  << "    w (E_F)     "   <<  "        real        "  <<  "         imaginary" << endl;
            } else {  // quantum boson units
                textFile << "#     q (k_B)    "  << "    w (E_B)     "   <<  "        real        "  <<  "         imaginary" << endl;
            }
            textFile << "#" << endl;
            
            double q = 0.0;
            double w = 0.0;
            for ( long iq = 1; iq < qDiv; iq++ ) {
                q = params.get_qValue(iq);
                for ( long iw = 0 ; iw < wDiv; iw++ ) {
                    w = params.get_wValue(iw);
                    
                    chiRPA = get_T0chiRPA(unit, s, t, q, w);  // unit conversions will be done inside the function
                    
                    textFile << fixed << setw(15) << setprecision(9) << q << ' ';
                    textFile << fixed << setw(15) << setprecision(9) << w << ' ';
                    textFile << scientific << setw(24) << setprecision(15) << chiRPA.real() << ' ';
                    textFile << scientific << setw(24) << setprecision(15) << chiRPA.imag() << endl;
                }
                
                textFile << endl;
                
            }
            
            // close the text file for (s, t)-dielectric function
            textFile.close();
        } // t-loop
    } // s-loop
}