//
//  Potential.h
//  rpa-plasma
//
//  Created by In-Gee Kim on 7/9/14.
//  Copyright (c) 2014 In-Gee Kim. All rights reserved.
//

#ifndef __rpa_plasma__Potential__
#define __rpa_plasma__Potential__

// system header(s)
#include <cmath>
using std::atan;
using std::pow;
using std::sqrt;

#include <vector>
using std::vector;

// user defined header(s)

class Potential {
public:
    // default constructor
    Potential(void);
    // default destructor
    ~Potential(void);
    
    double constant(double, vector<double> &);
    double yukawa(double, vector<double> &);
    double lennard_jones(double, vector<double> &);
    double modified_coulomb(double, vector<double> &);
    
private:
    const double _Ehe;  // statvolts in E_H / e
};

#endif /* defined(__rpa_plasma__Potential__) */
