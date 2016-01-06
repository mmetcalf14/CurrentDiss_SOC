//
//  Distribution.hpp
//  CurrentDissipation_1DRing
//
//  Created by mekena McGrew on 1/6/16.
//  Copyright © 2016 Mekena Metcalf. All rights reserved.
//

#ifndef Distribution_hpp
#define Distribution_hpp

#include <stdio.h>
//
//  RootFinder.h
//  CurrentDissipation_1DRing
//
//  Created by mekena McGrew on 12/11/15.
//  Copyright © 2015 Mekena Metcalf. All rights reserved.
//

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>


struct TotalNumber_params;

class Distribution
{
    friend class EQ_Current;
    friend class DissipatedIntegral;
    
private:
    double Temp;
    double Eo;
    double w_soc;
    double To;
    double dt;
    double N;
    
    int Nl;
    int Nt;
    int Tit;
    
    TotalNumber_params * p;
    
    std::vector<double> mu;
public:
    //Distribution(TotalNumber_params input_paras) {my_paras = input_paras;};
    Distribution(int _Nl, double _N, double _T, double _E, double _w, double _dt, double _To);
    
    double TotalNumber(double x, void *par);
    void Set_Params(TotalNumber_params * xp) {p = xp;};
    void ChemicalPotential();
    void FD_Distribution( int time_rev, int species, std::ofstream &fout, double ** FD_Mat);
};


struct TotalNumber_params
{
    int  Nl;//I want all of these in the class
    double N, T, Phi, Eo, w;//Initial conditions to solve function
    
    Distribution * pt_Dist;
};

double gslClassWrapper(double x, void * pp)
{
    TotalNumber_params *p = (TotalNumber_params *)pp;
    return p->pt_Dist->TotalNumber(x,p);
}






#endif /* Distribution_hpp */
