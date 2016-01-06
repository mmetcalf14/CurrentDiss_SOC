//
//  Diss_Integral.hpp
//  CurrentDissipation_1DRing
//
//  Created by mekena McGrew on 1/6/16.
//  Copyright © 2016 Mekena Metcalf. All rights reserved.
//

#ifndef Diss_Integral_hpp
#define Diss_Integral_hpp

#include <stdio.h>
#include "EQ_Current.hpp"
class DissipatedIntegral
{

private:

    int Nint;
    int Tint;
    double Tau;
    double dt;
    //how can I put the dynamic 2D array in here?

    std::vector<double> Phi_Diss_Table;
    std::vector<double> PC_Diss_Table;


public:

    DissipatedIntegral(const Distribution&) {};

    inline void GetRelaxTime(double tt) {Tau = tt;};
    double ** DissIntegral(const Distribution &FD, double ** FD_Pos, double ** FD_Neg);//I need to call this twice. Once for spin up and once for spin down
    //out put FD func to main and reinsert in current function
    void DissipatedCurrent(const Distribution &FD, std::ofstream& PC_Diss_out, double jo, double ** FDu_Diss, double ** FDd_Diss);
    inline void NintVal() {Nint = 4;};
    void CurveArea();

};

#endif /* Diss_Integral_hpp */
