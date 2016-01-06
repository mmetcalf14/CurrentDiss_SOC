//
//  EQ_Current.hpp
//  CurrentDissipation_1DRing
//
//  Created by mekena McGrew on 1/6/16.
//  Copyright © 2016 Mekena Metcalf. All rights reserved.
//

#ifndef EQ_Current_hpp
#define EQ_Current_hpp
#include "Distribution.hpp"
#include <stdio.h>
//no real need for a header file here

class EQ_Current
{

private:

    double jn_up;
    double jn_dn;
    double jn;

    double Tot_Current;

    std::vector<double> Phi_Table;
    std::vector<double> PC_EQ_Table;
    std::vector<double> PC_NEQ_Table;

public:
    EQ_Current(const Distribution&) {};//no idea if this is correct


    void Full_EQ_PC(const Distribution &FD,  double jo, double ** FDu_Mat, double ** FDd_Mat, double ** FDu_Neg, double ** FDd_Neg);
    void Full_NEQ_PC(const Distribution &FD,  double jo, double ** FDu_Mat, double ** FDd_Mat, double ** FDu_Neg, double ** FDd_Neg);
    void CurrentOut(const Distribution &FD, std::ofstream &fout);


};


#endif /* EQ_Current_hpp */
