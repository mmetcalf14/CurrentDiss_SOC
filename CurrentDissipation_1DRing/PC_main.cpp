//
//  main.cpp
//  CurrentDissipation_1DRing
//
//  Created by mekena McGrew on 12/11/15.
//  Copyright © 2015 Mekena Metcalf. All rights reserved.
//

#include <iostream>
//#include "Distribution.h"
//#include "EQ_Current.h"
#include "Diss_Integral.hpp" //all header files link through this one
//trying to avoid multiple compiling errors
using namespace std;

double ** GetArray(std::ifstream& fin, int Td, int Nt );

int main(int argc, const char * argv[])
{
    
    
    double N = 100.;
    double Temp = 0.;//Tf = 10000 =Ef
    
    double Eo = 1.0;
    int Nl = 400;
    double jo = 1.0;
    double w_so = 0.;
    
    double dt = 0.001;
    double To = 1.;
    double tau = 1.*To;
    int Tm = To/dt;
    int Ntable = (2*Nl)+1;
    
    int Tdim1 = Tm +1;
    
    
    //Arrays for Fermi Function with positive and negative flux
    double ** Fn_up_pos;
    Fn_up_pos = new double *[Tdim1];
    for(int i = 0; i <= Tdim1; i++)
        Fn_up_pos[i] = new double[Ntable];
    
    double ** Fn_dn_pos;
    Fn_dn_pos = new double *[Tdim1];
    for(int i = 0; i <= Tdim1; i++)
        Fn_dn_pos[i] = new double[Ntable];
    
    double ** Fn_up_neg;
    Fn_up_neg = new double *[Tdim1];
    for(int i = 0; i <= Tdim1; i++)
        Fn_up_neg[i] = new double[Ntable];
    
    double ** Fn_dn_neg;
    Fn_dn_neg = new double *[Tdim1];
    for(int i = 0; i <= Tdim1; i++)
        Fn_dn_neg[i] = new double[Ntable];
    
    double ** Fn_up_Relaxed;  //it would be better if they were contained in the class
    Fn_up_Relaxed = new double *[Tdim1];
    for(int i = 0; i <= Tdim1; i++)
        Fn_up_Relaxed[i] = new double[Ntable];
    
    double ** Fn_dn_Relaxed;
    Fn_dn_Relaxed = new double *[Tdim1];
    for(int i = 0; i <= Tdim1; i++)
        Fn_dn_Relaxed[i] = new double[Ntable];
    
    
    //input and output files for distribution functions and current
    ifstream fiup;
    fiup.open("FnArray_Sup_pos_N100_Nl400_T0_To1_w0_010416.dat");
    fiup.precision(8);
    fiup.setf(ios::fixed);
    fiup.setf(ios::showpoint);
    
    ifstream fidn;
    fidn.open("FnArray_Sdn_pos_N100_Nl400_T0_To1_w0_010416.dat");
    fidn.precision(8);
    fidn.setf(ios::fixed);
    fidn.setf(ios::showpoint);
    
    ifstream fiup_neg;
    fiup_neg.open("FnArray_Sup_neg_N100_Nl400_T0_To1_w0_010416.dat");
    fiup_neg.precision(8);
    fiup_neg.setf(ios::fixed);
    fiup_neg.setf(ios::showpoint);
    
    ifstream fidn_neg;
    fidn_neg.open("FnArray_Sdn_neg_N100_Nl400_T0_To1_w0_010416.dat");
    fidn_neg.precision(8);
    fidn_neg.setf(ios::fixed);
    fidn_neg.setf(ios::showpoint);
    
    ofstream fout_up_pos;
    fout_up_pos.precision(8);
    fout_up_pos.setf(ios::fixed);
    fout_up_pos.setf(ios::showpoint);
    
    ofstream fout_dn_pos;
    fout_dn_pos.precision(8);
    fout_dn_pos.setf(ios::fixed);
    fout_dn_pos.setf(ios::showpoint);
    
    ofstream fout_up_neg;
    fout_up_neg.precision(8);
    fout_up_neg.setf(ios::fixed);
    fout_up_neg.setf(ios::showpoint);
    
    ofstream fout_dn_neg;
    fout_dn_neg.precision(8);
    fout_dn_neg.setf(ios::fixed);
    fout_dn_neg.setf(ios::showpoint);
    
    ofstream EQ_PC_Out;
    EQ_PC_Out.open("EquilibriumCurrent_SOC_N100_Nl400_T0_To1_w0_010416.dat");
    EQ_PC_Out.precision(8);
    EQ_PC_Out.setf(ios::fixed);
    EQ_PC_Out.setf(ios::showpoint);
    
    ofstream Diss_PC_Out;
    Diss_PC_Out.open("DissipatedCurrent_SOC_N100_Nl400_T0_tau0_To1_w0_010416.dat");
    Diss_PC_Out.precision(8);
    Diss_PC_Out.setf(ios::fixed);
    Diss_PC_Out.setf(ios::showpoint);
    
    //Distribution Dist(Nl, N, Temp, Eo, w_so, dt, To);
    //EQ_Current EQ_PC(Dist);
    //DissipatedIntegral Diss_PC(Dist);
    
    
    if (fiup.is_open() && fidn.is_open())
    {
        cout << "Files exists \n";
        Fn_up_pos = GetArray(fiup, Tdim1, Ntable);
        Fn_dn_pos = GetArray(fidn, Tdim1, Ntable);
        Fn_up_neg = GetArray(fiup_neg, Tdim1, Ntable);
        Fn_dn_neg = GetArray(fidn_neg, Tdim1, Ntable);
        
    }
    else
    {
        
        fout_up_pos.open ("FnArray_Sup_pos_N100_Nl400_T0_To1_w0_010416.dat");
        fout_dn_pos.open ("FnArray_Sdn_pos_N100_Nl400_T0_To1_w0_010416.dat");
        fout_up_neg.open ("FnArray_Sup_neg_N100_Nl400_T0_To1_w0_010416.dat");
        fout_dn_neg.open ("FnArray_Sdn_neg_N100_Nl400_T0_To1_w0_010416.dat");
        
        //        Dist.ChemicalPotential();//I don't like params as an input
        //        Dist.FD_Distribution(0, 1, fout_up_pos, Fn_up_pos);//postive phi spin up
        //        Dist.FD_Distribution(0, 2, fout_dn_pos, Fn_dn_pos);//postive phi spin down
        //        Dist.FD_Distribution(1, 1, fout_up_neg, Fn_up_neg);//negative phi spin up
        //        Dist.FD_Distribution(1, 2, fout_dn_neg, Fn_dn_neg);//negative phi spin down
        
        //        EQ_PC.Full_EQ_PC(Dist, jo, Fn_up_pos, Fn_dn_pos, Fn_up_neg, Fn_dn_neg);
        //        EQ_PC.Full_NEQ_PC(Dist, jo, Fn_up_pos, Fn_dn_pos, Fn_up_neg, Fn_dn_neg);
        //        EQ_PC.CurrentOut(Dist, EQ_PC_Out);
        
        fout_up_pos.close();
        fout_dn_pos.close();
        fout_up_neg.close();
        fout_dn_neg.close();
        EQ_PC_Out.close();
        
    }
    
    //Diss_PC.GetRelaxTime(tau);
    
    //    Fn_up_Relaxed = Diss_PC.DissIntegral(Dist, Fn_up_pos, Fn_up_pos);
    //    Fn_dn_Relaxed = Diss_PC.DissIntegral(Dist, Fn_dn_pos, Fn_dn_pos);
    //
    //
    //    Diss_PC.DissipatedCurrent(Dist, Diss_PC_Out, jo, Fn_up_Relaxed, Fn_dn_Relaxed);
    
    
    Diss_PC_Out.close();
    
    
    return 0;
}

double ** GetArray(std::ifstream& fin, int Td, int Nt ){
    
    double ** Array;
    Array = new double *[Td];
    for(int i = 0; i <= Td; i++)
        Array[i] = new double[Nt];
    
    for (int row = 0; row < Td; row++)
    {
        for ( int col=0; col < Nt; col++)
        {
            fin >> Array[row][col];
        }
    }
    
    fin.close();
    
    return Array;
}
//
//void Distribution::ChemicalPotential()
//{
//
//    double t = 0.;
//    double Phi = 0;
//
//
//    for(int i = 0; i <= Tit; i++)
//    {
//        Phi = t/(2*To);
//        int status;
//        int iter = 0, max_iter = 1000;
//        const gsl_root_fsolver_type *T;
//        gsl_root_fsolver *s;
//        double x_lo = -2000.0;
//        double x_hi = 3000.0;
//
//        gsl_function F;
//        struct TotalNumber_params params = { Nl, N, Temp, Phi, Eo, w_soc};
//
//        Set_Params(&params);
//
//        p -> pt_Dist = this;
//
//
//        F.params = &p;//double check with chen-yen on this
//        //Should I call the structure like I did with the function above?
//        //F.species = 1; can I do my species thing?
//        F.function = &gslClassWrapper;
//
//        T = gsl_root_fsolver_falsepos;
//        s = gsl_root_fsolver_alloc(T);
//        gsl_root_fsolver_set(s, &F, x_lo, x_hi);
//        //solving has been initialized
//
//        //Algorithm to find root
//        do
//        {
//            iter++;
//            status = gsl_root_fsolver_iterate(s);
//
//            mu.push_back(gsl_root_fsolver_root(s));
//
//            x_lo = gsl_root_fsolver_x_lower(s);
//            x_hi = gsl_root_fsolver_x_upper(s);
//            status = gsl_root_test_interval(x_lo, x_hi, 0, 0.0001);
//
//            if (status == GSL_SUCCESS)
//            {cout << "Converged" << endl;}
//        } while (status == GSL_CONTINUE && iter < max_iter);
//
//        gsl_root_fsolver_free(s); //delete memory
//        t += dt;
//    }
//
//    cout << "Chemical potential: \n";
//    for (int i = 0; i<=Tit; i++)
//    {
//        cout << mu[i] << endl;
//    }
//    
//}


