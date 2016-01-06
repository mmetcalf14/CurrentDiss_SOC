//
//  Distribution.cpp
//  CurrentDissipation_1DRing
//
//  Created by mekena McGrew on 1/6/16.
//  Copyright © 2016 Mekena Metcalf. All rights reserved.
//

#include "Distribution.hpp"
using namespace std;

Distribution::Distribution(int _Nl, double _N, double _T, double _E, double _w, double _dt, double _To)
{
    Nl = _Nl;
    N = _N;
    Temp = _T;
    Eo = _E;
    w_soc = _w;
    dt = _dt;
    To = _To;
    Tit = To/dt;
    Nt = (2*Nl)+1;
    
}


double Distribution::TotalNumber(double x, void *par)//this species thing won't work will it?
{
    TotalNumber_params *p = (TotalNumber_params*) par;//would I be better off making this a
    //full class?
    //point to stucture variables with pointer *p. Keeps writing simple
    //how to call structure in function
    
    double N = p->N;
    double T = p->T;
    double Phi = p->Phi;
    double Eo = p->Eo;
    double w = p->w;
    int Nl = p->Nl;
    
    int Ntable = (2*Nl)+1;
    double TotalSum = 0.0;
    double NumberFunction;
    double E_up = 0.0;
    double E_dn = 0.0;
    //double E = 0.0;
    //THis creates the variables in the function which point to definitions
    //within the structure
    //cout << endl;
    
    for (int i = 0; i < Ntable; i++)
    {
        int n = i-Nl;
        //        if(species == 0)// REGULAR CURRENT
        //        {
        //            E = Eo*(n-Phi)*(n-Phi);
        //        }
        //        if (species == 1)// RASHBA SOC
        //        {
        E_up = Eo*((n - Phi + (1./2.)*(1-sqrt(1 + (w*w))))*(n - Phi + (1./2.)*(1-sqrt(1 + (w*w)))));
        TotalSum += 1 / (exp((E_up - x) / T) + 1);
        
        E_dn = Eo*((n - Phi + (1./2.)*(1+sqrt(1 + (w*w))))*(n - Phi + (1./2.)*(1+sqrt(1 + (w*w)))));
        TotalSum += 1 / (exp((E_dn - x) / T) + 1);
        //        }
        
    }
    
    NumberFunction = TotalSum - (2.*N); //how is this function infinite?
    //cout << Phi << " " << x << " " << TotalSum << " " << N << " " << NumberFunction << endl;
    return NumberFunction;
    
    
}

void Distribution::ChemicalPotential()
{
    
    double t = 0.;
    double Phi = 0;
    
    
    for(int i = 0; i <= Tit; i++)
    {
        Phi = t/(2*To);
        int status;
        int iter = 0, max_iter = 1000;
        const gsl_root_fsolver_type *T;
        gsl_root_fsolver *s;
        double x_lo = -2000.0;
        double x_hi = 3000.0;
        
        gsl_function F;
        struct TotalNumber_params params = { Nl, N, Temp, Phi, Eo, w_soc};
        
        Set_Params(&params);
        
        p -> pt_Dist = this;
        
        
        F.params = &p;//double check with chen-yen on this
        //Should I call the structure like I did with the function above?
        //F.species = 1; can I do my species thing?
        F.function = &gslClassWrapper;
        
        T = gsl_root_fsolver_falsepos;
        s = gsl_root_fsolver_alloc(T);
        gsl_root_fsolver_set(s, &F, x_lo, x_hi);
        //solving has been initialized
        
        //Algorithm to find root
        do
        {
            iter++;
            status = gsl_root_fsolver_iterate(s);
            
            mu.push_back(gsl_root_fsolver_root(s));
            
            x_lo = gsl_root_fsolver_x_lower(s);
            x_hi = gsl_root_fsolver_x_upper(s);
            status = gsl_root_test_interval(x_lo, x_hi, 0, 0.0001);
            
            if (status == GSL_SUCCESS)
            {cout << "Converged" << endl;}
        } while (status == GSL_CONTINUE && iter < max_iter);
        
        gsl_root_fsolver_free(s); //delete memory
        t += dt;
    }
    
    cout << "Chemical potential: \n";
    for (int i = 0; i<=Tit; i++)
    {
        cout << mu[i] << endl;
    }
    
}

void Distribution::FD_Distribution(int time_rev, int species, std::ofstream &fout, double **FD_Mat)
{
    double t = 0.0;
    double phi = 0.0;
    double E = 0.0;
    double FermiFunc;
    
    
    
    for (int i = 0; i <= Tit; i++)//do I need two arrays? one for up and one for down
    {
        if(time_rev == 0)
        {
            phi = t/(2*To);
        }
        else if (time_rev == 1)//this solve the problem of SOC. Will need for FD_Mat in main
        {
            phi = (-t)/(2*To);
        }
        for (int j = 0; j < Nt; j++)
        {int n = j - Nl;
            
            if (species == 0)//reg energy
            {
                E = Eo*(n-phi)*(n-phi);
                
            }
            else if (species == 1)//spin up energy
            {
                E = Eo*((n - phi + (1./2.)*(1-sqrt(1 + (w_soc*w_soc))))*(n - phi + (1./2.)*(1-sqrt(1 + (w_soc*w_soc)))));
            }
            else if  (species == 2)//spin down energy
            {
                
                E = Eo*((n - phi + (1./2.)*(1+sqrt(1 + (w_soc*w_soc))))*(n - phi + (1./2.)*(1+sqrt(1 + (w_soc*w_soc)))));
            }
            
            FermiFunc = 1/(exp((E-mu[i])/Temp)+1);
            
            FD_Mat[i][j] = FermiFunc;
            
        }
        
        t+= dt;
        
    }
    
}

