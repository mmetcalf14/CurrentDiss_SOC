//
//  EQ_Current.cpp
//  CurrentDissipation_1DRing
//
//  Created by mekena McGrew on 1/6/16.
//  Copyright © 2016 Mekena Metcalf. All rights reserved.
//

#include "EQ_Current.hpp"
void EQ_Current::Full_EQ_PC(const Distribution &FD,  double jo, double ** FDu_Mat, double ** FDd_Mat, double ** FDu_Neg, double ** FDd_Neg)
{

    double Phi = 0.;
    double t = 0.;
    int m = 0;

    for (int i = 0; i <= (2*FD.Tit); i++) //determine current in full equilibrium
        //changed from Tm to 2Tm to incorporate negative phi
    {

        if ( i == (FD.Tit/2)+1 || i == (FD.Tit +1) || i == ((3*FD.Tit)/2)+1)
        {
            m = 1;
        }
        if (i <= (FD.Tit/2))
            Phi = t/(2.*FD.To);
        else if (i > (FD.Tit/2) && i <=((3*FD.Tit)/2))
            Phi = 0.25 -((t-(FD.To/2.))/(2.*FD.To));
        else
            Phi = -0.25 + ((t-((3*FD.To)/2.))/(2.*FD.To));

        Phi_Table.push_back(Phi);

        Tot_Current = 0.;
        // cout << Phi << " " << m <<endl;
        int n;
        int k = (FD.Tit/2) - m;// reads row back
        //cout << m << " " << k << endl;

        for ( int j = 0; j < FD.Nt; j++)
        {

//            //read table from bottom row
//            int l = (FD.Nt-1) - j; //read table from right column
            n = j - FD.Nl;
            jn_up = jo*(n - Phi + (1./2.)*(1-sqrt(1 + (FD.w_soc*FD.w_soc))));
            jn_dn = jo*(n - Phi + (1./2.)*(1+sqrt(1 + (FD.w_soc*FD.w_soc))));

            if (i <= (FD.Tit/2))
            {
                Tot_Current += jn_up * FDu_Mat[m][j];//m and k stay the same
                Tot_Current += jn_dn * FDd_Mat[m][j];

            }
            else if (i > (FD.Tit/2) && i <= FD.Tit)
            {
                Tot_Current += jn_up * FDu_Mat[k][j];//reading time back is fine
                Tot_Current += jn_dn * FDd_Mat[k][j];

            }
            else if (i > FD.Tit && i <= ((3*FD.Tit)/2))//following loops are for negative phi since
                //f_n(\phi) \neq f_{-n}(\phi)
            {

                Tot_Current += jn_up * FDu_Neg[m][j];
                Tot_Current += jn_dn * FDd_Neg[m][j];

            }
            else
            {
                Tot_Current += jn_up * FDu_Neg[k][j];
                Tot_Current += jn_dn * FDd_Neg[k][j];
            }
            //cout << "Last Part \n";
            //cout << Phi << " " << Total_Current_Eq <<endl;

        }

        PC_EQ_Table.push_back(Tot_Current);

        std::cout << Phi_Table[i] <<" "<< PC_EQ_Table[i] <<std::endl;


        t += FD.dt;
        m++;

    }
//end loop

}

void EQ_Current::Full_NEQ_PC(const Distribution &FD,  double jo, double ** FDu_Mat, double ** FDd_Mat, double ** FDu_Neg, double ** FDd_Neg)
{
    double t = -FD.To;
    double Phi = 0;

    for(int i = 0; i <= (2*FD.Tit); i++) //determine current without dissapation
    {
        Phi = t/(2.*FD.To); //added to span whole spectrum
        Tot_Current = 0;
        int n;
        for ( int j = 0; j < FD.Nt; j++)
        {
            n = j - FD.Nl;
            //jn = jo*(n-Phi_Table[i]);
            jn_up = jo*(n - Phi + (1./2.)*(1-sqrt(1 + (FD.w_soc*FD.w_soc))));
            jn_dn = jo*(n - Phi + (1./2.)*(1+sqrt(1 + (FD.w_soc*FD.w_soc))));


            if (i <= (FD.Tit/2))
            {
                Tot_Current += jn_up * FDu_Mat[0][j];
                Tot_Current += jn_dn * FDd_Mat[0][j];

            }
            else if (i > (FD.Tit/2) && i <= FD.Tit)
            {
                Tot_Current += jn_up * FDu_Mat[0][j];
                Tot_Current += jn_dn * FDd_Mat[0][j];

            }
            else if (i > FD.Tit && i <= ((3*FD.Tit)/2))
            {
                Tot_Current += jn_up * FDu_Neg[0][j]; //problem is array isn't in i dim I need a new indices
                Tot_Current += jn_dn * FDd_Neg[0][j];
                //cout << n << " " << Fn_Array[i][l] <<endl;
            }
            else
            {
                Tot_Current += jn_up * FDu_Neg[0][j];
                Tot_Current += jn_dn * FDd_Neg[0][j];
            }



        }
        PC_NEQ_Table.push_back(Tot_Current);

        t += FD.dt;

    }
    //end loop
}

void EQ_Current::CurrentOut(const Distribution &FD, std::ofstream &fout)
{
    for ( int i = 0; i <= (2*FD.Tit); i++)
    {
        fout << Phi_Table[i] << " "  << PC_EQ_Table[i] << " " << PC_NEQ_Table[i] << std::endl;
    }
}



