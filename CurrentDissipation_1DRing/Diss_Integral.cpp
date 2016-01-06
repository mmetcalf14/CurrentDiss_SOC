//
//  Diss_Integral.cpp
//  CurrentDissipation_1DRing
//
//  Created by mekena McGrew on 1/6/16.
//  Copyright © 2016 Mekena Metcalf. All rights reserved.
//

#include "Diss_Integral.hpp"

double ** DissipatedIntegral::DissIntegral(const Distribution &FD, double ** FD_Pos, double ** FD_Neg)
{
    Tint = FD.Tit/Nint;
    int Tdim = (3*Tint)+1;//this could be an input

    double num1 = (14./45.);
    double num2 = (64./45.);
    double num3 = (24./45.);

    double sum1;
    double sum2;
    double sum3;
    double sum4;
    double sum5;
    double sum6;

    double t=0.;
    double Tf=0.;
    dt = FD.dt;

    double F_bin[Nint+1];
    double **Fermi_Func; //I dont know what to do about this array
    Fermi_Func = new double*[Tdim];
    for (int i = 0; i < Tdim; i++)
        Fermi_Func[i] = new double[FD.Nt];


    std::cout << "Beginning big loop\n";
    //This loop was constructed this way to remove a weird time dependence which in turn can eliminate error
    //prior loop construction is commented out below
    for ( int i = 0; i <= (3*Tint); i++)
    {
        int ip = Nint*i;
        Tf = ip*FD.dt;



        for ( int j =0; j < FD.Nt; j++)
        {
            //n = j - Nl; //energy levels
            int m = 0;
            int a = 125; //read table from bottom row
//            int b = (Ntable-1) - j; //read table from right column


            sum1 = 0;
            sum2 = 0;
            sum3 = 0;
            sum4 = 0;
            sum5 = 0;
            sum6 = 0;

            if ( i <= (Tint/2))//should this be less than and not leq?
            {
                for(int k = 0; k < i; k++)
                {int kp = Nint*k;
                    t = kp*FD.dt;

                    for(int l = 0; l <= Nint; l++)
                    {
                        F_bin[l] = FD_Pos[kp+l][j];
                    }
                    sum1 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));

                }

                Fermi_Func[i][j] = FD_Pos[0][j]*exp(-Tf/Tau) + sum1;
                //The array Fermi_Func is under a new time incriment of .005 rather than .001
                //Therefore in the current we take points every .005 which is the incriment of Tf
            }
            else if( i > (Tint/2) && i <= Tint)
            {
                for(int k = 0; k < i; k++)
                {
                    int kp = Nint*k;
                    if ( k == (Tint/2))
                    {
                        a = 125;
                    }
                    int ap = a*Nint;
                    t = kp*dt;

                    if ( k < (Tint/2))
                    {

                        for(int l = 0; l < (Nint+1); l++)
                        {
                            F_bin[l] = FD_Pos[kp+l][j];

                        }
                        sum1 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));
                    }
                    else
                    {  for(int l = 0; l < (Nint+1); l++)
                    {
                        F_bin[l] = FD_Pos[ap-l][j];
                    }
                        sum2 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));
                    }

                    //                    if (j==0)
                    //                    {
                    //                        cout << Tf << " " << t << " "  << sum1 << " " <<sum2 << endl;
                    //                    }

                    a--;
                }

                Fermi_Func[i][j] = FD_Pos[0][j]*exp(-Tf/Tau) + sum1 + sum2;
            }
            else if (i > Tint && i <= (3*Tint)/2)
            {
                for(int k = 0; k < i; k++)//k<= i -> k < i
                {
                    if ( k == (Tint/2) || k == (Tint ))
                    {
                        a = 125;//changed from 99 to 100
                    }
                    if ( k == (Tint/2) || k == (Tint ))
                    {
                        m = 0;//changed from 1 to 0 to account for first group
                        //check if this matches up
                    }

                    int kp = Nint*k;
                    int ap = Nint*a;
                    int mp = Nint*m;
                    t = kp*dt;


                    if ( k < (Tint/2))//kp, j
                    {
                        for(int l = 0; l < (Nint+1); l++)
                        {
                            F_bin[l] = FD_Pos[kp+l][j];

                        }

                        sum1 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));
                    }
                    else if (k >= (Tint/2) && k < Tint)//ap, j
                    {
                        for(int l = 0; l < (Nint+1); l++)
                        {
                            F_bin[l] = FD_Pos[ap-l][j];

                        }
                        sum2 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));
                    }
                    else//mp, b
                    { for(int l = 0; l < (Nint+1); l++)
                    {
                        F_bin[l] = FD_Neg[mp+l][j];

                    }
                        sum3 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));
                    }

                    a--;
                    m++;


                }

                Fermi_Func[i][j] = FD_Pos[0][j]*exp(-Tf/Tau) + sum1 + sum2 + sum3;
            }
            else if (i > (3*Tint)/2 && i <= (2*Tint))
            { for(int k = 0; k < i; k++)
            {
                if ( k == (Tint/2) || k == (Tint ) || k == ((3*Tint)/2) || k == ((2*Tint)))
                {
                    a = 125;
                }
                if ( k == (Tint/2) || k == (Tint ) || k == ((3*Tint)/2) || k == ((2*Tint)))
                {
                    m = 0;
                }
                int kp = Nint*k;
                int ap = Nint*a;
                int mp = Nint*m;
                t = kp*dt;

                if ( k < (Tint/2))
                {
                    for(int l = 0; l < (Nint+1); l++)
                    {
                        F_bin[l] = FD_Pos[kp+l][j];
                    }
                    sum1 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));
                }
                else if (k >= (Tint/2) && k < Tint)
                {
                    for(int l = 0; l < (Nint+1); l++)
                    {
                        F_bin[l] = FD_Pos[ap-l][j];
                    }
                    sum2 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));

                }
                else if (k >= Tint && k <(3*Tint)/2)
                { for(int l = 0; l < (Nint+1); l++)
                {
                    F_bin[l] = FD_Neg[mp+l][j];
                }
                    sum3 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));
                }
                else //ap, b
                {
                    for(int l = 0; l < (Nint+1); l++)
                    {
                        F_bin[l] = FD_Neg[ap-l][j];
                    }
                    sum4 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));
                }
                a--;
                m++;

            }

                Fermi_Func[i][j] = FD_Pos[0][j]*exp(-Tf/Tau) + sum1 + sum2 + sum3 + sum4;
            }

            else if ( i > (2*Tint) && i <= (5*Tint)/2)
            {
                for(int k = 0; k < i; k++)
                {
                    if ( k == (Tint/2) || k == (Tint) || k == ((3*Tint)/2) || k == ((2*Tint))|| k == ((5*Tint)/2))
                    {
                        a = 125;
                    }
                    if ( k == (Tint/2) || k == (Tint) || k == ((3*Tint)/2) || k == ((2*Tint))|| k == ((5*Tint)/2) )
                    {
                        m = 0;
                    }
                    int kp = Nint*k;
                    int ap = Nint*a;
                    int mp = Nint*m;
                    t = kp*dt;

                    if ( k < (Tint/2))
                    {
                        for(int l = 0; l < (Nint+1); l++)
                        {
                            F_bin[l] = FD_Pos[kp+l][j];
                        }
                        sum1 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));
                    }
                    else if (k >= (Tint/2) && k < Tint)
                    {    for(int l = 0; l < (Nint+1); l++)
                    {
                        F_bin[l] = FD_Pos[ap-l][j];
                    }
                        sum2 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));

                    }
                    else if (k >= Tint && k <(3*Tint)/2)
                    {
                        for(int l = 0; l < (Nint+1); l++)
                        {
                            F_bin[l] = FD_Neg[mp+l][j];
                        }
                        sum3 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));
                    }
                    else if (k >=(3*Tint)/2 && k < (2*Tint))
                    {
                        for(int l = 0; l < (Nint+1); l++)
                        {
                            F_bin[l] = FD_Neg[ap-l][j];
                        }
                        sum4 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));
                    }
                    else //mp, j
                    {
                        for(int l = 0; l < (Nint+1); l++)
                        {
                            F_bin[l] = FD_Pos[mp+l][j];
                        }
                        sum5 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));
                    }

                    a--;
                    m++;

                }

                Fermi_Func[i][j] = FD_Pos[0][j]*exp(-Tf/Tau) + sum1 + sum2 + sum3 + sum4 + sum5;
            }
            else
            {
                for(int k = 0; k < i; k++)
                {
                    if ( k == (Tint/2) || k == (Tint) || k == ((3*Tint)/2) || k == ((2*Tint))|| k == ((5*Tint)/2))
                    {
                        a = 125;
                    }
                    if ( k == (Tint/2) || k == (Tint) || k == ((3*Tint)/2) || k == ((2*Tint)) || k == ((5*Tint)/2))
                    {
                        m = 0;
                    }

                    int kp = Nint*k;
                    int ap = Nint*a;
                    int mp = Nint*m;
                    t = kp*dt;
                    if ( k < (Tint/2))
                    {
                        for(int l = 0; l < (Nint+1); l++)
                        {
                            F_bin[l] = FD_Pos[kp+l][j];

                        }
                        sum1 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));
                    }
                    else if (k >= (Tint/2) && k < Tint)
                    {
                        for(int l = 0; l < (Nint+1); l++)
                        {
                            F_bin[l] = FD_Pos[ap-l][j];

                        }
                        sum2 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));

                    }
                    else if (k >= Tint && k <(3*Tint)/2)
                    {
                        for(int l = 0; l < (Nint+1); l++)
                        {
                            F_bin[l] = FD_Neg[mp+l][j];

                        }
                        sum3 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));
                    }
                    else if (k >= (3*Tint)/2 && k < (2*Tint))
                    {
                        for(int l = 0; l < (Nint+1); l++)
                        {
                            F_bin[l] = FD_Neg[ap-l][j];

                        }
                        sum4 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));
                    }
                    else if (k >= (2*Tint) && k < (5*Tint)/2)
                    {
                        for(int l = 0; l < (Nint+1); l++)
                        {
                            F_bin[l] = FD_Pos[mp+l][j];

                        }
                        sum5 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));
                    }
                    else // ap, j
                    {
                        for(int l = 0; l < (Nint+1); l++)
                        {
                            F_bin[l] = FD_Pos[ap-l][j];

                        }
                        sum6 += dt*((num1*((1./Tau)*exp(-(Tf-(t))/Tau))*F_bin[0])+(num2*((1./Tau)*exp(-(Tf-(t+dt))/Tau))*F_bin[1])+(num3*((1./Tau)*exp(-(Tf-(t+(2.*dt)))/Tau))*F_bin[2])+(num2*((1./Tau)*exp(-(Tf-(t+(3.*dt)))/Tau))*F_bin[3])+(num1*((1./Tau)*exp(-(Tf-(t+(4.*dt)))/Tau))*F_bin[4]));
                    }
                    a--;
                    m++;

                }
                //            if (j==0)
                //            {cout << Tf << " " << sum1 + sum2 + sum3+ sum4 + sum5 + sum6<<endl;
                //            }
                Fermi_Func[i][j] = FD_Pos[0][j]*exp(-Tf/Tau) + sum1 + sum2 + sum3 + sum4 + sum5 + sum6;
            }
        }
        //cout << endl;

    }

    return Fermi_Func;
}

void DissipatedIntegral::DissipatedCurrent(const Distribution &FD, std::ofstream& PC_Diss_out, double jo, double ** FDu_Diss, double ** FDd_Diss)
{

    double t= 0.;
    //double dt = FD.dt;
    double Phi;
    double Tot_Diss_Current;
    double jn_up;
    double jn_dn;



    for (int i =0; i <= (3*Tint); i++) //determine total current for each phi value after dissapation
    {
        int ip = Nint * i;
        t = ip*dt;
        //cout << t << endl;
        if (i <= (Tint/2))
            Phi = t/(2.*FD.To);
        else if (i > (Tint/2) && i <=((3*Tint)/2))
            Phi = 0.25 -((t-(FD.To/2.))/(2.*FD.To));
        else if ( i >((3*Tint)/2) && i <= ((5*Tint)/2) )
            Phi = -0.25 + ((t-((3*FD.To)/2.))/(2.*FD.To));
        else
            Phi = 0.25 -((t-((5*FD.To)/2.))/(2.*FD.To));

        //Phi_Diss_Table[i] = Phi;

        Tot_Diss_Current = 0;

        for (int j = 0; j < FD.Nt; j++)
        {
            int n = j - FD.Nl;

            jn_up = jo*(n - Phi + (1./2.)*(1-sqrt(1 + (FD.w_soc*FD.w_soc))));
            jn_dn = jo*(n - Phi + (1./2.)*(1+sqrt(1 + (FD.w_soc*FD.w_soc))));
            Tot_Diss_Current += jn_up * FDu_Diss[i][j];
            Tot_Diss_Current += jn_dn * FDd_Diss[i][j];

        }
        PC_Diss_Table.push_back(Tot_Diss_Current);


        PC_Diss_out << Phi << " " << Tot_Diss_Current << std::endl;

    }

    CurveArea();

}

void DissipatedIntegral::CurveArea()
{
    double area1 = 0.;
    double area2 = 0.;
    double area = 0.;

    for(int i = ((Tint/2)-1); i < (5*Tint)/2; i++)
    {

        if(i <= (3*Tint)/2)
        {
            area1 += Nint*dt*PC_Diss_Table[i];

        }
        else
        {
            area2 += Nint*dt*PC_Diss_Table[i];

        }
    }

    area = area1 - area2;
    std::cout << "This is the area of the curve: "<< area << std::endl;
}


