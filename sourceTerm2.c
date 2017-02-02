/* Soot Method of Moments with Interpolative Closure Model For ANSYS FLUENT 17.2

***********************************************
*   Anamol Pundle                             *
*   Department of Mechanical Engineering      *
*   University of Washington, Seattle         *
***********************************************

This code implements the Method of Moments with Interpolative Closure by Frenkalch et. al 
(as given in M.Frenklach, Method of Moments with Interpolative Closure, Chem. Eng. Sci.
Volume 57, Issue 12, June 2002, Pages 2229–2239) for use with ANSYS FLUENT 17.2. The 
implementation is for three moments through FLUENT UDF's. The aggregation model is not yet
implemented.

This code can currently only be used for turbulent non-premixed flames with the 
flamelet model. This file contains:
1. Source term UDF's for each moment.
    Moment-0 : nucleation, coagulation
    Moment-1 : nucleation, HACA surface growth, oxidation via O2 and OH
    Moment-2 : nucleation, coagulation, HACA surface growth, oxidation via O2 and OH

2. UDF for calculating effective turbulent diffusivity. 
 

RUNNING THE CODE FOR STEADY STATE SOLUTION:
1. Declare three user defined scalars(UDS) and x user defined memories in Fluent.
2. Compile the UDF's in this file.
3. Select only the nucleation source terms for each UDS.
4. Set the effective turbulent diffusivity of each UDS equal to the eff. turb. diff. UDF.
5. Solve only the three UDS equations until convergence (switch other equations off).
6. Add the HACA source term to the first and second moment UDS's and the coagulation source 
   term to the zeroth moment UDS. Solve until convergence, starting with low URF's.
7. Add the oxidation source terms to the first and second moment UDS's. Solve until 
   convergence.
8. Add the coagulation term to the second moment UDS. Solve until convergence.



USER DEFINED MEMORY ALLOCATION
UDM 0 : Nucleation source term for zeroth moment  
UDM 1 : Nucleation source term for first moment
UDM 2 : Nucleation source term for second moment 
UDM 3 : C2H2 surface growth term for first moment
UDM 4 : C2H2 surface growth term for second moment
UDM 5 : O2 oxidation term for first moment
UDM 6 : OH oxidation term for first moment
UDM 7 : O2 oxidation term for second moment
UDM 8 : OH oxidation term for second moment
*/


#include "udf.h"

/* universal constants and environmental conditions */
#define avogad 6.024e26 /* kmole^-1 */
#define kB 1.3806e-23 /* J/K */
#define R 8.314 /* kJ kmole^-1 Kelvin^-1 */
#define pi 3.14159
#define oneThird 1.0/3.0
#define oneHalf 1.0/2.0
#define twoThird 2.0/3.0
#define Patm 1.01e5 /* Pa */
#define ScT 0.7 /* Turbulent Schmidt Number */
#define PrT 0.85 /* Turbulent Prandtl Number */ 

/* molecular weights in kg/kmol */
#define C2H2MW 26.04 
#define H2MW 2.01 
#define HMW 1.007 
#define O2MW 32 
#define OHMW 13 
#define H2OMW 16 
#define pyrMW 202.25 
#define benzeneMW 78.11

/* masses and diameters of atoms */
#define mC 1.9944e-26 /* kg/atom */
#define mOH 2.8232e-26 /* kg/atom */
#define dAir 3.6e-10 /* m */
#define dC 2.2e-10 /* m */
#define dPyr 2.5*1.395e-10*sqrt(3) /* m */
#define dBenzene 1.395e-10 /* m */

/* number of C atoms in species */
#define NPyr 16
#define NBenzene 6

/* MOMIC parameters */
#define gammaPyr 0.025
#define gammaBenzene 0.001
#define vDW 2.2
#define siteDensity 2.3e19 /* m^2 */
#define rhoSoot 1800 /* kg/m^3 */
#define normParameter 1e15
#define a1 12.65
#define a2 -0.00563
#define b1 -1.38
#define b2 0.00068

/* MOMIC kinetic parameters SI units */
#define A1f 4.2e10
#define Ea1f 13000*4.18
#define A1r 3.9e9
#define Ea1r 11000*4.18
#define A2f 1.0e7
#define n2f 0.734
#define Ea2f 1430*4.18
#define A2r 3.68e5
#define n2r 1.139
#define Ea2r 1710*4.18
#define A3 2.0e10
#define A4 8.0e4
#define n4 1.56
#define Ea4 3800*4.18
#define A5 2.2e9
#define Ea5 7500*4.18
#define gammaOH 0.13

/* Fudge Factors */
#define fudgeFac1 0.8

real lagInterp3mom(real p, real m0, real m1, real m2); /* Lagragian interpolation */
real logInterp(real p, real a, real b, real M, real N); /* logarithmic interpolation */
real pyrNucSource(cell_t c,Thread *t);
real benzNucSource(cell_t c,Thread *t);
real chiSootCalc(cell_t c,Thread *t, real cellTemp, real cellRho);

DEFINE_SOURCE (m_0_NucSource,c,t,dS,eqn)
{

    /* Source for pyrene nucleation */
    real sourcePyr = pyrNucSource(c,t);
    /* Pyrene nucleation source calculation ends */

    /* Source for benzene nucleation */
    /*real sourceBenz = benzSource(c,t); */
    /* Benzene nucleation source ends */

    dS[eqn] = 0.0;
    C_UDMI(c,t,0) = sourcePyr;
    return sourcePyr;
}

DEFINE_SOURCE (m_1_NucSource,c,t,dS,eqn)
{

    /* Source for pyrene nucleation */
    real sourcePyr = (2*NPyr)*pyrNucSource(c,t);
    /* Pyrene nucleation source calculation ends */

    dS[eqn] = 0.0;
    C_UDMI(c,t,1) = sourcePyr;
    return sourcePyr;
}

DEFINE_SOURCE (m_2_NucSource,c,t,dS,eqn)
{

    /* Source for pyrene nucleation */
    real sourcePyr = (2*NPyr)*(2*NPyr)*pyrNucSource(c,t);
    /* Pyrene nucleation source calculation ends */

    dS[eqn] = 0.0;
    C_UDMI(c,t,2) = sourcePyr;
    return sourcePyr;
}

DEFINE_SOURCE (m_1_C2H2Source,c,t,dS,eqn)
{
    real source;
    real cellRho = C_R(c,t);
    real cellTemp = C_T(c,t);
    real cellM0 = C_UDSI(c,t,0);
    real cellM1 = C_UDSI(c,t,1);
    real cellM2 = C_UDSI(c,t,2);
    real a = a1 + a2*cellTemp;
    real b = b1 + b2*cellTemp;
    real alpha = tanh(a/log10(cellM1/cellM0)+b); 
    if (alpha < 0 ) { alpha = 0.01; }

    Material *mat = THREAD_MATERIAL(t);

    int idC2H2 = mixture_specie_index(mat, "c2h2");
    real C2H2Mf = Pdf_Yi(c,t,idC2H2);
    real C2H2Mc = C2H2Mf*cellRho/C2H2MW;
    
    real k4 = A4*pow(cellTemp, n4)*exp(-Ea4/(R*cellTemp));
    real chiSoot = chiSootCalc(c,t,cellTemp,cellRho);
    real Cs = pow((6*mC/(pi*rhoSoot)),oneThird);
    real MtwoThird = cellM0*pow(10, lagInterp3mom(twoThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));

    source = 2*k4*C2H2Mc*alpha*chiSoot*pi*pow(Cs,2)*MtwoThird;

    C_UDMI(c,t,3) = source;

    dS[eqn] = 0.0;
    return source;
}

DEFINE_SOURCE (m_2_C2H2Source,c,t,dS,eqn)
{
    real source;
    real cellRho = C_R(c,t);
    real cellTemp = C_T(c,t);
    real cellM0 = C_UDSI(c,t,0);
    real cellM1 = C_UDSI(c,t,1);
    real cellM2 = C_UDSI(c,t,2);
    real a = a1 + a2*cellTemp;
    real b = b1 + b2*cellTemp;
    real alpha = tanh(a/log10(cellM1/cellM0)+b); 
    if (alpha < 0 ) { alpha = 0.01; }

    Material *mat = THREAD_MATERIAL(t);

    int idC2H2 = mixture_specie_index(mat, "c2h2");
    real C2H2Mf = Pdf_Yi(c,t,idC2H2);
    real C2H2Mc = C2H2Mf*cellRho/C2H2MW;
    
    real k4 = A4*pow(cellTemp, n4)*exp(-Ea4/(R*cellTemp));
    real chiSoot = chiSootCalc(c,t,cellTemp,cellRho);
    real Cs = pow((6*mC/(pi*rhoSoot)),oneThird);

    real MtwoThird = cellM0*pow(10, lagInterp3mom(twoThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real MfiveThird = cellM0*pow(10, lagInterp3mom(5*oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));

    source = k4*C2H2Mc*alpha*chiSoot*pi*pow(Cs,2)*(4*MtwoThird + 4*MfiveThird);

    C_UDMI(c,t,4) = source;

    dS[eqn] = 0.0;
    return source;
}

DEFINE_SOURCE(m_1_OxSource,c,t,dS,eqn)
{
    real source;
    real sourceO2;
    real sourceOH;
    real cellRho = C_R(c,t);
    real cellTemp = C_T(c,t);
    real cellM0 = C_UDSI(c,t,0);
    real cellM1 = C_UDSI(c,t,1);
    real cellM2 = C_UDSI(c,t,2);
    real a = a1 + a2*cellTemp;
    real b = b1 + b2*cellTemp;
    real alpha = tanh(a/log10(cellM1/cellM0)+b); 
    if (alpha < 0 ) { alpha = 0.01; }

    Material *mat = THREAD_MATERIAL(t);

    int idO2 = mixture_specie_index(mat, "o2");
    int idOH = mixture_specie_index(mat, "oh");

    real O2Mf = Pdf_Yi(c,t,idO2);
    real OHMf = Pdf_Yi(c,t,idOH);

    real O2Mc = O2Mf*cellRho/O2MW;
    real OHMc = OHMf*cellRho/OHMW;

    real k5 = A5*exp(-Ea5/(R*cellTemp));
    real chiSoot = chiSootCalc(c,t,cellTemp,cellRho);
    real Cs = pow((6*mC/(pi*rhoSoot)),oneThird);

    real MtwoThird = cellM0*pow(10, lagInterp3mom(twoThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));

    sourceO2 = -2*k5*O2Mc*alpha*chiSoot*pi*pow(Cs,2)*MtwoThird;
    sourceOH = -2*gammaOH*OHMc*sqrt(pi*kB*cellTemp/(2*mOH))*pi*pow(Cs,2)*MtwoThird;
    source = sourceO2 + sourceOH;

    C_UDMI(c,t,5) = sourceO2;
    C_UDMI(c,t,6) = sourceOH;

    dS[eqn] = 0.0;
    return source*fudgeFac1;
}

DEFINE_SOURCE(m_2_OxSource,c,t,dS,eqn)
{
    real source;
    real sourceO2;
    real sourceOH;
    real cellRho = C_R(c,t);
    real cellTemp = C_T(c,t);
    real cellM0 = C_UDSI(c,t,0);
    real cellM1 = C_UDSI(c,t,1);
    real cellM2 = C_UDSI(c,t,2);
    real a = a1 + a2*cellTemp;
    real b = b1 + b2*cellTemp;
    real alpha = tanh(a/log10(cellM1/cellM0)+b); 
    if (alpha < 0 ) { alpha = 0.01; }

    Material *mat = THREAD_MATERIAL(t);

    int idO2 = mixture_specie_index(mat, "o2");
    int idOH = mixture_specie_index(mat, "oh");

    real O2Mf = Pdf_Yi(c,t,idO2);
    real OHMf = Pdf_Yi(c,t,idOH);

    real O2Mc = O2Mf*cellRho/O2MW;
    real OHMc = OHMf*cellRho/OHMW;

    real k5 = A5*exp(-Ea5/(R*cellTemp));
    real chiSoot = chiSootCalc(c,t,cellTemp,cellRho);
    real Cs = pow((6*mC/(pi*rhoSoot)),oneThird);

    real MtwoThird = cellM0*pow(10, lagInterp3mom(twoThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real MfiveThird = cellM0*pow(10, lagInterp3mom(5*oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));

    sourceO2 = k5*O2Mc*alpha*chiSoot*pi*pow(Cs,2)*(4*MtwoThird - 4*MfiveThird);
    sourceOH = gammaOH*OHMc*sqrt(pi*kB*cellTemp/(2*mOH))*pi*pow(Cs,2)*(4*MtwoThird - 4*MfiveThird);
    source = sourceO2 + sourceOH;

    C_UDMI(c,t,7) = sourceO2;
    C_UDMI(c,t,8) = sourceOH;

    dS[eqn] = 0.0;
    return source*fudgeFac1;
}

DEFINE_SOURCE (m_0_CoagSource,c,t,dS,eqn)
{
    real source;
    real cellPressure = C_P(c,t) + Patm;
    real cellTemp = C_T(c,t);
    real cellM0 = C_UDSI(c,t,0);
    real cellM1 = C_UDSI(c,t,1);
    real cellM2 = C_UDSI(c,t,2);
    real lamVisc = C_MU_L(c,t);

    real meanFreePath = kB*cellTemp/(sqrt(2)*pi*pow(dAir,2)*cellPressure);
    real avgDia = logInterp(oneThird, 0.0, 1.0, cellM0, cellM1)/cellM0*dC; 
    C_UDMI(c,t,9) = avgDia;
    real Kn = 2*meanFreePath/avgDia;
    C_UDMI(c,t,10) = Kn;
    real Kc = 2*kB*cellTemp/(3*lamVisc);
    real KcPrime = 2.514*meanFreePath*pow(pi*rhoSoot/6,oneThird);
    real Kf = vDW*sqrt(6*kB*cellTemp/rhoSoot)*pow(3/(4*pi*rhoSoot),(1/6));

    real muOneThird = logInterp(oneThird, 0.0, 1.0, (cellM0/cellM0), (cellM1/cellM0));
    real muTwoThird = logInterp(twoThird, 0.0, 1.0, (cellM0/cellM0), (cellM1/cellM0));
    real muFourThird = logInterp(2*twoThird, 1.0, 2.0, (cellM1/cellM0), (cellM2/cellM0));
    real muFiveThird = logInterp(5*oneThird, 1.0, 2.0, (cellM1/cellM0), (cellM2/cellM0));
    real muSevenThird = pow(10, lagInterp3mom(7*oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muEightThird = pow(10, lagInterp3mom(8*oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muMinusOneThird = pow(10, lagInterp3mom(-oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muMinusTwoThird = pow(10, lagInterp3mom(-twoThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));

    real f_0_0_0 = 2*muTwoThird + 2*muOneThird*muOneThird;
    real f_1_0_0 = 2*muFiveThird + 2*muTwoThird*cellM1/cellM0 + 4*muFourThird*muOneThird;
    real f_2_0_0 = 2*muEightThird + 2*muTwoThird*cellM2/cellM0 + 4*muFiveThird*cellM1/cellM0 + 4*muSevenThird*muOneThird + 4*pow(muFourThird,2);
    real f_oneHalf_0_0 = pow(10, lagInterp3mom(oneHalf, log10(f_0_0_0), log10(f_1_0_0), log10(f_2_0_0)));

    real Gc = -Kc*(1 + muOneThird*muMinusOneThird + KcPrime*(muMinusOneThird + muOneThird*muMinusTwoThird))*pow(cellM0,2);

    real Gf = -0.5*Kf*pow(cellM0,2)*f_oneHalf_0_0;

    C_UDMI(c,t,11) = Gc;
    C_UDMI(c,t,12) = Gf;


    if (Kn < 0.1) { source = Gc; }

    else if (Kn > 8.0) { source = Gf; }

    else { source = Gc*Gf/(Gc + Gf); }

    dS[eqn] = 0.0;

    C_UDMI(c,t,13) = source;
    return source;

}

DEFINE_SOURCE (m_2_CoagSource,c,t,dS,eqn)
{

    real source;
    real cellPressure = C_P(c,t) + Patm;
    real cellTemp = C_T(c,t);
    real cellM0 = C_UDSI(c,t,0);
    real cellM1 = C_UDSI(c,t,1);
    real cellM2 = C_UDSI(c,t,2);
    real lamVisc = C_MU_L(c,t);

    real meanFreePath = kB*cellTemp/(sqrt(2)*pi*pow(dAir,2)*cellPressure);
    real avgDia = logInterp(oneThird, 0.0, 1.0, cellM0, cellM1)/cellM0*dC; 
    real Kn = 2*meanFreePath/avgDia;
    real Kc = 2*kB*cellTemp/(3*lamVisc);
    real KcPrime = 2.514*meanFreePath*pow(pi*rhoSoot/6,oneThird);
    real Kf = vDW*sqrt(6*kB*cellTemp/rhoSoot)*pow(3/(4*pi*rhoSoot),(1/6));

    real muOneThird = logInterp(oneThird, 0.0, 1.0, (cellM0/cellM0), (cellM1/cellM0));
    real muTwoThird = logInterp(twoThird, 0.0, 1.0, (cellM0/cellM0), (cellM1/cellM0));
    real muFourThird = logInterp(2*twoThird, 1.0, 2.0, (cellM1/cellM0), (cellM2/cellM0));
    real muFiveThird = logInterp(5*oneThird, 1.0, 2.0, (cellM1/cellM0), (cellM2/cellM0));
    real muSevenThird = pow(10, lagInterp3mom(7*oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muEightThird = pow(10, lagInterp3mom(8*oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));

    real f_0_1_1 = 2*muFiveThird*cellM1/cellM0 + 2*muFourThird*muFourThird;
    real f_1_1_1 = 2*muEightThird*cellM1/cellM0 + 2*muTwoThird*cellM1/cellM0 + 4*muSevenThird*muFourThird;
    real f_oneHalf_1_1 = logInterp(oneHalf,0.0,1.0,f_0_1_1,f_1_1_1);

    real Gc = Kc*(2*cellM1*cellM1/(cellM0*cellM0) + 2*muFourThird*muTwoThird + KcPrime*(2*muTwoThird*muOneThird + 2*muFourThird*muOneThird))*pow(cellM0,2);

    real Gf = Kf*pow(cellM0,2)*f_oneHalf_1_1;

    C_UDMI(c,t,14) = Gc;
    C_UDMI(c,t,15) = Gf;

    if (Kn < 0.1) { source = Gc; }

    else if (Kn > 8.0) { source = Gf; }

    else { source = Gc*Gf/(Gc + Gf); }

    dS[eqn] = 0.0;
    C_UDMI(c,t,16) = source;

    return source;
}


real logInterp(real p, real a, real b, real M, real N)
{
    real f = (p-a)/(b-a);
    return pow(N,f) * pow(M,1-f);
}

real lagInterp3mom(real p, real m0, real m1, real m2)
{
    real mp = 0.5*(p-1)*(p-2)*m0 - p*(p-2)*m1 + 0.5*p*(p-1)*m2;

    return mp;

}

DEFINE_DIFFUSIVITY (uds_diff,c,t,i)
{
    real diff = C_MU_EFF(c,t)/(ScT*PrT); 

    return diff;
}

real pyrNucSource(cell_t c,Thread *t)
{
    Material *mat = THREAD_MATERIAL(t); 
    int idPyr = mixture_specie_index(mat, "a4"); 
    real pyrMf = Pdf_Yi(c,t,idPyr); 
    real pyrMc = pyrMf*C_R(c,t)/pyrMW; 
    real CNucPyr = vDW*sqrt(4*pi*kB/(mC*NPyr))*pow(dPyr*avogad,2);
    real sourcePyrNuc = gammaPyr*CNucPyr*sqrt(C_T(c,t))*pyrMc*pyrMc; 
    sourcePyrNuc = sourcePyrNuc/normParameter;
    return sourcePyrNuc;
}

real benzNucSource(cell_t c,Thread *t)
{
    Material *mat = THREAD_MATERIAL(t);
    int idBenzene = mixture_specie_index(mat, "c6h6");
    real benzeneMf = Pdf_Yi(c,t,idBenzene);
    real benzeneMc = benzeneMf*C_R(c,t)/benzeneMW;
    real CNucBenzene = vDW*sqrt(4*pi*kB/(mC*NBenzene))*pow(dBenzene*avogad,2);
    real sourceBenzeneNuc = gammaBenzene*CNucBenzene*sqrt(C_T(c,t))*benzeneMc*benzeneMc;
    sourceBenzeneNuc = sourceBenzeneNuc/normParameter;
    return sourceBenzeneNuc;
}

real chiSootCalc(cell_t c,Thread *t, real cellTemp, real cellRho)
{

    Material *mat = THREAD_MATERIAL(t);

    int idC2H2 = mixture_specie_index(mat, "c2h2");
    int idH2 = mixture_specie_index(mat, "h2");
    int idH = mixture_specie_index(mat, "h");
    int idO2 = mixture_specie_index(mat, "o2");
    int idOH = mixture_specie_index(mat, "oh");
    int idH2O = mixture_specie_index(mat, "h2o");

    real C2H2Mf = Pdf_Yi(c,t,idC2H2);
    real H2Mf = Pdf_Yi(c,t,idH2);
    real HMf = Pdf_Yi(c,t,idH);
    real O2Mf = Pdf_Yi(c,t,idO2);
    real OHMf = Pdf_Yi(c,t,idOH);
    real H2OMf = Pdf_Yi(c,t,idH2O);

    real C2H2Mc = C2H2Mf*cellRho/C2H2MW;
    real H2Mc = H2Mf*cellRho/H2MW;
    real HMc = HMf*cellRho/HMW;
    real O2Mc = O2Mf*cellRho/O2MW;
    real OHMc = OHMf*cellRho/OHMW;
    real H2OMc = H2OMf*cellRho/H2OMW;
    

    real k1f = A1f*exp(-Ea1f/(R*cellTemp));
    real k1r = A1r*exp(-Ea1r/(R*cellTemp));
    real k2f = A2f*pow(cellTemp, n2f)*exp(-Ea2f/(R*cellTemp));
    real k2r = A2r*pow(cellTemp, n2r)*exp(-Ea2r/(R*cellTemp));
    real k3 = A3;
    real k4 = A4*pow(cellTemp, n4)*exp(-Ea4/(R*cellTemp));
    real k5 = A5*exp(-Ea5/(R*cellTemp));
    real chiSoot = (k1f*HMc + k2f*OHMc)/(k1r*H2Mc + k2r*H2OMc + k3*HMc + k4*C2H2Mc + k5*O2Mc) * siteDensity;
    return chiSoot;
}


