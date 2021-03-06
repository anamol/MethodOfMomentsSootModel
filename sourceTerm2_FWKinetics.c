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
    Moment-1 : nucleation, surface growth via HACA and PAH condensation, oxidation via O2 and OH
    Moment-2 : nucleation, coagulation, surface growth via HACA and PAH condensation, oxidation via O2 and OH

2. UDF for calculating effective turbulent diffusivity. 

3. UDF's for calculating the absorption coefficient of sooty gas, taken from Widmann (2003)
    and Sazhin (1994). 
 

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

Source terms:
UDM 0 : Nucleation source term for zeroth moment  
UDM 1 : Nucleation source term for first moment
UDM 2 : Nucleation source term for second moment 
UDM 3 : C2H2 surface growth term for first moment
UDM 4 : C2H2 surface growth term for second moment
UDM 5 : O2 oxidation term for first moment
UDM 6 : OH oxidation term for first moment
UDM 7 : Total oxidation term for first moment
UDM 8 : O2 oxidation term for second moment
UDM 9 : OH oxidation term for second moment
UDM 10: Total oxidation term for second moment
UDM 11: PAH surface growth source term for first moment
UDM 12: PAH surface growth source term for second moment
UDM 13: Continuum coagulation source term for zeroth moment
UDM 14: Free molecular coagulation source term for zeroth moment
UDM 15: Actual coagulation source term for zeroth moment 
UDM 16: Continuum coagulation source term for second moment
UDM 17: Free molecular coagulation source term for second moment
UDM 18: Actual coagulation source term for second moment 

Others:
UDM 19: alpha (fraction of surface sites available for reaction)
UDM 20: Mean free path of gas
UDM 21: Average diameter of soot particle
UDM 22: Knudsen number
*/


#include "udf.h"

/* universal constants and environmental conditions */
#define avogad 6.024e26 /* kmole^-1 */
#define kB 1.3806e-23 /* J/K */
#define R 8.314 /* kJ kmole^-1 Kelvin^-1 */
#define pi 3.14159
#define oneThird 1.0/3.0
#define oneHalf 1.0/2.0
#define oneSixth 1.0/6.0
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
#define naphMW 128.17
#define phenanMW 178.234

/* masses and diameters of atoms */
#define mC 1.9944e-26 /* kg/atom */
#define mOH 2.8232e-26 /* kg/atom */
#define dAir 3.6e-10 /* m */
#define dC 2.2e-10 /* m */
#define dA 1.395*sqrt(3)
#define dPyr 2.5*1.395e-10*sqrt(3) /* m */
#define dBenzene 1.395e-10 /* m */

/* number of C atoms in species */
#define NPyr 16
#define NBenzene 6
#define NPhenan 14
#define NNaph 10

/* MOMIC parameters */
#define gammaPyr 0.001
#define gammaBenzene 0.001
#define vDW 4.0
#define siteDensity 2.3e19 /* m^2 */
#define rhoSoot 1800 /* kg/m^3 */
#define normParameter 1e15
#define a1 12.65
#define a2 -56.3e-4
#define b1 -1.38
#define b2 6.8e-4

/* MOMIC kinetic parameters SI units */
#define A1f 2.5e11
#define Ea1f 66900
#define A1r 3.9e9
#define Ea1r 39000
#define A2 1.0e11
#define A3 8.4e8
#define n3 0.4
#define Ea3 35100
#define A4 2.2e9
#define Ea4 31300
#define gammaOH 0.13

/* Fudge Factors */
#define fudgeFac1 1.0
#define fudgeFac2 1.0

real lagInterp3mom(real p, real m0, real m1, real m2); /* Lagragian interpolation */
real logInterp(real p, real a, real b, real M, real N); /* logarithmic interpolation */
real pyrNucSource(cell_t c,Thread *t);
real benzNucSource(cell_t c,Thread *t);
real chiSootCalc(cell_t c,Thread *t, real cellTemp, real cellRho);

DEFINE_SOURCE (m_0_NucSourcePyr,c,t,dS,eqn)
{


    real sourcePyr = pyrNucSource(c,t); 

    dS[eqn] = 0.0;
    C_UDMI(c,t,0) = sourcePyr;
    return sourcePyr;
}

DEFINE_SOURCE (m_0_NucSourceBenz,c,t,dS,eqn)
{



    real sourceBenz = benzNucSource(c,t); 

    dS[eqn] = 0.0;
    C_UDMI(c,t,0) = sourceBenz;
    return sourceBenz;
}

DEFINE_SOURCE (m_1_NucSourcePyr,c,t,dS,eqn)
{

    real sourcePyr = (2*NPyr)*pyrNucSource(c,t); 
    real source = sourcePyr; 


    dS[eqn] = 0.0;
    C_UDMI(c,t,1) = source;
    return source;
}

DEFINE_SOURCE (m_1_NucSourceBenz,c,t,dS,eqn)
{

    real sourceBenz = (2*NBenzene)*benzNucSource(c,t); 
    real source = sourceBenz;

    dS[eqn] = 0.0;
    C_UDMI(c,t,1) = source;
    return source;
}

DEFINE_SOURCE (m_2_NucSourcePyr,c,t,dS,eqn)
{

    real sourcePyr = (2*NPyr)*(2*NPyr)*pyrNucSource(c,t);
    real source = sourcePyr; 

    dS[eqn] = 0.0;
    C_UDMI(c,t,2) = source;
    return source;
}

DEFINE_SOURCE (m_2_NucSourceBenz,c,t,dS,eqn)
{

    real sourceBenz = (2*NBenzene)*(2*NBenzene)*benzNucSource(c,t); 
    real source = sourceBenz;

    dS[eqn] = 0.0;
    C_UDMI(c,t,2) = source;
    return source;
}

DEFINE_SOURCE (m_1_C2H2Source,c,t,dS,eqn)
{
    real source;
    real cellRho = C_R(c,t);
    real cellTemp = C_T(c,t);
    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2);
    real a = a1 + a2*cellTemp;
    real b = b1 + b2*cellTemp;
    /*real alpha = tanh(a/log10(cellM1/cellM0)+b); 
    if (alpha < 0 ) { alpha = 0.01; } */

    real alpha = 1.0;

    C_UDMI(c,t,19) = alpha;

    Material *mat1 = THREAD_MATERIAL(t);

    int idC2H2 = mixture_specie_index(mat1, "c2h2");
    real C2H2Mf = Pdf_Yi(c,t,idC2H2);
    real C2H2Mc = C2H2Mf*cellRho/C2H2MW;

    
    real k3 = A3*pow(cellTemp, n3)*exp(-Ea3/(R*cellTemp));
    real chiSoot = chiSootCalc(c,t,cellTemp,cellRho); 
    real Cs = pow((6*mC/(pi*rhoSoot)),oneThird);

    real muTwoThird = pow(10, lagInterp3mom(twoThird, 0.0, log10(cellM1/cellM0), log10(cellM2/cellM0))); 

    source = 2*k3*C2H2Mc*alpha*chiSoot*pi*pow(Cs,2)*muTwoThird*cellM0*fudgeFac2;
    dS[eqn] = 2*k3*C2H2Mc*alpha*chiSoot*pi*pow(Cs,2)*pow(cellM0,twoThird/3)*pow(cellM2,-oneThird/3)*pow(cellM1,-oneThird/3)*8.0/9.0*fudgeFac2;

    C_UDMI(c,t,3) = source;

    return source;
}

DEFINE_SOURCE (m_2_C2H2Source,c,t,dS,eqn)
{
    real source;
    real cellRho = C_R(c,t);
    real cellTemp = C_T(c,t);
    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2);
    real a = a1 + a2*cellTemp;
    real b = b1 + b2*cellTemp;
    /*real alpha = tanh(a/log10(cellM1/cellM0)+b); 
    if (alpha < 0 ) { alpha = 0.01; } */

    real alpha = 1.0;

    Material *mat = THREAD_MATERIAL(t);

    int idC2H2 = mixture_specie_index(mat, "c2h2");
    real C2H2Mf = Pdf_Yi(c,t,idC2H2);
    real C2H2Mc = C2H2Mf*cellRho/C2H2MW;
    
    real k3 = A3*pow(cellTemp, n3)*exp(-Ea3/(R*cellTemp));
    real chiSoot = chiSootCalc(c,t,cellTemp,cellRho);
    real Cs = pow((6*mC/(pi*rhoSoot)),oneThird);

    real MtwoThird = cellM0*pow(10, lagInterp3mom(twoThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real MfiveThird = cellM0*pow(10, lagInterp3mom(5*oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));

    source = k3*C2H2Mc*alpha*chiSoot*pi*pow(Cs,2)*(4*MtwoThird + 4*MfiveThird)*fudgeFac2;

    dS[eqn] = k3*C2H2Mc*alpha*chiSoot*pi*pow(Cs,2)*(-4*pow(cellM0,twoThird/3)*pow(cellM2,-11*oneThird/3)*\
        pow(cellM1,-8*oneThird/3)/9.0+4*5.0/9.0*pow(cellM0,-oneThird/3)*pow(cellM2,-4*oneThird/3)*pow(cellM1,5*oneThird/3))*fudgeFac2;

    C_UDMI(c,t,4) = source;

    return source;
}

DEFINE_SOURCE(m_1_OxSource,c,t,dS,eqn)
{
    real source;
    real sourceO2;
    real sourceOH;
    real cellRho = C_R(c,t);
    real cellTemp = C_T(c,t);
    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2);
    real a = a1 + a2*cellTemp;
    real b = b1 + b2*cellTemp;
    /*real alpha = tanh(a/log10(cellM1/cellM0)+b); 
    if (alpha < 0 ) { alpha = 0.01; }*/

    real alpha = 1.0;

    Material *mat = THREAD_MATERIAL(t);

    int idO2 = mixture_specie_index(mat, "o2");
    int idOH = mixture_specie_index(mat, "oh");

    real O2Mf = Pdf_Yi(c,t,idO2);
    real OHMf = Pdf_Yi(c,t,idOH);

    real O2Mc = O2Mf*cellRho/O2MW;
    real OHMc = OHMf*cellRho/OHMW;

    real k4 = A4*exp(-Ea4/(R*cellTemp));
    real chiSoot = chiSootCalc(c,t,cellTemp,cellRho);
    real Cs = pow((6*mC/(pi*rhoSoot)),oneThird);

    real MtwoThird = cellM0*pow(10, lagInterp3mom(twoThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));

    sourceO2 = -2*k4*O2Mc*alpha*chiSoot*pi*pow(Cs,2)*MtwoThird;
    sourceOH = -2*gammaOH*OHMc*sqrt(pi*kB*cellTemp/(2*mOH))*pi*pow(Cs,2)*MtwoThird;
    source = sourceO2 + sourceOH;

    C_UDMI(c,t,5) = sourceO2;
    C_UDMI(c,t,6) = sourceOH;

    C_UDMI(c,t,7) = source*fudgeFac1;

    dS[eqn] = -2*k4*O2Mc*alpha*chiSoot*pi*pow(Cs,2)*pow(cellM0,twoThird/3)*pow(cellM2,-oneThird/3)*pow(cellM1,-oneThird/3)*8.0/9.0;
    return source*fudgeFac1;
}

DEFINE_SOURCE(m_2_OxSource,c,t,dS,eqn)
{
    real source;
    real sourceO2;
    real sourceOH;
    real cellRho = C_R(c,t);
    real cellTemp = C_T(c,t);
    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2);
    real a = a1 + a2*cellTemp;
    real b = b1 + b2*cellTemp;
    /*real alpha = tanh(a/log10(cellM1/cellM0)+b); 
    if (alpha < 0 ) { alpha = 0.01; } */

    real alpha = 1.0;

    Material *mat = THREAD_MATERIAL(t);

    int idO2 = mixture_specie_index(mat, "o2");
    int idOH = mixture_specie_index(mat, "oh");

    real O2Mf = Pdf_Yi(c,t,idO2);
    real OHMf = Pdf_Yi(c,t,idOH);

    real O2Mc = O2Mf*cellRho/O2MW;
    real OHMc = OHMf*cellRho/OHMW;

    real k4 = A4*exp(-Ea4/(R*cellTemp));
    real chiSoot = chiSootCalc(c,t,cellTemp,cellRho);
    real Cs = pow((6*mC/(pi*rhoSoot)),oneThird);

    real MtwoThird = cellM0*pow(10, lagInterp3mom(twoThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real MfiveThird = cellM0*pow(10, lagInterp3mom(5*oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));

    sourceO2 = k4*O2Mc*alpha*chiSoot*pi*pow(Cs,2)*(4*MtwoThird - 4*MfiveThird);
    sourceOH = gammaOH*OHMc*sqrt(pi*kB*cellTemp/(2*mOH))*pi*pow(Cs,2)*(4*MtwoThird - 4*MfiveThird);
    source = sourceO2 + sourceOH;

    C_UDMI(c,t,8) = sourceO2;
    C_UDMI(c,t,9) = sourceOH;

    dS[eqn] = k4*O2Mc*alpha*chiSoot*pi*pow(Cs,2)*(-4*pow(cellM0,twoThird/3)*pow(cellM2,-11*oneThird/3)*\
        pow(cellM1,-8*oneThird/3)/9.0+4*5.0/9.0*pow(cellM0,-oneThird/3)*pow(cellM2,-4*oneThird/3)*pow(cellM1,5*oneThird/3));


    C_UDMI(c,t,10) = source*fudgeFac1;
    return source*fudgeFac1;
}

DEFINE_SOURCE (m_1_PAHsource,c,t,dS,eqn)
{
    real source;
    real cellRho = C_R(c,t);
    real cellTemp = C_T(c,t);
    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2);

    Material *mat = THREAD_MATERIAL(t);

    int ida2 = mixture_specie_index(mat,"a2");
    int ida3 = mixture_specie_index(mat,"a3");
    int ida4 = mixture_specie_index(mat,"a4");

    real a2Mf = Pdf_Yi(c,t,ida2);
    real a3Mf = Pdf_Yi(c,t,ida3);
    real a4Mf = Pdf_Yi(c,t,ida4);

    /*real PAHM0 = cellRho*avogad*(a2Mf/naphMW + a3Mf/phenanMW + a4Mf/pyrMW)/normParameter;
    real PAHM1 = cellRho*avogad*(a2Mf/naphMW*NNaph + a3Mf/phenanMW*NPhenan + a4Mf/pyrMW*NPyr)/normParameter;
    real PAHM2 = cellRho*avogad*(a2Mf/naphMW*pow(NNaph,2) + a3Mf/phenanMW*pow(NPhenan,2) + a4Mf/pyrMW*pow(NPyr,2))/normParameter;*/
    real PAHM0 = cellRho*avogad*(a4Mf/pyrMW)/normParameter;
    real PAHM1 = PAHM0*NPyr;
    real PAHM2 = PAHM1*NPyr; 

    C_UDMI(c,t,23) = PAHM0;
    C_UDMI(c,t,24) = PAHM1;
    C_UDMI(c,t,25) = PAHM2;

    real MoneThird = cellM0*pow(10, lagInterp3mom(oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real MtwoThird = cellM0*pow(10, lagInterp3mom(twoThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));

    C_UDMI(c,t,28) = MoneThird;
    C_UDMI(c,t,29) = MtwoThird;

    real MPAHoneHalf = 0.5*(PAHM0 + PAHM1);
    real MPAHthreeHalf = 0.5*(PAHM1 + PAHM2);

    /*real MPAHoneHalf = pow(PAHM0,-0.625)*pow(PAHM1,0.75)*pow(PAHM2,-0.125);
    real MPAHthreeHalf = PAHM0*pow(10, lagInterp3mom(3*oneHalf, log10(PAHM0/PAHM0), log10(PAHM1/PAHM0), log10(PAHM2/PAHM0)));
    C_UDMI(c,t,26) = MPAHoneHalf;
    real MPAHoneHalf = PAHM0*pow(10, lagInterp3mom(oneHalf, log10(PAHM0/PAHM0), log10(PAHM1/PAHM0), log10(PAHM2/PAHM0))); */

    C_UDMI(c,t,26) = MPAHoneHalf;
    C_UDMI(c,t,27) = MPAHthreeHalf;

    real Ch = dA*sqrt(twoThird);
    real Cs = pow((6*mC/(pi*rhoSoot)),oneThird);

    source = 0.5*vDW*sqrt(pi*kB*cellTemp/(2*mC))*(pow(Ch,2)*MPAHthreeHalf*cellM0 + 2*Ch*Cs*PAHM1*MoneThird + pow(Cs,2)*MPAHoneHalf*MtwoThird);

    dS[eqn] = 0.5*vDW*sqrt(pi*kB*cellTemp/(2*mC))*(2*Ch*Cs*PAHM1*5.0/9.0*pow(cellM0,5.0/9.0)*pow(cellM1,-4.0/9.0)*pow(cellM2,-1.0/9.0) +\
        pow(Cs,2)*MPAHoneHalf*2.0/9.0*pow(cellM0,8.0/9.0)*pow(cellM1,-7.0/9.0)*pow(cellM2,-1.0/9.0));

    C_UDMI(c,t,11) = source;
    return source;
}

DEFINE_SOURCE (m_2_PAHsource,c,t,dS,eqn)
{
    real source;
    real cellRho = C_R(c,t);
    real cellTemp = C_T(c,t);
    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2);

    Material *mat = THREAD_MATERIAL(t);

    int ida2 = mixture_specie_index(mat,"a2");
    int ida3 = mixture_specie_index(mat,"a3");
    int ida4 = mixture_specie_index(mat,"a4");

    real a2Mf = Pdf_Yi(c,t,ida2);
    real a3Mf = Pdf_Yi(c,t,ida3);
    real a4Mf = Pdf_Yi(c,t,ida4);

    /*real PAHM0 = cellRho*avogad*(a2Mf/naphMW + a3Mf/phenanMW + a4Mf/pyrMW)/normParameter;
    real PAHM1 = cellRho*avogad*(a2Mf/naphMW*NNaph + a3Mf/phenanMW*NPhenan + a4Mf/pyrMW*NPyr)/normParameter;
    real PAHM2 = cellRho*avogad*(a2Mf/naphMW*pow(NNaph,2) + a3Mf/phenanMW*pow(NPhenan,2) + a4Mf/pyrMW*pow(NPyr,2))/normParameter; */

    real PAHM0 = cellRho*avogad*(a4Mf/pyrMW)/normParameter;
    real PAHM1 = PAHM0*NPyr;
    real PAHM2 = PAHM1*NPyr;

    real MoneThird = cellM0*pow(10, lagInterp3mom(oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real MtwoThird = cellM0*pow(10, lagInterp3mom(twoThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real MfourThird = cellM0*pow(10, lagInterp3mom(2*twoThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real MfiveThird = cellM0*pow(10, lagInterp3mom(5*oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));

    /*real MPAHthreeHalf = PAHM0*pow(10, lagInterp3mom(3*oneHalf, log10(PAHM0/PAHM0), log10(PAHM1/PAHM0), log10(PAHM2/PAHM0)));
    real MPAHfiveHalf = PAHM0*pow(10, lagInterp3mom(5*oneHalf, log10(PAHM0/PAHM0), log10(PAHM1/PAHM0), log10(PAHM2/PAHM0)));
    real MPAHoneHalf = PAHM0*pow(10, lagInterp3mom(oneHalf, log10(PAHM0/PAHM0), log10(PAHM1/PAHM0), log10(PAHM2/PAHM0))); */

    real MPAHoneHalf = 0.5*(PAHM0 + PAHM1);
    real MPAHthreeHalf = 0.5*(PAHM1 + PAHM2);
    real MPAHfiveHalf = 1.5*PAHM2 - 0.5*PAHM1;

    real Ch = dA*sqrt(twoThird);
    real Cs = pow((6*mC/(pi*rhoSoot)),oneThird);

    source = 0.5*vDW*sqrt(pi*kB*cellTemp/(2*mC))*(pow(Ch,2)*MPAHfiveHalf*cellM0 + 2*Ch*Cs*PAHM2*MoneThird + pow(Cs,2)*MPAHthreeHalf*MtwoThird +\
        2*(pow(Ch,2)*MPAHthreeHalf*cellM1 + 2*Ch*Cs*PAHM1*MfourThird + pow(Cs,2)*MPAHoneHalf*MfiveThird));

    dS[eqn] = 0.5*vDW*sqrt(pi*kB*cellTemp/(2*mC))*(-2*Ch*Cs*PAHM2*1.0/9.0*pow(cellM0,5.0/9.0)*pow(cellM1,5.0/9.0)*pow(cellM2,-10.0/9.0) - \
        pow(Cs,2)*MPAHthreeHalf*1.0/9.0*pow(cellM0,8.0/9.0)*pow(cellM1,2.0/9.0)*pow(cellM2,-10.0/9.0) + \
        2*(2*Ch*Cs*PAHM1*2.0/9.0*pow(cellM0,-1.0/9.0)*pow(cellM1,8.0/9.0)*pow(cellM2,-4.0/9.0)) + \
        pow(Cs,2)*MPAHoneHalf*5.0/9.0*pow(cellM0,-1.0/9.0)*pow(cellM1,5.0/9.0)*pow(cellM2,-4.0/9.0));

    C_UDMI(c,t,12) = source;
    return source;
}


DEFINE_SOURCE (m_0_CoagSource,c,t,dS,eqn)
{
    real source;
    real cellPressure = C_P(c,t) + Patm;
    real cellRho = C_R(c,t);
    real cellTemp = C_T(c,t);
    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2);
    real lamVisc = C_MU_L(c,t);

    real muOneThird = pow(10, lagInterp3mom(oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));

    real meanFreePath = kB*cellTemp/(sqrt(2)*pi*pow(dAir,2)*cellPressure);
    real avgDia = pow(10, lagInterp3mom(oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)))*dC; 
    real Kn = 2*meanFreePath/avgDia;

    C_UDMI(c,t,20) = meanFreePath; 
    C_UDMI(c,t,21) = avgDia;
    C_UDMI(c,t,22) = Kn;

    real Kc = 2*kB*cellTemp/(3*lamVisc);
    real KcPrime = 2.514*meanFreePath*pow(pi*rhoSoot/6,oneThird);
    real Kf = vDW*sqrt(6*kB*cellTemp/rhoSoot)*pow(3*mC/(4*pi*rhoSoot),oneSixth);
    
    real muTwoThird = pow(10, lagInterp3mom(twoThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muFourThird = pow(10, lagInterp3mom(4*oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muFiveThird = pow(10, lagInterp3mom(5*oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));

    real muSevenThird = pow(10, lagInterp3mom(7*oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muEightThird = pow(10, lagInterp3mom(8*oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muMinusOneThird = pow(10, lagInterp3mom(-oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muMinusTwoThird = pow(10, lagInterp3mom(-twoThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));

    real muMinusOneHalf = pow(10, lagInterp3mom(-oneHalf, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muMinusOneSixth = pow(10, lagInterp3mom(-oneSixth, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muOneSixth = pow(10, lagInterp3mom(oneSixth, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muOneHalf = pow(10, lagInterp3mom(oneHalf, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muThreeHalf = pow(10, lagInterp3mom(3*oneHalf, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muFiveSixth = pow(10, lagInterp3mom(5*oneSixth, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muSevenSixth = pow(10, lagInterp3mom(7*oneSixth, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muElevenSixth = pow(10, lagInterp3mom(11*oneSixth, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muThirteenSixth = pow(10, lagInterp3mom(13*oneSixth, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));



    real f_0_0_0 = 2*(muMinusOneHalf*muOneSixth + pow(muMinusOneSixth,2));
    real f_1_0_0 = 2*(muMinusOneHalf*muSevenSixth + 2*muMinusOneSixth*muFiveSixth + muOneSixth*muOneHalf);
    real f_2_0_0 = 2*(muMinusOneHalf*muThirteenSixth + 2*muMinusOneSixth*muElevenSixth + muOneSixth*muThreeHalf + 2*muOneHalf*muSevenSixth + 2*pow(muFiveSixth,2));

    real f_oneHalf_0_0 = pow(10, lagInterp3mom(oneHalf, log10(f_0_0_0), log10(f_1_0_0), log10(f_2_0_0)));

    real Gc = -Kc*(1 + muOneThird*muMinusOneThird + KcPrime*(muMinusOneThird + muOneThird*muMinusTwoThird))*pow(cellM0,2)*normParameter;
    real Gf = -0.5*Kf*pow(cellM0,2)*f_oneHalf_0_0*normParameter;

    
    C_UDMI(c,t,13) = Gc;
    C_UDMI(c,t,14) = Gf;


    if (Kn < 0.1) { source = Gc; }

    else if (Kn > 8.0) { source = Gf; }

    else { source = Gc*Gf/(Gc + Gf); }

    dS[eqn] = 0.0;

    C_UDMI(c,t,15) = source;
    return source;

}

DEFINE_SOURCE (m_2_CoagSource,c,t,dS,eqn)
{

    real source;
    real cellPressure = C_P(c,t) + Patm;
    real cellRho = C_R(c,t);
    real cellTemp = C_T(c,t);
    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2);
    real lamVisc = C_MU_L(c,t);

    real muOneThird = pow(10, lagInterp3mom(oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));

    real meanFreePath = kB*cellTemp/(sqrt(2)*pi*pow(dAir,2)*cellPressure);
    real avgDia = pow(10, lagInterp3mom(oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)))*dC;
    real Kn = 2*meanFreePath/avgDia;
    real Kc = 2*kB*cellTemp/(3*lamVisc);
    real KcPrime = 2.514*meanFreePath*pow(pi*rhoSoot/6,oneThird);
    real Kf = vDW*sqrt(6*kB*cellTemp/rhoSoot)*pow(3*mC/(4*pi*rhoSoot),oneSixth);
    
    real muTwoThird = pow(10, lagInterp3mom(twoThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muFourThird = pow(10, lagInterp3mom(4*oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muFiveThird = pow(10, lagInterp3mom(5*oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muSevenThird = pow(10, lagInterp3mom(7*oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muEightThird = pow(10, lagInterp3mom(8*oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muMinusOneHalf = pow(10, lagInterp3mom(-oneHalf, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muMinusOneSixth = pow(10, lagInterp3mom(-oneSixth, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muOneSixth = pow(10, lagInterp3mom(oneSixth, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muOneHalf = pow(10, lagInterp3mom(oneHalf, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muThreeHalf = pow(10, lagInterp3mom(3*oneHalf, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muFiveSixth = pow(10, lagInterp3mom(5*oneSixth, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muSevenSixth = pow(10, lagInterp3mom(7*oneSixth, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muElevenSixth = pow(10, lagInterp3mom(11*oneSixth, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muThirteenSixth = pow(10, lagInterp3mom(13*oneSixth, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muNineteenSixth = pow(10, lagInterp3mom(19*oneSixth, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muSeventeenSixth = pow(10, lagInterp3mom(17*oneSixth, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));
    real muFiveHalf = pow(10, lagInterp3mom(5*oneHalf, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));

    real f_0_1_1 = 2*muSevenSixth*muOneHalf + 4*muFiveSixth*muFiveSixth;
    real f_1_1_1 = 2*muThirteenSixth*muOneHalf + 4*muElevenSixth*muFiveSixth + 2*muThreeHalf*muSevenSixth;
    real f_2_1_1 = 2*muNineteenSixth*muOneHalf + 4*muSeventeenSixth*muFiveSixth + 2*muFiveHalf*muSevenSixth + 4*muThirteenSixth*muThreeHalf + 4*muElevenSixth*muElevenSixth;
    real f_oneHalf_1_1 = pow(10, lagInterp3mom(oneHalf, log10(f_0_1_1), log10(f_1_1_1), log10(f_2_1_1)));

    real Gc = Kc*(2*cellM1*cellM1/(cellM0*cellM0) + 2*muFourThird*muTwoThird + KcPrime*(2*muTwoThird*muOneThird + 2*muFourThird*muOneThird))*pow(cellM0,2)*normParameter;

    real Gf = Kf*pow(cellM0,2)*f_oneHalf_1_1*normParameter;

    C_UDMI(c,t,16) = Gc;
    C_UDMI(c,t,17) = Gf;

    if (Kn < 0.1) { source = Gc; }

    else if (Kn > 8.0) { source = Gf; }

    else { source = Gc*Gf/(Gc + Gf); }

    dS[eqn] = 0.0;
    
    C_UDMI(c,t,18) = source;

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
    int idBenzene = mixture_specie_index(mat, "a1");
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
    real k2 = A2; 
    real k3 = A3*pow(cellTemp, n3)*exp(-Ea3/(R*cellTemp));
    real k4 = A4*exp(-Ea4/(R*cellTemp));

    real numerator = k1f*HMc;
    real denominator = k1r*H2Mc + k2*HMc + k3*C2H2Mc + k4*O2Mc;

    real chiSoot = numerator/(denominator + 1.0) * siteDensity;

    return chiSoot;
}

DEFINE_DIFFUSIVITY(uds_diff,c,t,i)
{
    real diff = C_MU_EFF(c,t)/(ScT*PrT); 

    return diff;
}

DEFINE_PROPERTY(abs_coeff_Widmann,c,t)
{
    real cellTemp = C_T(c,t);
    real cellRho = C_R(c,t);
    real cellM1 = C_UDSI(c,t,1);
    real sootfv = cellM1*mC*cellRho/rhoSoot*normParameter;

    real absCoeff = 2370*cellTemp*sootfv;
    return absCoeff;
}

DEFINE_PROPERTY(abs_coeff_Sazhin,c,t)
{
    real cellTemp = C_T(c,t);
    real cellRho = C_R(c,t);
    real cellM1 = C_UDSI(c,t,1);
    real sootfv = cellM1*mC*cellRho/rhoSoot*normParameter;

    real absCoeff = 1232.4*rhoSoot*sootfv*(1+4.8e-4*(cellTemp-2000));
    return absCoeff;
}

DEFINE_PROPERTY(thermal_conductivity,c,t)
{
    real cellTemp = C_T(c,t);
    real Tr = cellTemp/132.52;


    real thermCond = 4.358e-3*(33.972*pow(Tr,-1) - 164.702*pow(Tr,-twoThird) + 262.108*pow(Tr,-oneThird) 
        - 21.534 - 443.455*pow(Tr,oneThird) + 607.339*pow(Tr, twoThird) - 368.79*Tr + 111.296*pow(Tr, 2*twoThird) 
        - 13.412*pow(Tr, 5*oneThird));

    return thermCond;
}