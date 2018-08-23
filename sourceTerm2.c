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
#define gammaAcet 1e-10
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
#define A1f 4.2e10
#define Ea1f 13000*4.18
#define A1r 3.9e9
#define Ea1r 11000*4.18
#define A2f 1.0e7
#define n2f 0.734
#define Ea2f 1430*4.18
#define A2r 3.68e5
#define n2r 1.139
#define Ea2r 17100*4.18
#define A3 2.0e10
#define A4 8.0e4
#define n4 1.56
#define Ea4 3800*4.18
#define A5 2.2e9
#define Ea5 7500*4.18
#define gammaOH 0.13

/* Fudge Factors */
#define fudgeFac1 1.0
#define fudgeFac2 1.0

real lagInterp3mom(real p, real m0, real m1, real m2); /* Lagragian interpolation */
real logInterp(real p, real a, real b, real M, real N); /* logarithmic interpolation */
real pyrNucSource(cell_t c,Thread *t);
real benzNucSource(cell_t c,Thread *t);
real acetNucSource(cell_t c,Thread *t);
real chiSootCalc(cell_t c,Thread *t, real cellTemp, real cellRho);
real pdfpyrNucSource(cell_t c,Thread *t);
real pdfChiSootCalc(real Temp, real Pstat, real RGas, real HMf, real OHMf, real H2Mf, real H2OMf, real C2H2Mf, real O2Mf);
real HACAalphaCalc(real Temp, real cellM0, real cellM1);
real fractionalMoments(real p, real m0, real m1, real m2);
real HACAIntegrandCalc(real Temp, real TempAvg, real TempRMS, real C2H2Mf, real C2H2MfAvg, real C2H2MfRMS, real alpha, real ChiSoot);
real O2IntegrandCalc(real Temp, real TempAvg, real TempRMS, real O2Mf, real O2MfAvg, real O2MfRMS, real alpha, real ChiSoot);
real gaussianMonoPDFCalc(real Var, real VarMean, real VarStd);
real fracMomFirstDiffM0(real p, real m0, real m1, real m2);
real fracMomFirstDiffM1(real p, real m0, real m1, real m2);
real fracMomFirstDiffM2(real p, real m0, real m1, real m2);
real factorialCalc(n);

DEFINE_SOURCE (m_0_NucSourceAcet,c,t,dS,eqn)
{


    real sourceAcet = acetNucSource(c,t); 

    dS[eqn] = 0.0;
    C_UDMI(c,t,0) = sourceAcet;
    return sourceAcet;
}


DEFINE_SOURCE (m_0_NucSourcePyr,c,t,dS,eqn)
{


    real sourcePyr = pyrNucSource(c,t); 

    dS[eqn] = 0.0;
    C_UDMI(c,t,0) = sourcePyr;
    return sourcePyr;
}

DEFINE_SOURCE (pdf_m_0_NucSourcePyr,c,t,dS,eqn)
{


    real sourcePyr = pdfpyrNucSource(c,t); 

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

DEFINE_SOURCE (pdf_m_1_NucSourcePyr,c,t,dS,eqn)
{

    real sourcePyr = (2*NPyr)*pdfpyrNucSource(c,t); 
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

DEFINE_SOURCE (pdf_m_2_NucSourcePyr,c,t,dS,eqn)
{

    real sourcePyr = (2*NPyr)*(2*NPyr)*pdfpyrNucSource(c,t);
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
    real alpha = tanh(a/log10(cellM1/cellM0)+b); 
    if (alpha < 0 ) { alpha = 0.01; }

    C_UDMI(c,t,19) = alpha;

    Material *mat1 = THREAD_MATERIAL(t);

    int idC2H2 = mixture_specie_index(mat1, "c2h2");
    real C2H2Mf = Pdf_Yi(c,t,idC2H2);
    real C2H2Mc = C2H2Mf*cellRho/C2H2MW;

    
    real k4 = A4*pow(cellTemp, n4)*exp(-Ea4/(R*cellTemp));
    real chiSoot = chiSootCalc(c,t,cellTemp,cellRho); 
    real Cs = pow((6*mC/(pi*rhoSoot)),oneThird);

    real muTwoThird = pow(10, lagInterp3mom(twoThird, 0.0, log10(cellM1/cellM0), log10(cellM2/cellM0))); 

    source = 2*k4*C2H2Mc*alpha*chiSoot*pi*pow(Cs,2)*muTwoThird*cellM0*fudgeFac2;
    dS[eqn] = 2*k4*C2H2Mc*alpha*chiSoot*pi*pow(Cs,2)*pow(cellM0,twoThird/3)*pow(cellM2,-oneThird/3)*pow(cellM1,-oneThird/3)*8.0/9.0*fudgeFac2;

    C_UDMI(c,t,3) = source;

    return source;
}

DEFINE_SOURCE (pdf_m_1_C2H2Source,c,t,dS,eqn)
{
    real source = 0.0;
    real integrand = 0.0;
    real pdfIntegrand = 0.0;
    real temp;
    real C2H2;
    real C2H2_new;
    real RGas = C_RGAS(c,t) * 1e-3;
    real cellRho = C_R(c,t);

    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2);

    Material *mat = THREAD_MATERIAL(t);
    int idC2H2 = mixture_specie_index(mat, "c2h2");

    real cellPressure = (C_P(c,t) + Patm)/1e3;
    real cellTempAvg = C_T(c,t);
    real cellTempRMS = C_UDMI(c,t,36);
    real C2H2MfAvg = Pdf_Yi(c,t,idC2H2);
    real C2H2MfRMS = C_UDMI(c,t,38);

    int idH2 = mixture_specie_index(mat, "h2");
    int idH = mixture_specie_index(mat, "h");
    int idO2 = mixture_specie_index(mat, "o2");
    int idOH = mixture_specie_index(mat, "oh");
    int idH2O = mixture_specie_index(mat, "h2o");

    real H2Mf = Pdf_Yi(c,t,idH2);
    real HMf = Pdf_Yi(c,t,idH);
    real O2Mf = Pdf_Yi(c,t,idO2);
    real OHMf = Pdf_Yi(c,t,idOH);
    real H2OMf = Pdf_Yi(c,t,idH2O);

    real Cs = pow((6*mC/(pi*rhoSoot)),oneThird);
    real muTwoThird = fractionalMoments(twoThird, cellM0, cellM1, cellM2); 

    real multConstant = 2 * A4 * cellPressure/(RGas*C2H2MW)*pow(Cs, 2) * pi * muTwoThird * cellM0;
    real tempRatio = cellTempRMS/(cellTempAvg);
    real C2H2MfRatio = C2H2MfRMS/(C2H2MfAvg + 1e-7);

    C_UDMI(c,t,33) = tempRatio;
    C_UDMI(c,t,34) = C2H2MfRatio;
    

    if (tempRatio < 0.02 && C2H2MfRatio < 0.02)
    {
        real alpha = HACAalphaCalc(cellTempAvg, cellM0, cellM1);
        real chiSoot = pdfChiSootCalc(cellTempAvg, cellPressure, RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);
        integrand = pow(cellTempAvg, n4) / cellTempAvg * exp(-Ea4/(R*cellTempAvg)) * alpha * chiSoot * C2H2MfAvg;
        pdfIntegrand = 1.0;
        C_UDMI(c,t,32) = 5.0;

    }

    else if (tempRatio >= 0.02 && C2H2MfRatio < 0.02)
    {

        real tempLower;

        real tempUpper = cellTempAvg + 3.5 * cellTempRMS;

        if (cellTempAvg - 3.5 * cellTempRMS > 300) { tempLower = cellTempAvg - 3.5 * cellTempRMS; }
        else { tempLower = 300; }

        real tempInt = 0.5 * cellTempRMS;

        real alpha_low = HACAalphaCalc(tempLower, cellM0, cellM1);
        real alpha_high = HACAalphaCalc(tempUpper, cellM0, cellM1);

        real ChiSoot_low = pdfChiSootCalc(tempLower, cellPressure, RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);
        real ChiSoot_high = pdfChiSootCalc(tempUpper, cellPressure, RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);

        real f_low = pow(tempLower, n4)/ tempLower * exp(-Ea4/(R*tempLower)) * alpha_low * ChiSoot_low * \
        gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS);

        real P_low = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS);

        real f_high = pow(tempUpper, n4) / tempLower * exp(-Ea4/(R*tempUpper)) * alpha_high * ChiSoot_high * \
        gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS);

        real P_high = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS);

        real f_xi = 0.0; real P_xi = 0.0;
        

        int ctr;

        for (ctr = 1; ctr < 13; ctr ++)
        {
            temp = tempLower + ctr * tempInt;
            real alpha = HACAalphaCalc(temp, cellM0, cellM1);
            real ChiSoot = pdfChiSootCalc(temp, cellPressure,RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);

            f_xi = f_xi + pow(temp, n4 - 1.0) * exp(-Ea4/(R*temp)) * alpha * ChiSoot * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
            P_xi = P_xi + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
        }

        pdfIntegrand = (P_low + P_high + 2 * P_xi) * tempInt / 2.0;
        integrand = C2H2MfAvg * (f_low + f_high + 2 * f_xi) / pdfIntegrand;
        
        C_UDMI(c,t,32) = 4.0;
    }

    else if (tempRatio < 0.02 && C2H2MfRatio >= 0.02)
    {
        real alpha = HACAalphaCalc(cellTempAvg, cellM0, cellM1);
        real chiSoot = pdfChiSootCalc(cellTempAvg, cellPressure, RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);
        real C2H2Lower; real C2H2Upper;

        if (C2H2MfAvg - 3.5 * C2H2MfRMS > 0) { C2H2Lower = C2H2MfAvg - 3.5 * C2H2MfRMS; }
        else { C2H2Lower = 0.0; }
        if (C2H2MfAvg + 3.5 * C2H2MfRMS < 1.0) { C2H2Upper = C2H2MfAvg + 3.5 * C2H2MfRMS; }
        else { C2H2Upper = 1.0; }

        real C2H2Int = 0.5 * C2H2MfRMS;

        int ctr;

        real f_low = C2H2Lower * gaussianMonoPDFCalc(C2H2Lower, C2H2MfAvg, C2H2MfRMS);
        real P_low = gaussianMonoPDFCalc(C2H2Lower, C2H2MfAvg, C2H2MfRMS);

        real f_high = C2H2Upper * gaussianMonoPDFCalc(C2H2Upper, C2H2MfAvg, C2H2MfRMS);
        real P_high = gaussianMonoPDFCalc(C2H2Upper, C2H2MfAvg, C2H2MfRMS);

        real f_yi = 0.0; real P_yi = 0.0;

        for (ctr = 1; ctr < 13; ctr ++)
        {
            C2H2 = C2H2Lower + C2H2Int * ctr;
            f_yi = f_yi + C2H2 * gaussianMonoPDFCalc(C2H2, C2H2MfAvg, C2H2MfRMS);
            P_yi = P_yi + gaussianMonoPDFCalc(C2H2, C2H2MfAvg, C2H2MfRMS);
        }

        pdfIntegrand = (P_high + P_low + 2 * P_yi) * C2H2Int/2.0;
        integrand = pow(cellTempAvg, n4 - 1.0) * exp(-Ea4/(R*cellTempAvg)) * alpha * chiSoot * (f_low + f_high + 2 * f_yi) / pdfIntegrand;
        C_UDMI(c,t,32) = 3.0;
    }


    else 
    {
        real tempUpper; real tempLower;
        real C2H2Upper; real C2H2Lower;

        tempUpper = cellTempAvg + 3.5 * cellTempRMS;
        if (cellTempAvg - 3.5 * cellTempRMS > 300) { tempLower = cellTempAvg - 3.5 * cellTempRMS; }
        else { tempLower = 300; }

        if (C2H2MfAvg - 3.5 * C2H2MfRMS > 0) { C2H2Lower = C2H2MfAvg - 3.5 * C2H2MfRMS; }
        else { C2H2Lower = 0.0; }

        if (C2H2MfAvg + 3.5 * C2H2MfRMS < 1.0) { C2H2Upper = C2H2MfAvg + 3.5 * C2H2MfRMS; }
        else { C2H2Upper = 1.0; }

        real tempInt = 0.5 * cellTempRMS; real C2H2Int = 0.5 * C2H2MfRMS;

        real alpha_low = HACAalphaCalc(tempLower, cellM0, cellM1);
        real alpha_high = HACAalphaCalc(tempUpper, cellM0, cellM1);

        real ChiSoot_low = pdfChiSootCalc(tempLower, cellPressure,RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);
        real ChiSoot_high = pdfChiSootCalc(tempUpper, cellPressure,RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);

        real f_low_low = HACAIntegrandCalc(tempLower, cellTempAvg, cellTempRMS, C2H2Lower, C2H2MfAvg, C2H2MfRMS, alpha_low, ChiSoot_low);
        real P_low_low = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(C2H2Lower, C2H2MfAvg, C2H2MfRMS);

        real f_high_low = HACAIntegrandCalc(tempUpper, cellTempAvg, cellTempRMS, C2H2Lower, C2H2MfAvg, C2H2MfRMS, alpha_high, ChiSoot_high);
        real P_high_low = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(C2H2Lower, C2H2MfAvg, C2H2MfRMS);

        real f_low_high = HACAIntegrandCalc(tempLower, cellTempAvg, cellTempRMS, C2H2Upper, C2H2MfAvg, C2H2MfRMS, alpha_low, ChiSoot_low);
        real P_low_high = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(C2H2Upper, C2H2MfAvg, C2H2MfRMS);

        real f_high_high = HACAIntegrandCalc(tempUpper, cellTempAvg, cellTempRMS, C2H2Upper, C2H2MfAvg, C2H2MfRMS, alpha_high, ChiSoot_high);
        real P_high_high = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(C2H2Upper, C2H2MfAvg, C2H2MfRMS);

        real f_xi_low = 0.0; real P_xi_low = 0.0;
        real f_xi_high = 0.0; real P_xi_high = 0.0;
        real f_low_yi = 0.0; real P_low_yi = 0.0;
        real f_high_yi = 0.0; real P_high_yi = 0.0;
        real f_xi_yi = 0.0; real P_xi_yi = 0.0;

        int ctr;
        int ctr1;

        for (ctr = 1; ctr < 13; ctr++)
        {
            temp = tempLower + ctr * tempInt; C2H2 = C2H2Lower + ctr * C2H2Int;
            real alpha = HACAalphaCalc(temp, cellM0, cellM1);
            real ChiSoot = pdfChiSootCalc(temp, cellPressure,RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);
            f_xi_low = f_xi_low + HACAIntegrandCalc(temp, cellTempAvg, cellTempRMS, C2H2Lower, C2H2MfAvg, C2H2MfRMS, alpha, ChiSoot);
            P_xi_low = P_xi_low + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(C2H2Lower, C2H2MfAvg, C2H2MfRMS);

            f_xi_high = f_xi_high + HACAIntegrandCalc(temp, cellTempAvg, cellTempRMS, C2H2Upper, C2H2MfAvg, C2H2MfRMS, alpha, ChiSoot);
            P_xi_high = P_xi_high + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(C2H2Upper, C2H2MfAvg, C2H2MfRMS);

            f_low_yi = f_low_yi + HACAIntegrandCalc(tempLower, cellTempAvg, cellTempRMS, C2H2, C2H2MfAvg, C2H2MfRMS, alpha_low, ChiSoot_low);
            P_low_yi = P_low_yi + gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(C2H2, C2H2MfAvg, C2H2MfRMS);

            f_high_yi = f_high_yi + HACAIntegrandCalc(tempUpper, cellTempAvg, cellTempRMS, C2H2, C2H2MfAvg, C2H2MfRMS, alpha_high, ChiSoot_high);
            P_high_yi = P_high_yi + gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(C2H2, C2H2MfAvg, C2H2MfRMS);

            for (ctr1 = 1; ctr1 < 13; ctr1++)
            {
                C2H2_new = C2H2Lower + ctr1 * C2H2Int;
                f_xi_yi = f_xi_yi + HACAIntegrandCalc(temp, cellTempAvg, cellTempRMS, C2H2_new, C2H2MfAvg, C2H2MfRMS, alpha, ChiSoot);
                P_xi_yi = P_xi_yi + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(C2H2_new, C2H2MfAvg, C2H2MfRMS);
            }
        }
        pdfIntegrand = 0.25 * tempInt * C2H2Int * (P_low_high + P_high_high + P_high_low + P_low_low + \
            2 * (P_low_yi + P_high_yi + P_xi_low + P_xi_high) + 4 * P_xi_yi);

        integrand = 0.25 * tempInt * C2H2Int * (f_low_high + f_high_high + f_high_low + f_low_low + \
            2 * (f_low_yi + f_high_yi + f_xi_low + f_xi_high) + 4 * f_xi_yi);

        integrand = integrand / pdfIntegrand;
        C_UDMI(c,t,32) = 2.0;

    }

    C_UDMI(c,t,33) = pdfIntegrand;
    

    source = multConstant * integrand;

    real muTwoThirdDiff = fracMomFirstDiffM1(twoThird, cellM0, cellM1, cellM2);

    dS[eqn] = source / muTwoThird * muTwoThirdDiff ;

    /* dS[eqn] = 2 * A4 * cellPressure/(RGas*C2H2MW)*pow(Cs, 2) * pi * muTwoThirdDiff * integrand; */

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

    source = k4*C2H2Mc*alpha*chiSoot*pi*pow(Cs,2)*(4*MtwoThird + 4*MfiveThird)*fudgeFac2;

    dS[eqn] = k4*C2H2Mc*alpha*chiSoot*pi*pow(Cs,2)*(-4*pow(cellM0,twoThird/3)*pow(cellM2,-11*oneThird/3)*\
        pow(cellM1,-8*oneThird/3)/9.0+4*5.0/9.0*pow(cellM0,-oneThird/3)*pow(cellM2,-4*oneThird/3)*pow(cellM1,5*oneThird/3))*fudgeFac2;

    C_UDMI(c,t,4) = source;

    return source;
}

DEFINE_SOURCE (pdf_m_2_C2H2Source,c,t,dS,eqn)
{
    real source = 0.0;
    real integrand = 0.0;
    real pdfIntegrand = 0.0;
    real temp;
    real C2H2;
    real C2H2_new;
    real RGas = C_RGAS(c,t) * 1e-3;
    real cellRho = C_R(c,t);

    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2);

    Material *mat = THREAD_MATERIAL(t);
    int idC2H2 = mixture_specie_index(mat, "c2h2");

    real cellPressure = (C_P(c,t) + Patm)/1e3;
    real cellTempAvg = C_T(c,t);
    real cellTempRMS = C_UDMI(c,t,36);
    real C2H2MfAvg = Pdf_Yi(c,t,idC2H2);
    real C2H2MfRMS = C_UDMI(c,t,38);

    int idH2 = mixture_specie_index(mat, "h2");
    int idH = mixture_specie_index(mat, "h");
    int idO2 = mixture_specie_index(mat, "o2");
    int idOH = mixture_specie_index(mat, "oh");
    int idH2O = mixture_specie_index(mat, "h2o");

    real H2Mf = Pdf_Yi(c,t,idH2);
    real HMf = Pdf_Yi(c,t,idH);
    real O2Mf = Pdf_Yi(c,t,idO2);
    real OHMf = Pdf_Yi(c,t,idOH);
    real H2OMf = Pdf_Yi(c,t,idH2O);

    real Cs = pow((6*mC/(pi*rhoSoot)),oneThird);
    real muTwoThird = fractionalMoments(twoThird, cellM0, cellM1, cellM2);
    real muFiveThird = fractionalMoments(5*oneThird, cellM0, cellM1, cellM2); 

    real multConstant = A4 * cellPressure/(RGas*C2H2MW)*pow(Cs, 2) * pi * (4 * muTwoThird + 4 * muFiveThird) * cellM0;
    real tempRatio = cellTempRMS/(cellTempAvg);
    real C2H2MfRatio = C2H2MfRMS/(C2H2MfAvg + 1e-7);

    C_UDMI(c,t,33) = tempRatio;
    C_UDMI(c,t,34) = C2H2MfRatio;
    

    if (tempRatio < 0.02 && C2H2MfRatio < 0.02)
    {
        real alpha = HACAalphaCalc(cellTempAvg, cellM0, cellM1);
        real chiSoot = pdfChiSootCalc(cellTempAvg, cellPressure, RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);
        integrand = pow(cellTempAvg, n4) / cellTempAvg * exp(-Ea4/(R*cellTempAvg)) * alpha * chiSoot * C2H2MfAvg;
        pdfIntegrand = 1.0;
        C_UDMI(c,t,32) = 5.0;

    }

    else if (tempRatio >= 0.02 && C2H2MfRatio < 0.02)
    {

        real tempLower;

        real tempUpper = cellTempAvg + 3.5 * cellTempRMS;

        if (cellTempAvg - 3.5 * cellTempRMS > 300) { tempLower = cellTempAvg - 3.5 * cellTempRMS; }
        else { tempLower = 300; }

        real tempInt = 0.5 * cellTempRMS;

        real alpha_low = HACAalphaCalc(tempLower, cellM0, cellM1);
        real alpha_high = HACAalphaCalc(tempUpper, cellM0, cellM1);

        real ChiSoot_low = pdfChiSootCalc(tempLower, cellPressure, RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);
        real ChiSoot_high = pdfChiSootCalc(tempUpper, cellPressure, RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);

        real f_low = pow(tempLower, n4)/ tempLower * exp(-Ea4/(R*tempLower)) * alpha_low * ChiSoot_low * \
        gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS);

        real P_low = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS);

        real f_high = pow(tempUpper, n4) / tempLower * exp(-Ea4/(R*tempUpper)) * alpha_high * ChiSoot_high * \
        gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS);

        real P_high = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS);

        real f_xi = 0.0; real P_xi = 0.0;
        

        int ctr;

        for (ctr = 1; ctr < 13; ctr ++)
        {
            temp = tempLower + ctr * tempInt;
            real alpha = HACAalphaCalc(temp, cellM0, cellM1);
            real ChiSoot = pdfChiSootCalc(temp, cellPressure,RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);

            f_xi = f_xi + pow(temp, n4 - 1.0) * exp(-Ea4/(R*temp)) * alpha * ChiSoot * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
            P_xi = P_xi + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
        }

        pdfIntegrand = (P_low + P_high + 2 * P_xi) * tempInt / 2.0;
        integrand = C2H2MfAvg * (f_low + f_high + 2 * f_xi) / pdfIntegrand;
        
        C_UDMI(c,t,32) = 4.0;
    }

    else if (tempRatio < 0.02 && C2H2MfRatio >= 0.02)
    {
        real alpha = HACAalphaCalc(cellTempAvg, cellM0, cellM1);
        real chiSoot = pdfChiSootCalc(cellTempAvg, cellPressure, RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);
        real C2H2Lower; real C2H2Upper;

        if (C2H2MfAvg - 3.5 * C2H2MfRMS > 0) { C2H2Lower = C2H2MfAvg - 3.5 * C2H2MfRMS; }
        else { C2H2Lower = 0.0; }
        if (C2H2MfAvg + 3.5 * C2H2MfRMS < 1.0) { C2H2Upper = C2H2MfAvg + 3.5 * C2H2MfRMS; }
        else { C2H2Upper = 1.0; }

        real C2H2Int = 0.5 * C2H2MfRMS;

        int ctr;

        real f_low = C2H2Lower * gaussianMonoPDFCalc(C2H2Lower, C2H2MfAvg, C2H2MfRMS);
        real P_low = gaussianMonoPDFCalc(C2H2Lower, C2H2MfAvg, C2H2MfRMS);

        real f_high = C2H2Upper * gaussianMonoPDFCalc(C2H2Upper, C2H2MfAvg, C2H2MfRMS);
        real P_high = gaussianMonoPDFCalc(C2H2Upper, C2H2MfAvg, C2H2MfRMS);

        real f_yi = 0.0; real P_yi = 0.0;

        for (ctr = 1; ctr < 13; ctr ++)
        {
            C2H2 = C2H2Lower + C2H2Int * ctr;
            f_yi = f_yi + C2H2 * gaussianMonoPDFCalc(C2H2, C2H2MfAvg, C2H2MfRMS);
            P_yi = P_yi + gaussianMonoPDFCalc(C2H2, C2H2MfAvg, C2H2MfRMS);
        }

        pdfIntegrand = (P_high + P_low + 2 * P_yi) * C2H2Int/2.0;
        integrand = pow(cellTempAvg, n4 - 1.0) * exp(-Ea4/(R*cellTempAvg)) * alpha * chiSoot * (f_low + f_high + 2 * f_yi) / pdfIntegrand;
        C_UDMI(c,t,32) = 3.0;
    }


    else 
    {
        real tempUpper; real tempLower;
        real C2H2Upper; real C2H2Lower;

        tempUpper = cellTempAvg + 3.5 * cellTempRMS;
        if (cellTempAvg - 3.5 * cellTempRMS > 300) { tempLower = cellTempAvg - 3.5 * cellTempRMS; }
        else { tempLower = 300; }

        if (C2H2MfAvg - 3.5 * C2H2MfRMS > 0) { C2H2Lower = C2H2MfAvg - 3.5 * C2H2MfRMS; }
        else { C2H2Lower = 0.0; }

        if (C2H2MfAvg + 3.5 * C2H2MfRMS < 1.0) { C2H2Upper = C2H2MfAvg + 3.5 * C2H2MfRMS; }
        else { C2H2Upper = 1.0; }

        real tempInt = 0.5 * cellTempRMS; real C2H2Int = 0.5 * C2H2MfRMS;

        real alpha_low = HACAalphaCalc(tempLower, cellM0, cellM1);
        real alpha_high = HACAalphaCalc(tempUpper, cellM0, cellM1);

        real ChiSoot_low = pdfChiSootCalc(tempLower, cellPressure,RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);
        real ChiSoot_high = pdfChiSootCalc(tempUpper, cellPressure,RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);

        real f_low_low = HACAIntegrandCalc(tempLower, cellTempAvg, cellTempRMS, C2H2Lower, C2H2MfAvg, C2H2MfRMS, alpha_low, ChiSoot_low);
        real P_low_low = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(C2H2Lower, C2H2MfAvg, C2H2MfRMS);

        real f_high_low = HACAIntegrandCalc(tempUpper, cellTempAvg, cellTempRMS, C2H2Lower, C2H2MfAvg, C2H2MfRMS, alpha_high, ChiSoot_high);
        real P_high_low = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(C2H2Lower, C2H2MfAvg, C2H2MfRMS);

        real f_low_high = HACAIntegrandCalc(tempLower, cellTempAvg, cellTempRMS, C2H2Upper, C2H2MfAvg, C2H2MfRMS, alpha_low, ChiSoot_low);
        real P_low_high = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(C2H2Upper, C2H2MfAvg, C2H2MfRMS);

        real f_high_high = HACAIntegrandCalc(tempUpper, cellTempAvg, cellTempRMS, C2H2Upper, C2H2MfAvg, C2H2MfRMS, alpha_high, ChiSoot_high);
        real P_high_high = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(C2H2Upper, C2H2MfAvg, C2H2MfRMS);

        real f_xi_low = 0.0; real P_xi_low = 0.0;
        real f_xi_high = 0.0; real P_xi_high = 0.0;
        real f_low_yi = 0.0; real P_low_yi = 0.0;
        real f_high_yi = 0.0; real P_high_yi = 0.0;
        real f_xi_yi = 0.0; real P_xi_yi = 0.0;

        int ctr;
        int ctr1;

        for (ctr = 1; ctr < 13; ctr++)
        {
            temp = tempLower + ctr * tempInt; C2H2 = C2H2Lower + ctr * C2H2Int;
            real alpha = HACAalphaCalc(temp, cellM0, cellM1);
            real ChiSoot = pdfChiSootCalc(temp, cellPressure,RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);
            f_xi_low = f_xi_low + HACAIntegrandCalc(temp, cellTempAvg, cellTempRMS, C2H2Lower, C2H2MfAvg, C2H2MfRMS, alpha, ChiSoot);
            P_xi_low = P_xi_low + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(C2H2Lower, C2H2MfAvg, C2H2MfRMS);

            f_xi_high = f_xi_high + HACAIntegrandCalc(temp, cellTempAvg, cellTempRMS, C2H2Upper, C2H2MfAvg, C2H2MfRMS, alpha, ChiSoot);
            P_xi_high = P_xi_high + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(C2H2Upper, C2H2MfAvg, C2H2MfRMS);

            f_low_yi = f_low_yi + HACAIntegrandCalc(tempLower, cellTempAvg, cellTempRMS, C2H2, C2H2MfAvg, C2H2MfRMS, alpha_low, ChiSoot_low);
            P_low_yi = P_low_yi + gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(C2H2, C2H2MfAvg, C2H2MfRMS);

            f_high_yi = f_high_yi + HACAIntegrandCalc(tempUpper, cellTempAvg, cellTempRMS, C2H2, C2H2MfAvg, C2H2MfRMS, alpha_high, ChiSoot_high);
            P_high_yi = P_high_yi + gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(C2H2, C2H2MfAvg, C2H2MfRMS);

            for (ctr1 = 1; ctr1 < 13; ctr1++)
            {
                C2H2_new = C2H2Lower + ctr1 * C2H2Int;
                f_xi_yi = f_xi_yi + HACAIntegrandCalc(temp, cellTempAvg, cellTempRMS, C2H2_new, C2H2MfAvg, C2H2MfRMS, alpha, ChiSoot);
                P_xi_yi = P_xi_yi + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(C2H2_new, C2H2MfAvg, C2H2MfRMS);
            }
        }
        pdfIntegrand = 0.25 * tempInt * C2H2Int * (P_low_high + P_high_high + P_high_low + P_low_low + \
            2 * (P_low_yi + P_high_yi + P_xi_low + P_xi_high) + 4 * P_xi_yi);

        integrand = 0.25 * tempInt * C2H2Int * (f_low_high + f_high_high + f_high_low + f_low_low + \
            2 * (f_low_yi + f_high_yi + f_xi_low + f_xi_high) + 4 * f_xi_yi);

        integrand = integrand / pdfIntegrand;
        C_UDMI(c,t,32) = 2.0;

    }

    C_UDMI(c,t,33) = pdfIntegrand;
    

    source = multConstant * integrand;

    real muTwoThirdDiff = fracMomFirstDiffM2(twoThird, cellM0, cellM1, cellM2);
    real muFiveThirdDiff = fracMomFirstDiffM2(5 * oneThird, cellM0, cellM1, cellM2);

    dS[eqn] = source / (4 * muTwoThird + 4 * muFiveThird) * (4 * muTwoThirdDiff + 4 * muFiveThirdDiff);

    /* dS[eqn] = A4 * cellPressure/(RGas*C2H2MW)*pow(Cs, 2) * pi * (4 * muTwoThirdDiff + 4 * muFiveThirdDiff) * integrand; */

    C_UDMI(c,t,3) = source;

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
    sourceOH = -gammaOH*OHMc*sqrt(pi*kB*cellTemp/(2*mOH))*pi*pow(Cs,2)*MtwoThird*avogad;
    source = sourceO2 + sourceOH;

    C_UDMI(c,t,5) = sourceO2;
    C_UDMI(c,t,6) = sourceOH;

    C_UDMI(c,t,7) = source*fudgeFac1;

    dS[eqn] = -2*k5*O2Mc*alpha*chiSoot*pi*pow(Cs,2)*pow(cellM0,twoThird/3)*pow(cellM2,-oneThird/3)*pow(cellM1,-oneThird/3)*8.0/9.0;
    return source*fudgeFac1;
}

DEFINE_SOURCE(pdf_m_1_OxSource,c,t,dS,eqn)
{
    real sourceOH = 0.0;
    real integrandOH = 0.0;
    real pdfIntegrandOH = 0.0;

    real sourceO2 = 0.0;
    real integrandO2 = 0.0;
    real pdfIntegrandO2 = 0.0;

    real temp;
    real OH;
    real OH_new;
    real O2;
    real O2_new;

    real RGas = C_RGAS(c,t) * 1e-3;
    real cellRho = C_R(c,t);
    real cellPressure = (C_P(c,t) + Patm)/1e3;

    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2);

    Material *mat = THREAD_MATERIAL(t);
    int idOH = mixture_specie_index(mat, "oh");
    int idO2 = mixture_specie_index(mat, "o2");

    
    real cellTempAvg = C_T(c,t);
    real cellTempRMS = C_UDMI(c,t,36);
    real TARMS = C_UDMI(c,t,31);
    real TAAvg = C_UDMI(c,t,30);
    real OHMfAvg = Pdf_Yi(c,t,idOH);
    real OHMfRMS = C_UDMI(c,t,39);
    real O2MfAvg = Pdf_Yi(c,t,idO2);
    real O2MfRMS = C_UDMI(c,t,40);

    real ta_max = C_UDMI(c,t,41);

    int idH2 = mixture_specie_index(mat, "h2");
    int idH = mixture_specie_index(mat, "h");
    int idH2O = mixture_specie_index(mat, "h2o");
    int idC2H2 = mixture_specie_index(mat, "c2h2");

    real H2Mf = Pdf_Yi(c,t,idH2);
    real HMf = Pdf_Yi(c,t,idH);
    real O2Mf = Pdf_Yi(c,t,idO2);
    real H2OMf = Pdf_Yi(c,t,idH2O);
    real C2H2Mf = Pdf_Yi(c,t,idC2H2);

    real Cs = pow((6*mC/(pi*rhoSoot)),oneThird);
    real muTwoThird = fractionalMoments(twoThird, cellM0, cellM1, cellM2); 

    real OHmultConstant =  -gammaOH * cellPressure * pow(Cs, 2) * sqrt(pi * kB/(2 * mOH)) * cellM0 * muTwoThird /(OHMW * RGas);
    real O2multConstant = -2 * A5 * O2Mf * cellPressure * pi * pow(Cs, 2) * muTwoThird * cellM0 / RGas;

    real tempRatio = cellTempRMS/(cellTempAvg);
    real OHMfRatio = OHMfRMS/(OHMfAvg + 1e-8);
    real TARatio = TARMS/TAAvg;
    real O2MfRatio = O2MfRMS / O2MfAvg;

    C_UDMI(c,t,33) = tempRatio;
    C_UDMI(c,t,34) = OHMfRatio;
    C_UDMI(c,t,35) = TARatio;

    if (tempRatio < 0.02 && OHMfRatio < 0.02)
    {

        integrandOH = OHMfAvg /sqrt(cellTempAvg);
        pdfIntegrandOH = 1.0;

    }

    else if (tempRatio >= 0.02 && OHMfRatio < 0.02)
    {

        real tempLower;

        real tempUpper = cellTempAvg + 3.5 * cellTempRMS;

        if (cellTempAvg - 3.5 * cellTempRMS > 300) { tempLower = cellTempAvg - 3.5 * cellTempRMS; }
        else { tempLower = 300; }

        real tempInt = 0.5 * cellTempRMS;


        real f_low = 1/sqrt(tempLower) * gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS);
        real P_low = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS);

        real f_high = 1/sqrt(tempUpper) * gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS);
        real P_high = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS);

        real f_xi = 0.0; real P_xi = 0.0;
        

        int ctr;

        for (ctr = 1; ctr < 13; ctr ++)
        {
            temp = tempLower + ctr * tempInt;

            f_xi = f_xi + 1/sqrt(temp) * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
            P_xi = P_xi + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
        }

        pdfIntegrandOH = (P_low + P_high + 2 * P_xi) * tempInt / 2.0;
        integrandOH = OHMfAvg * (f_low + f_high + 2 * f_xi) / pdfIntegrandOH;
        
    }

    else if (tempRatio < 0.02 && OHMfRatio >= 0.02)
    {
        real OHLower; real OHUpper;

        if (OHMfAvg - 3.5 * OHMfRMS > 0) { OHLower = OHMfAvg - 3.5 * OHMfRMS; }
        else { OHLower = 0.0; }
        if (OHMfAvg + 3.5 * OHMfRMS < 1.0) { OHUpper = OHMfAvg + 3.5 * OHMfRMS; }
        else { OHUpper = 1.0; }

        real OHInt = 0.5 * OHMfRMS;

        int ctr;

        real f_low = OHLower * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);
        real P_low = gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);

        real f_high = OHUpper * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);
        real P_high = gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);

        real f_yi = 0.0; real P_yi = 0.0;

        for (ctr = 1; ctr < 13; ctr ++)
        {
            OH = OHLower + OHInt * ctr;
            f_yi = f_yi + OH * gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);
            P_yi = P_yi + gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);
        }

        pdfIntegrandOH = (P_high + P_low + 2 * P_yi) * OHInt/2.0;
        integrandOH = 1 / sqrt(cellTempAvg) * (f_low + f_high + 2 * f_yi) / pdfIntegrandOH;

    }

    else 
    {
        real tempUpper; real tempLower;
        real OHUpper; real OHLower;

        tempUpper = cellTempAvg + 3.5 * cellTempRMS;
        if (cellTempAvg - 3.5 * cellTempRMS > 300) { tempLower = cellTempAvg - 3.5 * cellTempRMS; }
        else { tempLower = 300; }

        if (OHMfAvg - 3.5 * OHMfRMS > 0) { OHLower = OHMfAvg - 3.5 * OHMfRMS; }
        else { OHLower = 0.0; }

        if (OHMfAvg + 3.5 * OHMfRMS < 1.0) { OHUpper = OHMfAvg + 3.5 * OHMfRMS; }
        else { OHUpper = 1.0; }

        real tempInt = 0.5 * cellTempRMS; real OHInt = 0.5 * OHMfRMS;


        real f_low_low = OHLower / sqrt(tempLower);
        real P_low_low = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);

        real f_high_low = OHLower / sqrt(tempUpper);
        real P_high_low = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);

        real f_low_high = OHUpper / sqrt(tempLower);
        real P_low_high = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);

        real f_high_high = OHUpper / sqrt(tempUpper);
        real P_high_high = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);

        real f_xi_low = 0.0; real P_xi_low = 0.0;
        real f_xi_high = 0.0; real P_xi_high = 0.0;
        real f_low_yi = 0.0; real P_low_yi = 0.0;
        real f_high_yi = 0.0; real P_high_yi = 0.0;
        real f_xi_yi = 0.0; real P_xi_yi = 0.0;

        int ctr;
        int ctr1;

        for (ctr = 1; ctr < 13; ctr++)
        {
            temp = tempLower + ctr * tempInt; OH = OHLower + ctr * OHInt;

            f_xi_low = f_xi_low + OHLower / sqrt(temp);
            P_xi_low = P_xi_low + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);

            f_xi_high = f_xi_high + OHUpper / sqrt(temp);
            P_xi_high = P_xi_high + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);

            f_low_yi = f_low_yi + OH / sqrt(tempLower);
            P_low_yi = P_low_yi + gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            f_high_yi = f_high_yi + OH / sqrt(tempUpper);
            P_high_yi = P_high_yi + gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            for (ctr1 = 1; ctr1 < 13; ctr1++)
            {
                OH_new = OHLower + ctr1 * OHInt;
                f_xi_yi = f_xi_yi + OH / sqrt(temp);
                P_xi_yi = P_xi_yi + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OH_new, OHMfAvg, OHMfRMS);
            }
        }
        pdfIntegrandOH = 0.25 * tempInt * OHInt * (P_low_high + P_high_high + P_high_low + P_low_low + \
            2 * (P_low_yi + P_high_yi + P_xi_low + P_xi_high) + 4 * P_xi_yi);

        integrandOH = 0.25 * tempInt * OHInt * (f_low_high + f_high_high + f_high_low + f_low_low + \
            2 * (f_low_yi + f_high_yi + f_xi_low + f_xi_high) + 4 * f_xi_yi);

        integrandOH = integrandOH / pdfIntegrandOH;

    }

    /* Start pdf O2 oxidation  */

    if (tempRatio < 0.02 && O2MfRatio < 0.02)
    {
        real alpha = HACAalphaCalc(cellTempAvg, cellM0, cellM1);
        real chiSoot = pdfChiSootCalc(cellTempAvg, cellPressure, RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
        integrandO2 = 1 / cellTempAvg * exp(-Ea5/(R*cellTempAvg)) * alpha * chiSoot * O2MfAvg;
        pdfIntegrandO2 = 1.0;

    }

    else if (tempRatio >= 0.02 && O2MfRatio < 0.02)
    {

        real tempLower;

        real tempUpper = cellTempAvg + 3.5 * cellTempRMS;

        if (cellTempAvg - 3.5 * cellTempRMS > 300) { tempLower = cellTempAvg - 3.5 * cellTempRMS; }
        else { tempLower = 300; }

        real tempInt = 0.5 * cellTempRMS;

        real alpha_low = HACAalphaCalc(tempLower, cellM0, cellM1);
        real alpha_high = HACAalphaCalc(tempUpper, cellM0, cellM1);

        real ChiSoot_low = pdfChiSootCalc(tempLower, cellPressure, RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
        real ChiSoot_high = pdfChiSootCalc(tempUpper, cellPressure, RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);

        real f_low = 1 / tempLower * exp(-Ea5/(R*tempLower)) * alpha_low * ChiSoot_low * \
        gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS);

        real P_low = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS);

        real f_high = 1 / tempLower * exp(-Ea5/(R*tempUpper)) * alpha_high * ChiSoot_high * \
        gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS);

        real P_high = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS);

        real f_xi = 0.0; real P_xi = 0.0;
        

        int ctr;

        for (ctr = 1; ctr < 13; ctr ++)
        {
            temp = tempLower + ctr * tempInt;
            real alpha = HACAalphaCalc(temp, cellM0, cellM1);
            real ChiSoot = pdfChiSootCalc(temp, cellPressure,RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);

            f_xi = f_xi + pow(temp, - 1.0) * exp(-Ea5/(R*temp)) * alpha * ChiSoot * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
            P_xi = P_xi + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
        }

        pdfIntegrandO2 = (P_low + P_high + 2 * P_xi) * tempInt / 2.0;
        integrandO2 = O2MfAvg * (f_low + f_high + 2 * f_xi) / pdfIntegrandO2;
        
    }

    else if (tempRatio < 0.02 && O2MfRatio >= 0.02)
    {
        real alpha = HACAalphaCalc(cellTempAvg, cellM0, cellM1);
        real chiSoot = pdfChiSootCalc(cellTempAvg, cellPressure, RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
        real O2Lower; real O2Upper;

        if (O2MfAvg - 3.5 * O2MfRMS > 0) { O2Lower = O2MfAvg - 3.5 * O2MfRMS; }
        else { O2Lower = 0.0; }
        if (O2MfAvg + 3.5 * O2MfRMS < 1.0) { O2Upper = O2MfAvg + 3.5 * O2MfRMS; }
        else { O2Upper = 1.0; }

        real O2Int = 0.5 * O2MfRMS;

        int ctr;

        real f_low = O2Lower * gaussianMonoPDFCalc(O2Lower, O2MfAvg, O2MfRMS);
        real P_low = gaussianMonoPDFCalc(O2Lower, O2MfAvg, O2MfRMS);

        real f_high = O2Upper * gaussianMonoPDFCalc(O2Upper, O2MfAvg, O2MfRMS);
        real P_high = gaussianMonoPDFCalc(O2Upper, O2MfAvg, O2MfRMS);

        real f_yi = 0.0; real P_yi = 0.0;

        for (ctr = 1; ctr < 13; ctr ++)
        {
            O2 = O2Lower + O2Int * ctr;
            f_yi = f_yi + O2 * gaussianMonoPDFCalc(O2, O2MfAvg, O2MfRMS);
            P_yi = P_yi + gaussianMonoPDFCalc(O2, O2MfAvg, O2MfRMS);
        }

        pdfIntegrandO2 = (P_high + P_low + 2 * P_yi) * O2Int/2.0;
        integrandO2 = pow(cellTempAvg, - 1.0) * exp(-Ea5/(R*cellTempAvg)) * alpha * chiSoot * (f_low + f_high + 2 * f_yi) / pdfIntegrandO2;
    }


    else 
    {
        real tempUpper; real tempLower;
        real O2Upper; real O2Lower;

        tempUpper = cellTempAvg + 3.5 * cellTempRMS;
        if (cellTempAvg - 3.5 * cellTempRMS > 300) { tempLower = cellTempAvg - 3.5 * cellTempRMS; }
        else { tempLower = 300; }

        if (O2MfAvg - 3.5 * O2MfRMS > 0) { O2Lower = O2MfAvg - 3.5 * O2MfRMS; }
        else { O2Lower = 0.0; }

        if (O2MfAvg + 3.5 * O2MfRMS < 1.0) { O2Upper = O2MfAvg + 3.5 * O2MfRMS; }
        else { O2Upper = 1.0; }

        real tempInt = 0.5 * cellTempRMS; real O2Int = 0.5 * O2MfRMS;

        real alpha_low = HACAalphaCalc(tempLower, cellM0, cellM1);
        real alpha_high = HACAalphaCalc(tempUpper, cellM0, cellM1);

        real ChiSoot_low = pdfChiSootCalc(tempLower, cellPressure,RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
        real ChiSoot_high = pdfChiSootCalc(tempUpper, cellPressure,RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);

        real f_low_low = O2IntegrandCalc(tempLower, cellTempAvg, cellTempRMS, O2Lower, O2MfAvg, O2MfRMS, alpha_low, ChiSoot_low);
        real P_low_low = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(O2Lower, O2MfAvg, O2MfRMS);

        real f_high_low = O2IntegrandCalc(tempUpper, cellTempAvg, cellTempRMS, O2Lower, O2MfAvg, O2MfRMS, alpha_high, ChiSoot_high);
        real P_high_low = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(O2Lower, O2MfAvg, O2MfRMS);

        real f_low_high = O2IntegrandCalc(tempLower, cellTempAvg, cellTempRMS, O2Upper, O2MfAvg, O2MfRMS, alpha_low, ChiSoot_low);
        real P_low_high = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(O2Upper, O2MfAvg, O2MfRMS);

        real f_high_high = O2IntegrandCalc(tempUpper, cellTempAvg, cellTempRMS, O2Upper, O2MfAvg, O2MfRMS, alpha_high, ChiSoot_high);
        real P_high_high = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(O2Upper, O2MfAvg, O2MfRMS);

        real f_xi_low = 0.0; real P_xi_low = 0.0;
        real f_xi_high = 0.0; real P_xi_high = 0.0;
        real f_low_yi = 0.0; real P_low_yi = 0.0;
        real f_high_yi = 0.0; real P_high_yi = 0.0;
        real f_xi_yi = 0.0; real P_xi_yi = 0.0;

        int ctr;
        int ctr1;

        for (ctr = 1; ctr < 13; ctr++)
        {
            temp = tempLower + ctr * tempInt; O2 = O2Lower + ctr * O2Int;
            real alpha = HACAalphaCalc(temp, cellM0, cellM1);
            real ChiSoot = pdfChiSootCalc(temp, cellPressure,RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
            f_xi_low = f_xi_low + O2IntegrandCalc(temp, cellTempAvg, cellTempRMS, O2Lower, O2MfAvg, O2MfRMS, alpha, ChiSoot);
            P_xi_low = P_xi_low + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(O2Lower, O2MfAvg, O2MfRMS);

            f_xi_high = f_xi_high + O2IntegrandCalc(temp, cellTempAvg, cellTempRMS, O2Upper, O2MfAvg, O2MfRMS, alpha, ChiSoot);
            P_xi_high = P_xi_high + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(O2Upper, O2MfAvg, O2MfRMS);

            f_low_yi = f_low_yi + O2IntegrandCalc(tempLower, cellTempAvg, cellTempRMS, O2, O2MfAvg, O2MfRMS, alpha_low, ChiSoot_low);
            P_low_yi = P_low_yi + gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(O2, O2MfAvg, O2MfRMS);

            f_high_yi = f_high_yi + O2IntegrandCalc(tempUpper, cellTempAvg, cellTempRMS, O2, O2MfAvg, O2MfRMS, alpha_high, ChiSoot_high);
            P_high_yi = P_high_yi + gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(O2, O2MfAvg, O2MfRMS);

            for (ctr1 = 1; ctr1 < 13; ctr1++)
            {
                O2_new = O2Lower + ctr1 * O2Int;
                f_xi_yi = f_xi_yi + O2IntegrandCalc(temp, cellTempAvg, cellTempRMS, O2_new, O2MfAvg, O2MfRMS, alpha, ChiSoot);
                P_xi_yi = P_xi_yi + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(O2_new, O2MfAvg, O2MfRMS);
            }
        }
        pdfIntegrandO2 = 0.25 * tempInt * O2Int * (P_low_high + P_high_high + P_high_low + P_low_low + \
            2 * (P_low_yi + P_high_yi + P_xi_low + P_xi_high) + 4 * P_xi_yi);

        integrandO2 = 0.25 * tempInt * O2Int * (f_low_high + f_high_high + f_high_low + f_low_low + \
            2 * (f_low_yi + f_high_yi + f_xi_low + f_xi_high) + 4 * f_xi_yi);

        integrandO2 = integrandO2 / pdfIntegrandO2;

    }

    sourceOH = OHmultConstant * integrandOH;
    sourceO2 = O2multConstant * integrandO2;

    real muTwoThirdDiff = fracMomFirstDiffM1(twoThird, cellM0, cellM1, cellM2);
    /*real OHdiff = -gammaOH * cellPressure * pow(Cs, 2) * sqrt(pi * kB/(2 * mOH)) * cellM0 * avogad * muTwoThirdDiff /(OHMW * RGas);
    real O2diff = -2 * A5 * O2Mf * cellPressure * pi * pow(Cs, 2) * muTwoThirdDiff * cellM0 / RGas; */

    real OHdiff = sourceOH / muTwoThird * muTwoThirdDiff;
    real O2diff = sourceO2 / muTwoThird * muTwoThirdDiff;

    dS[eqn] = OHdiff + O2diff;
    real source = sourceOH + sourceO2;

    C_UDMI(c,t,5) = sourceOH;
    C_UDMI(c,t,6) = sourceO2;
    C_UDMI(c,t,7) = source;

    return source;
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

    C_UDMI(c,t,8) = sourceO2;
    C_UDMI(c,t,9) = sourceOH;

    dS[eqn] = k5*O2Mc*alpha*chiSoot*pi*pow(Cs,2)*(-4*pow(cellM0,twoThird/3)*pow(cellM2,-11*oneThird/3)*\
        pow(cellM1,-8*oneThird/3)/9.0+4*5.0/9.0*pow(cellM0,-oneThird/3)*pow(cellM2,-4*oneThird/3)*pow(cellM1,5*oneThird/3));


    C_UDMI(c,t,10) = source*fudgeFac1;
    return source*fudgeFac1;
}

DEFINE_SOURCE(pdf_m_2_OxSource,c,t,dS,eqn)
{
    real sourceOH = 0.0;
    real integrandOH = 0.0;
    real pdfIntegrandOH = 0.0;

    real sourceO2 = 0.0;
    real integrandO2 = 0.0;
    real pdfIntegrandO2 = 0.0;

    real temp;
    real OH;
    real OH_new;
    real O2;
    real O2_new;

    real RGas = C_RGAS(c,t) * 1e-3;
    real cellRho = C_R(c,t);
    real cellPressure = (C_P(c,t) + Patm)/1e3;

    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2);

    Material *mat = THREAD_MATERIAL(t);
    int idOH = mixture_specie_index(mat, "oh");
    int idO2 = mixture_specie_index(mat, "o2");

    real cellTempAvg = C_T(c,t);
    real cellTempRMS = C_UDMI(c,t,36);
    real TARMS = C_UDMI(c,t,31);
    real TAAvg = C_UDMI(c,t,30);
    real OHMfAvg = Pdf_Yi(c,t,idOH);
    real OHMfRMS = C_UDMI(c,t,39);
    real O2MfAvg = Pdf_Yi(c,t,idO2);
    real O2MfRMS = C_UDMI(c,t,40);

    real ta_max = C_UDMI(c,t,41);

    int idH2 = mixture_specie_index(mat, "h2");
    int idH = mixture_specie_index(mat, "h");
    int idH2O = mixture_specie_index(mat, "h2o");
    int idC2H2 = mixture_specie_index(mat, "c2h2");

    real H2Mf = Pdf_Yi(c,t,idH2);
    real HMf = Pdf_Yi(c,t,idH);
    real O2Mf = Pdf_Yi(c,t,idO2);
    real H2OMf = Pdf_Yi(c,t,idH2O);
    real C2H2Mf = Pdf_Yi(c,t,idC2H2);

    real Cs = pow((6*mC/(pi*rhoSoot)),oneThird);
    real muTwoThird = fractionalMoments(twoThird, cellM0, cellM1, cellM2);
    real muFiveThird = fractionalMoments(5.0 * oneThird, cellM0, cellM1, cellM2); 

    real OHmultConstant =  gammaOH * cellPressure * pow(Cs, 2) * sqrt(pi * kB/(2 * mOH)) * cellM0 * \
    (4 * muTwoThird - 4 * muFiveThird) * /(OHMW * RGas);

    real O2multConstant = A5 * O2Mf * cellPressure * pi * pow(Cs, 2) * (4 * muTwoThird - 4 * muFiveThird) * cellM0 / RGas;

    real tempRatio = cellTempRMS/(cellTempAvg);
    real OHMfRatio = OHMfRMS/(OHMfAvg + 1e-8);
    real TARatio = TARMS/TAAvg;
    real O2MfRatio = O2MfRMS / O2MfAvg;

    C_UDMI(c,t,33) = tempRatio;
    C_UDMI(c,t,34) = OHMfRatio;
    C_UDMI(c,t,35) = TARatio;

    if (tempRatio < 0.02 && OHMfRatio < 0.02)
    {

        integrandOH = OHMfAvg /sqrt(cellTempAvg);
        pdfIntegrandOH = 1.0;

    }

    else if (tempRatio >= 0.02 && OHMfRatio < 0.02)
    {

        real tempLower;

        real tempUpper = cellTempAvg + 3.5 * cellTempRMS;

        if (cellTempAvg - 3.5 * cellTempRMS > 300) { tempLower = cellTempAvg - 3.5 * cellTempRMS; }
        else { tempLower = 300; }

        real tempInt = 0.5 * cellTempRMS;


        real f_low = 1/sqrt(tempLower) * gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS);
        real P_low = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS);

        real f_high = 1/sqrt(tempUpper) * gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS);
        real P_high = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS);

        real f_xi = 0.0; real P_xi = 0.0;
        

        int ctr;

        for (ctr = 1; ctr < 13; ctr ++)
        {
            temp = tempLower + ctr * tempInt;

            f_xi = f_xi + 1/sqrt(temp) * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
            P_xi = P_xi + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
        }

        pdfIntegrandOH = (P_low + P_high + 2 * P_xi) * tempInt / 2.0;
        integrandOH = OHMfAvg * (f_low + f_high + 2 * f_xi) / pdfIntegrandOH;
        
    }

    else if (tempRatio < 0.02 && OHMfRatio >= 0.02)
    {
        real OHLower; real OHUpper;

        if (OHMfAvg - 3.5 * OHMfRMS > 0) { OHLower = OHMfAvg - 3.5 * OHMfRMS; }
        else { OHLower = 0.0; }
        if (OHMfAvg + 3.5 * OHMfRMS < 1.0) { OHUpper = OHMfAvg + 3.5 * OHMfRMS; }
        else { OHUpper = 1.0; }

        real OHInt = 0.5 * OHMfRMS;

        int ctr;

        real f_low = OHLower * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);
        real P_low = gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);

        real f_high = OHUpper * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);
        real P_high = gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);

        real f_yi = 0.0; real P_yi = 0.0;

        for (ctr = 1; ctr < 13; ctr ++)
        {
            OH = OHLower + OHInt * ctr;
            f_yi = f_yi + OH * gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);
            P_yi = P_yi + gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);
        }

        pdfIntegrandOH = (P_high + P_low + 2 * P_yi) * OHInt/2.0;
        integrandOH = 1 / sqrt(cellTempAvg) * (f_low + f_high + 2 * f_yi) / pdfIntegrandOH;

    }

    else 
    {
        real tempUpper; real tempLower;
        real OHUpper; real OHLower;

        tempUpper = cellTempAvg + 3.5 * cellTempRMS;
        if (cellTempAvg - 3.5 * cellTempRMS > 300) { tempLower = cellTempAvg - 3.5 * cellTempRMS; }
        else { tempLower = 300; }

        if (OHMfAvg - 3.5 * OHMfRMS > 0) { OHLower = OHMfAvg - 3.5 * OHMfRMS; }
        else { OHLower = 0.0; }

        if (OHMfAvg + 3.5 * OHMfRMS < 1.0) { OHUpper = OHMfAvg + 3.5 * OHMfRMS; }
        else { OHUpper = 1.0; }

        real tempInt = 0.5 * cellTempRMS; real OHInt = 0.5 * OHMfRMS;


        real f_low_low = OHLower / sqrt(tempLower);
        real P_low_low = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);

        real f_high_low = OHLower / sqrt(tempUpper);
        real P_high_low = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);

        real f_low_high = OHUpper / sqrt(tempLower);
        real P_low_high = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);

        real f_high_high = OHUpper / sqrt(tempUpper);
        real P_high_high = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);

        real f_xi_low = 0.0; real P_xi_low = 0.0;
        real f_xi_high = 0.0; real P_xi_high = 0.0;
        real f_low_yi = 0.0; real P_low_yi = 0.0;
        real f_high_yi = 0.0; real P_high_yi = 0.0;
        real f_xi_yi = 0.0; real P_xi_yi = 0.0;

        int ctr;
        int ctr1;

        for (ctr = 1; ctr < 13; ctr++)
        {
            temp = tempLower + ctr * tempInt; OH = OHLower + ctr * OHInt;

            f_xi_low = f_xi_low + OHLower / sqrt(temp);
            P_xi_low = P_xi_low + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);

            f_xi_high = f_xi_high + OHUpper / sqrt(temp);
            P_xi_high = P_xi_high + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);

            f_low_yi = f_low_yi + OH / sqrt(tempLower);
            P_low_yi = P_low_yi + gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            f_high_yi = f_high_yi + OH / sqrt(tempUpper);
            P_high_yi = P_high_yi + gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            for (ctr1 = 1; ctr1 < 13; ctr1++)
            {
                OH_new = OHLower + ctr1 * OHInt;
                f_xi_yi = f_xi_yi + OH / sqrt(temp);
                P_xi_yi = P_xi_yi + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OH_new, OHMfAvg, OHMfRMS);
            }
        }
        pdfIntegrandOH = 0.25 * tempInt * OHInt * (P_low_high + P_high_high + P_high_low + P_low_low + \
            2 * (P_low_yi + P_high_yi + P_xi_low + P_xi_high) + 4 * P_xi_yi);

        integrandOH = 0.25 * tempInt * OHInt * (f_low_high + f_high_high + f_high_low + f_low_low + \
            2 * (f_low_yi + f_high_yi + f_xi_low + f_xi_high) + 4 * f_xi_yi);

        integrandOH = integrandOH / pdfIntegrandOH;

    }

    /* Start pdf O2 oxidation  */

    if (tempRatio < 0.02 && O2MfRatio < 0.02)
    {
        real alpha = HACAalphaCalc(cellTempAvg, cellM0, cellM1);
        real chiSoot = pdfChiSootCalc(cellTempAvg, cellPressure, RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
        integrandO2 = 1 / cellTempAvg * exp(-Ea5/(R*cellTempAvg)) * alpha * chiSoot * O2MfAvg;
        pdfIntegrandO2 = 1.0;

    }

    else if (tempRatio >= 0.02 && O2MfRatio < 0.02)
    {

        real tempLower;

        real tempUpper = cellTempAvg + 3.5 * cellTempRMS;

        if (cellTempAvg - 3.5 * cellTempRMS > 300) { tempLower = cellTempAvg - 3.5 * cellTempRMS; }
        else { tempLower = 300; }

        real tempInt = 0.5 * cellTempRMS;

        real alpha_low = HACAalphaCalc(tempLower, cellM0, cellM1);
        real alpha_high = HACAalphaCalc(tempUpper, cellM0, cellM1);

        real ChiSoot_low = pdfChiSootCalc(tempLower, cellPressure, RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
        real ChiSoot_high = pdfChiSootCalc(tempUpper, cellPressure, RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);

        real f_low = 1 / tempLower * exp(-Ea5/(R*tempLower)) * alpha_low * ChiSoot_low * \
        gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS);

        real P_low = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS);

        real f_high = 1 / tempLower * exp(-Ea5/(R*tempUpper)) * alpha_high * ChiSoot_high * \
        gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS);

        real P_high = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS);

        real f_xi = 0.0; real P_xi = 0.0;
        

        int ctr;

        for (ctr = 1; ctr < 13; ctr ++)
        {
            temp = tempLower + ctr * tempInt;
            real alpha = HACAalphaCalc(temp, cellM0, cellM1);
            real ChiSoot = pdfChiSootCalc(temp, cellPressure,RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);

            f_xi = f_xi + pow(temp, - 1.0) * exp(-Ea5/(R*temp)) * alpha * ChiSoot * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
            P_xi = P_xi + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
        }

        pdfIntegrandO2 = (P_low + P_high + 2 * P_xi) * tempInt / 2.0;
        integrandO2 = O2MfAvg * (f_low + f_high + 2 * f_xi) / pdfIntegrandO2;
        
    }

    else if (tempRatio < 0.02 && O2MfRatio >= 0.02)
    {
        real alpha = HACAalphaCalc(cellTempAvg, cellM0, cellM1);
        real chiSoot = pdfChiSootCalc(cellTempAvg, cellPressure, RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
        real O2Lower; real O2Upper;

        if (O2MfAvg - 3.5 * O2MfRMS > 0) { O2Lower = O2MfAvg - 3.5 * O2MfRMS; }
        else { O2Lower = 0.0; }
        if (O2MfAvg + 3.5 * O2MfRMS < 1.0) { O2Upper = O2MfAvg + 3.5 * O2MfRMS; }
        else { O2Upper = 1.0; }

        real O2Int = 0.5 * O2MfRMS;

        int ctr;

        real f_low = O2Lower * gaussianMonoPDFCalc(O2Lower, O2MfAvg, O2MfRMS);
        real P_low = gaussianMonoPDFCalc(O2Lower, O2MfAvg, O2MfRMS);

        real f_high = O2Upper * gaussianMonoPDFCalc(O2Upper, O2MfAvg, O2MfRMS);
        real P_high = gaussianMonoPDFCalc(O2Upper, O2MfAvg, O2MfRMS);

        real f_yi = 0.0; real P_yi = 0.0;

        for (ctr = 1; ctr < 13; ctr ++)
        {
            O2 = O2Lower + O2Int * ctr;
            f_yi = f_yi + O2 * gaussianMonoPDFCalc(O2, O2MfAvg, O2MfRMS);
            P_yi = P_yi + gaussianMonoPDFCalc(O2, O2MfAvg, O2MfRMS);
        }

        pdfIntegrandO2 = (P_high + P_low + 2 * P_yi) * O2Int/2.0;
        integrandO2 = pow(cellTempAvg, - 1.0) * exp(-Ea5/(R*cellTempAvg)) * alpha * chiSoot * (f_low + f_high + 2 * f_yi) / pdfIntegrandO2;
    }


    else 
    {
        real tempUpper; real tempLower;
        real O2Upper; real O2Lower;

        tempUpper = cellTempAvg + 3.5 * cellTempRMS;
        if (cellTempAvg - 3.5 * cellTempRMS > 300) { tempLower = cellTempAvg - 3.5 * cellTempRMS; }
        else { tempLower = 300; }

        if (O2MfAvg - 3.5 * O2MfRMS > 0) { O2Lower = O2MfAvg - 3.5 * O2MfRMS; }
        else { O2Lower = 0.0; }

        if (O2MfAvg + 3.5 * O2MfRMS < 1.0) { O2Upper = O2MfAvg + 3.5 * O2MfRMS; }
        else { O2Upper = 1.0; }

        real tempInt = 0.5 * cellTempRMS; real O2Int = 0.5 * O2MfRMS;

        real alpha_low = HACAalphaCalc(tempLower, cellM0, cellM1);
        real alpha_high = HACAalphaCalc(tempUpper, cellM0, cellM1);

        real ChiSoot_low = pdfChiSootCalc(tempLower, cellPressure,RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
        real ChiSoot_high = pdfChiSootCalc(tempUpper, cellPressure,RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);

        real f_low_low = O2IntegrandCalc(tempLower, cellTempAvg, cellTempRMS, O2Lower, O2MfAvg, O2MfRMS, alpha_low, ChiSoot_low);
        real P_low_low = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(O2Lower, O2MfAvg, O2MfRMS);

        real f_high_low = O2IntegrandCalc(tempUpper, cellTempAvg, cellTempRMS, O2Lower, O2MfAvg, O2MfRMS, alpha_high, ChiSoot_high);
        real P_high_low = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(O2Lower, O2MfAvg, O2MfRMS);

        real f_low_high = O2IntegrandCalc(tempLower, cellTempAvg, cellTempRMS, O2Upper, O2MfAvg, O2MfRMS, alpha_low, ChiSoot_low);
        real P_low_high = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(O2Upper, O2MfAvg, O2MfRMS);

        real f_high_high = O2IntegrandCalc(tempUpper, cellTempAvg, cellTempRMS, O2Upper, O2MfAvg, O2MfRMS, alpha_high, ChiSoot_high);
        real P_high_high = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(O2Upper, O2MfAvg, O2MfRMS);

        real f_xi_low = 0.0; real P_xi_low = 0.0;
        real f_xi_high = 0.0; real P_xi_high = 0.0;
        real f_low_yi = 0.0; real P_low_yi = 0.0;
        real f_high_yi = 0.0; real P_high_yi = 0.0;
        real f_xi_yi = 0.0; real P_xi_yi = 0.0;

        int ctr;
        int ctr1;

        for (ctr = 1; ctr < 13; ctr++)
        {
            temp = tempLower + ctr * tempInt; O2 = O2Lower + ctr * O2Int;
            real alpha = HACAalphaCalc(temp, cellM0, cellM1);
            real ChiSoot = pdfChiSootCalc(temp, cellPressure,RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
            f_xi_low = f_xi_low + O2IntegrandCalc(temp, cellTempAvg, cellTempRMS, O2Lower, O2MfAvg, O2MfRMS, alpha, ChiSoot);
            P_xi_low = P_xi_low + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(O2Lower, O2MfAvg, O2MfRMS);

            f_xi_high = f_xi_high + O2IntegrandCalc(temp, cellTempAvg, cellTempRMS, O2Upper, O2MfAvg, O2MfRMS, alpha, ChiSoot);
            P_xi_high = P_xi_high + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(O2Upper, O2MfAvg, O2MfRMS);

            f_low_yi = f_low_yi + O2IntegrandCalc(tempLower, cellTempAvg, cellTempRMS, O2, O2MfAvg, O2MfRMS, alpha_low, ChiSoot_low);
            P_low_yi = P_low_yi + gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(O2, O2MfAvg, O2MfRMS);

            f_high_yi = f_high_yi + O2IntegrandCalc(tempUpper, cellTempAvg, cellTempRMS, O2, O2MfAvg, O2MfRMS, alpha_high, ChiSoot_high);
            P_high_yi = P_high_yi + gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(O2, O2MfAvg, O2MfRMS);

            for (ctr1 = 1; ctr1 < 13; ctr1++)
            {
                O2_new = O2Lower + ctr1 * O2Int;
                f_xi_yi = f_xi_yi + O2IntegrandCalc(temp, cellTempAvg, cellTempRMS, O2_new, O2MfAvg, O2MfRMS, alpha, ChiSoot);
                P_xi_yi = P_xi_yi + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(O2_new, O2MfAvg, O2MfRMS);
            }
        }
        pdfIntegrandO2 = 0.25 * tempInt * O2Int * (P_low_high + P_high_high + P_high_low + P_low_low + \
            2 * (P_low_yi + P_high_yi + P_xi_low + P_xi_high) + 4 * P_xi_yi);

        integrandO2 = 0.25 * tempInt * O2Int * (f_low_high + f_high_high + f_high_low + f_low_low + \
            2 * (f_low_yi + f_high_yi + f_xi_low + f_xi_high) + 4 * f_xi_yi);

        integrandO2 = integrandO2 / pdfIntegrandO2;

    }

    sourceOH = OHmultConstant * integrandOH;
    sourceO2 = O2multConstant * integrandO2;

    real muTwoThirdDiff = fracMomFirstDiffM2(twoThird, cellM0, cellM1, cellM2);
    real muFiveThirdDiff = fracMomFirstDiffM2(5.0 * oneThird, cellM0, cellM1, cellM2);
    /*real OHdiff = -gammaOH * cellPressure * pow(Cs, 2) * sqrt(pi * kB/(2 * mOH)) * cellM0 * avogad * muTwoThirdDiff /(OHMW * RGas);
    real O2diff = -2 * A5 * O2Mf * cellPressure * pi * pow(Cs, 2) * muTwoThirdDiff * cellM0 / RGas; */

    real OHdiff = sourceOH / (4 * muTwoThird - 4 * muFiveThird) * (4 * muTwoThirdDiff - 4 * muFiveThirdDiff);
    real O2diff = sourceO2 / (4 * muTwoThird - 4 * muFiveThird) * (4 * muTwoThirdDiff - 4 * muFiveThirdDiff);

    dS[eqn] = OHdiff + O2diff;
    real source = sourceOH + sourceO2;

    C_UDMI(c,t,5) = sourceOH;
    C_UDMI(c,t,6) = sourceO2;
    C_UDMI(c,t,7) = source;

    return source;
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
    real f_2_0_0 = 2*(muMinusOneHalf*muThirteenSixth + 2*muMinusOneSixth*muElevenSixth + muOneSixth*muThreeHalf + \
        2*muOneHalf*muSevenSixth + 2*pow(muFiveSixth,2));

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

DEFINE_SOURCE (new_m_0_CoagSource,c,t,dS,eqn)
{
    real source;
    real cellPressure = C_P(c,t) + Patm;
    real cellRho = C_R(c,t);
    real cellTemp = C_T(c,t);
    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2);
    real lamVisc = C_MU_L(c,t);


    real muOneThird = fractionalMoments(oneThird, cellM0, cellM1, cellM2);

    real meanFreePath = kB*cellTemp/(sqrt(2)*pi*pow(dAir,2)*cellPressure);
    real avgDia = muOneThird * dC; 
    real Kn = 2*meanFreePath/avgDia;

    C_UDMI(c,t,20) = meanFreePath; 
    C_UDMI(c,t,21) = avgDia;
    C_UDMI(c,t,22) = Kn;

    real Kc = 2*kB*cellTemp/(3*lamVisc);
    real KcPrime = 2.514*meanFreePath*pow(pi*rhoSoot/6,oneThird);
    real Kf = vDW*sqrt(6*kB*cellTemp/rhoSoot)*pow(3*mC/(4*pi*rhoSoot),oneSixth);
    
    /*real muTwoThird = pow(10, lagInterp3mom(twoThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));*/
    real muTwoThird = fractionalMoments(twoThird, cellM0, cellM1, cellM2);
    /*real muFourThird = pow(10, lagInterp3mom(4*oneThird, log10(cellM0/cellM0), log10(cellM1/cellM0), log10(cellM2/cellM0)));*/
    real muFourThird = fractionalMoments(4*oneThird, cellM0, cellM1, cellM2);
    real muFiveThird = fractionalMoments(5*oneThird, cellM0, cellM1, cellM2);

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
    real f_2_0_0 = 2*(muMinusOneHalf*muThirteenSixth + 2*muMinusOneSixth*muElevenSixth + muOneSixth*muThreeHalf + \
        2*muOneHalf*muSevenSixth + 2*pow(muFiveSixth,2));

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
    real f_2_1_1 = 2*muNineteenSixth*muOneHalf + 4*muSeventeenSixth*muFiveSixth + 2*muFiveHalf*muSevenSixth + \
    4*muThirteenSixth*muThreeHalf + 4*muElevenSixth*muElevenSixth;

    real f_oneHalf_1_1 = pow(10, lagInterp3mom(oneHalf, log10(f_0_1_1), log10(f_1_1_1), log10(f_2_1_1)));

    real Gc = Kc*(2*cellM1*cellM1/(cellM0*cellM0) + 2*muFourThird*muTwoThird + KcPrime*(2*muTwoThird*muOneThird + \
        2*muFourThird*muOneThird))*pow(cellM0,2)*normParameter;

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

real pdfpyrNucSource(cell_t c,Thread *t)
{
    Material *mat = THREAD_MATERIAL(t); 
    int idPyr = mixture_specie_index(mat, "a4"); 
    real pyrMf = Pdf_Yi(c,t,idPyr);
    real pyrMfMean = Pdf_Yi(c,t,idPyr);
    real CNucPyr = vDW*sqrt(4*pi*kB/(mC*NPyr))*pow(dPyr*avogad,2);
    real cellPressure = (C_P(c,t) + Patm)/1e3;
    real multConstant = gammaPyr*CNucPyr*pow(cellPressure/(C_RGAS(c,t)*1e-3*pyrMW),2);
    
    real tempMean = C_T(c,t);
    real tempVar = C_UDMI(c,t,36);
    real pyrMfVar = C_UDMI(c,t,37); 
    /*real tempVar = 50.0;
    real pyrMfVar = 1e-6;*/

    real tempPDFIntegral = 0.0;
    real tempOnlyPDFIntegral = 0.0;
    real pyrPDFIntegral = 0.0;
    real pyrOnlyPDFIntegral = 0.0;
    real lowerLim = 0.0;
    real upperLim = 0.0;
    int ctr = 1;
    
    if (tempVar/tempMean > 0.02)
    {
        if (tempMean - 3.5*tempVar > 300) { lowerLim = tempMean - 3.5*tempVar; }
        else { lowerLim = 300; }

        upperLim = tempMean + 4.0*tempVar;

        for (ctr = 1 ; ctr < 15; ctr = ctr + 1)
        {
            real temp = lowerLim + 0.5*ctr*tempVar;
            real temp_old = lowerLim + 0.5*(ctr-1.0)*tempVar;
            real temp0 = pow(temp_old,-3*oneHalf)*1/(sqrt(2.0*pi)*tempVar)*exp(-(pow(temp_old - tempMean, 2))/(2.0*tempVar*tempVar));
            
            real temp1 = pow(temp,-3*oneHalf)*1/(sqrt(2.0*pi)*tempVar)*exp(-(pow(temp - tempMean, 2))/(2.0*tempVar*tempVar));
            
            real tempdf0 = 1/(sqrt(2.0*pi)*tempVar)*exp(-(pow(temp_old - tempMean, 2))/(2.0*tempVar*tempVar));
            real tempdf1 = 1/(sqrt(2.0*pi)*tempVar)*exp(-(pow(temp - tempMean, 2))/(2.0*tempVar*tempVar));

            tempPDFIntegral = tempPDFIntegral + 0.5*(temp0 + temp1)*0.5*tempVar;
            
            tempOnlyPDFIntegral = tempOnlyPDFIntegral + 0.5*(tempdf0 + tempdf1)*0.5*tempVar;
            
        }

        tempPDFIntegral = tempPDFIntegral/(tempOnlyPDFIntegral);
    }
    
    if (tempVar/tempMean <= 0.02)
    {
        tempPDFIntegral = pow(tempMean, -3*oneHalf);
    }

    
    
    if (pyrMfVar/pyrMfMean > 0.02)
    {
        if (pyrMfMean - 3.5*pyrMfVar > 0) { lowerLim = pyrMfMean - 3.5*pyrMfVar; }
        else { lowerLim = 0; }

        if (pyrMfMean + 3.5*pyrMfVar < 1) { upperLim = pyrMfMean + 4*pyrMfVar; }
        else { upperLim = 1.0; }

        for (ctr = 1; ctr < 15; ctr = ctr + 1)
        {
            real pyr = lowerLim + 0.5*ctr*pyrMfVar;
            real pyr0 = pow(pyr - 0.5*pyrMfVar, 2)*1/(sqrt(2.0*pi)*pyrMfVar)*exp(-(pow(pyr - 0.5*pyrMfVar - pyrMfMean, 2))/(2.0*pyrMfVar*pyrMfVar));
            real pyr1 = pow(pyr, 2)*1/(sqrt(2.0*pi)*pyrMfVar)*exp(-(pow(pyr - pyrMfMean, 2))/(2.0*pyrMfVar*pyrMfVar));
            real pyrpdf0 = 1/(sqrt(2.0*pi)*pyrMfVar)*exp(-(pow(pyr - 0.5*pyrMfVar - pyrMfMean, 2))/(2.0*pyrMfVar*pyrMfVar));
            real pyrpdf1 = 1/(sqrt(2.0*pi)*pyrMfVar)*exp(-(pow(pyr - pyrMfMean, 2))/(2.0*pyrMfVar*pyrMfVar));

            pyrPDFIntegral = pyrPDFIntegral + 0.5*(pyr0 + pyr1)*0.5*pyrMfVar;
            pyrOnlyPDFIntegral = pyrOnlyPDFIntegral + 0.5*(pyrpdf0 + pyrpdf1)*0.5*pyrMfVar;
        }

        pyrPDFIntegral = pyrPDFIntegral/pyrOnlyPDFIntegral;
    
    }

    if (pyrMfVar/pyrMfMean <= 0.02)
    {
        pyrPDFIntegral = pyrMf*pyrMf;
    }
    
    real sourcePyrNuc = multConstant*tempPDFIntegral*pyrPDFIntegral/normParameter;
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
    real k2f = A2f*pow(cellTemp, n2f)*exp(-Ea2f/(R*cellTemp));
    real k2r = A2r*pow(cellTemp, n2r)*exp(-Ea2r/(R*cellTemp));
    real k3 = A3; 
    real k4 = A4*pow(cellTemp, n4)*exp(-Ea4/(R*cellTemp));
    real k5 = A5*exp(-Ea5/(R*cellTemp));

    real numerator = k1f*HMc + k2f*OHMc;
    real denominator = k1r*H2Mc + k2r*H2OMc + k3*HMc + k4*C2H2Mc + k5*O2Mc;

    real chiSoot = numerator/(denominator + 1.0) * siteDensity;

    return chiSoot;
}

real HACAalphaCalc(real Temp, real cellM0, real cellM1)
{
    real a = a1 + a2*Temp;
    real b = b1 + b2*Temp;
    real alpha = tanh(a/log10(cellM1/cellM0)+b); 
    if (alpha < 0 ) { alpha = 0.001; }

    return alpha;
}

real alphaOxCalculator(real ta, real ta_max)
{
    real retval = pow(ta_max/ta, 2) * exp(2*(1 - ta_max/ta));
    return retval;
}

real pdfChiSootCalc(real Temp, real Pstat, real RGas, real HMf, real OHMf, real H2Mf, real H2OMf, real C2H2Mf, real O2Mf)
{
    real Rho = Pstat/(RGas * Temp);

    real C2H2Mc = C2H2Mf*Rho/C2H2MW;
    real H2Mc = H2Mf*Rho/H2MW;
    real HMc = HMf*Rho/HMW;
    real O2Mc = O2Mf*Rho/O2MW;
    real OHMc = OHMf*Rho/OHMW;
    real H2OMc = H2OMf*Rho/H2OMW;

    real k1f = A1f*exp(-Ea1f/(R*Temp));
    real k1r = A1r*exp(-Ea1r/(R*Temp));
    real k2f = A2f*pow(Temp, n2f)*exp(-Ea2f/(R*Temp));
    real k2r = A2r*pow(Temp, n2r)*exp(-Ea2r/(R*Temp));
    real k3 = A3; 
    real k4 = A4*pow(Temp, n4)*exp(-Ea4/(R*Temp));
    real k5 = A5*exp(-Ea5/(R*Temp));

    real numerator = k1f*HMc + k2f*OHMc;
    real denominator = k1r*H2Mc + k2r*H2OMc + k3*HMc + k4*C2H2Mc + k5*O2Mc;

    real chiSoot = numerator/(denominator + 1.0) * siteDensity;

    return chiSoot;
}

real GcCalc(r, Kc, KcPrime, cellM0, cellM1, cellM2)
{
    int k;
    real retval = 0;
    for (k = 1; k < r ; k++)
    {
        real muk = fractionalMoments(k, cellM0, cellM1, cellM2);
        real murMinusk = fractionalMoments(r-k, cellM0, cellM1, cellM2);
        real mukPlusOneThird = fractionalMoments(k + oneThird, cellM0, cellM1, cellM2);
        real mukMinusOneThird = fractionalMoments(k - oneThird, cellM0, cellM1, cellM2);
        real murMinuskMinusOneThird = fractionalMoments(r - k - oneThird, cellM0, cellM1, cellM2);
        real murMinuskPlusOneThird = fractionalMoments(r - k + oneThird, cellM0, cellM1, cellM2);
        real murMinuskMinustwoThird = fractionalMoments(r - k - twoThird, cellM0, cellM1, cellM2);
        real mukMinusTwoThird = fractionalMoments(k - twoThird, cellM0, cellM1, cellM2);
        int comb = factorialCalc(r)/(factorialCalc(k) * factorialCalc(r-k));
        int KcPrimeTerm = KcPrime * (mukMinusOneThird*murMinusk + muk*murMinuskMinusOneThird + \
            mukPlusOneThird*murMinuskMinustwoThird + mukMinusTwoThird*murMinuskPlusOneThird);

        retval = retval + comb * (2*muk*murMinusk + mukPlusOneThird*murMinuskMinusOneThird + \
            mukMinusOneThird*murMinuskPlusOneThird + KcPrimeTerm);

    }
    retval = retval * 0.5 * Kc * cellM0 * cellM0;
    return retval;
}

real GcDiffM2Calc(r, Kc, KcPrime, cellM0, cellM1, cellM2)
{
    int k;
    real retval = 0;
    for (k = 1; k < r ; k++)
    {
        real muk = fracMomFirstDiffM2(k, cellM0, cellM1, cellM2);
        real murMinusk = fracMomFirstDiffM2(r-k, cellM0, cellM1, cellM2);
        real mukPlusOneThird = fracMomFirstDiffM2(k + oneThird, cellM0, cellM1, cellM2);
        real mukMinusOneThird = fracMomFirstDiffM2(k - oneThird, cellM0, cellM1, cellM2);
        real murMinuskMinusOneThird = fracMomFirstDiffM2(r - k - oneThird, cellM0, cellM1, cellM2);
        real murMinuskPlusOneThird = fracMomFirstDiffM2(r - k + oneThird, cellM0, cellM1, cellM2);
        real murMinuskMinustwoThird = fracMomFirstDiffM2(r - k - twoThird, cellM0, cellM1, cellM2);
        real mukMinusTwoThird = fracMomFirstDiffM2(k - twoThird, cellM0, cellM1, cellM2);
        int comb = factorialCalc(r)/(factorialCalc(k) * factorialCalc(r-k));
        int KcPrimeTerm = KcPrime * (mukMinusOneThird*murMinusk + muk*murMinuskMinusOneThird + \
            mukPlusOneThird*murMinuskMinustwoThird + mukMinusTwoThird*murMinuskPlusOneThird);

        retval = retval + comb * (2*muk*murMinusk + mukPlusOneThird*murMinuskMinusOneThird + \
            mukMinusOneThird*murMinuskPlusOneThird + KcPrimeTerm);
    }
    retval = retval * 0.5 * Kc * cellM0 * cellM0;
    return retval; 
}


real factorialCalc(n)
{
    int ctr;
    int retval = 1;
    for (ctr = n; ctr > 0; ctr --)
    {
        retval = retval * ctr;
    }
    return retval;
}

real HACAIntegrandCalc(real Temp, real TempAvg, real TempRMS, real C2H2Mf, real C2H2MfAvg, real C2H2MfRMS, real alpha, real ChiSoot)
{
    real P_Temp = gaussianMonoPDFCalc(Temp, TempAvg, TempRMS);
    real P_C2H2 = gaussianMonoPDFCalc(C2H2Mf, C2H2MfAvg, C2H2MfRMS);
    real integrand = pow(Temp, n4 - 1.0) * exp(-Ea4/(R*Temp)) * alpha * ChiSoot * C2H2Mf * P_Temp * P_C2H2;
    return integrand;
}

real O2IntegrandCalc(real Temp, real TempAvg, real TempRMS, real O2Mf, real O2MfAvg, real O2MfRMS, real alpha, real ChiSoot)
{
    real P_Temp = gaussianMonoPDFCalc(Temp, TempAvg, TempRMS);
    real P_O2 = gaussianMonoPDFCalc(O2Mf, O2MfAvg, O2MfRMS);
    real integrand = pow(Temp,- 1.0) * exp(-Ea5/(R*Temp)) * alpha * ChiSoot * O2Mf * P_Temp * P_O2;
    return integrand;
}

real fractionalMoments(real p, real m0, real m1, real m2)
{
    real fracMoment = pow(m0, 0.5*p*(p-3)) * pow(m1, -p*(p-2)) * pow(m2, 0.5*p*(p-1));
    return fracMoment;
}

real fracMomFirstDiffM0(real p, real m0, real m1, real m2)
{
    real fracMoment = 0.5 * p * (p-3) * pow(m0, 0.5*p*(p-3) - 1) * pow(m1, -p*(p-2)) * pow(m2, 0.5*p*(p-1));
    return fracMoment;
}

real fracMomFirstDiffM1(real p, real m0, real m1, real m2)
{
    real diff = -p * (p-2) * pow(m0, 0.5*p*(p-3)) * pow(m1, -p*(p-2)-1) * pow(m2, 0.5*p*(p-1));
    return diff;
}

real fracMomFirstDiffM2(real p, real m0, real m1, real m2)
{
    real diff = 0.5 * p * (p-1) * pow(m0, 0.5*p*(p-3)) * pow(m1, -p*(p-2)) * pow(m2, 0.5*p*(p-1) - 1);
    return diff;
}

real gaussianMonoPDFCalc(real Var, real VarMean, real VarStd)
{
    real value = 1/(VarStd * sqrt(2.0*pi)) * exp(-pow(Var - VarMean, 2)/(2.0 * VarStd * VarStd));
    return value;
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