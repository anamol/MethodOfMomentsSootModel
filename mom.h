
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

#ifndef mom_h
#define mom_h

real lagInterp3mom(real p, real m0, real m1, real m2); /* Lagragian interpolation */
real logInterp(real p, real a, real b, real M, real N); /* logarithmic interpolation */
real pyrNucSource(cell_t c,Thread *t);
real pdfpyrNucSource(cell_t c,Thread *t);
real benzNucSource(cell_t c,Thread *t);
real chiSootCalc(cell_t c,Thread *t, real cellTemp, real cellRho);
real pdfChiSootCalc(real Temp, real Pstat, real RGas, real HMf, real OHMf, real H2Mf, real H2OMf, real C2H2Mf, real O2Mf);
real HACAalphaCalc(real Temp, real cellM0, real cellM1);
real HACAIntegrandCalc(real Temp, real TempAvg, real TempRMS, real C2H2Mf, real C2H2MfAvg, real C2H2MfRMS, real alpha, real ChiSoot);
real O2IntegrandCalc(real Temp, real TempAvg, real TempRMS, real O2Mf, real O2MfAvg, real O2MfRMS, real alpha, real ChiSoot);
real fractionalMoments(real p, real m0, real m1, real m2);
real gaussianMonoPDFCalc(real Var, real VarMean, real VarStd);
real factorialCalc(n);
real alphaOxCalculator(real ta, real ta_max);
real GcCalc(r, Kc, KcPrime, cellM0, cellM1, cellM2);
real fracMomFirstDiffM0(real p, real m0, real m1, real m2);
real fracMomFirstDiffM1(real p, real m0, real m1, real m2);
real fracMomFirstDiffM2(real p, real m0, real m1, real m2);
real GcDiffM2Calc(r, Kc, KcPrime, cellM0, cellM1, cellM2);

#endif