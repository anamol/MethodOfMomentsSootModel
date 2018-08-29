#include "udf.h"
#include "mom.h"

real lagInterp3mom(real p, real m0, real m1, real m2)
{
    real mp = 0.5*(p-1)*(p-2)*m0 - p*(p-2)*m1 + 0.5*p*(p-1)*m2;

    return mp;

}

real logInterp(real p, real a, real b, real M, real N)
{
    real f = (p-a)/(b-a);
    return pow(N,f) * pow(M,1-f);
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

real HACAalphaCalc(real Temp, real cellM0, real cellM1)
{
    real a = a1 + a2*Temp;
    real b = b1 + b2*Temp;
    real alpha = tanh(a/log10(cellM1/cellM0)+b); 
    if (alpha < 0 ) { alpha = 0.001; }

    return alpha;
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

real gaussianMonoPDFCalc(real Var, real VarMean, real VarStd)
{
    real value = 1/(VarStd * sqrt(2.0*pi)) * exp(-pow(Var - VarMean, 2)/(2.0 * VarStd * VarStd));
    return value;
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

real alphaOxCalculator(real ta, real ta_max)
{
    real retval = pow(ta_max/ta, 2) * exp(2*(1 - ta_max/ta));
    return retval;
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


