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

real O2IntegrandCalc_TApdf(real Temp, real TempAvg, real TempRMS, real TA, real TAAvg, real TARMS, real ta_max, real ChiSoot)
{
    real P_Temp = gaussianMonoPDFCalc(Temp, TempAvg, TempRMS);
    real P_TA = gaussianMonoPDFCalc(TA, TAAvg, TARMS);
    real integrand = pow(Temp,- 1.0) * exp(-Ea5/(R*Temp)) * alphaOxCalculator(TA, ta_max) * ChiSoot * P_Temp * P_TA;
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
    if (ta == 0.0)
    {
        ta = 1e-5;
    }
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

real CoagDiffM0(m0,m1,m2,Kf)
{
    real diff = -(Kf*pow(m0,39.0/32)*pow(m2,17.0/96) * ( \
    262*pow(m1,58.0/9) + 181*pow(m0,32.0/9)*pow(m2,26.0/9) + 244*m0*pow(m1,40.0/9)*m2 + \
    572*pow(m0,2.0/9)*m1^6*pow(m2,2.0/9) + 254*m0^2*pow(m1,22.0/9)*m2^2 + \
    536*pow(m0,11.0/9)*m1^4*pow(m2,11.0/9) + 157*pow(m0,14.0/9)*m1^4*pow(m2,8.0/9) + \
    656*pow(m0,7.0/3)*pow(m1,22.0/9)*pow(m2,5.0/3) + 320*pow(m0,11.0/9)*pow(m1,14.0/3)*pow(m2,5.0/9) + \
    316*pow(m0,5.0/3)*pow(m1,28.0/9)*pow(m2,5.0/3) + 392*pow(m0,17.0/9)*pow(m1,10.0/3)*pow(m2,11.0/9) + \
    692*pow(m0,4.0/3)*pow(m1,40.0/9)*pow(m2,2.0/3) + 193*pow(m0,20.0/9)*pow(m1,8.0/3)*pow(m2,14.0/9) + \
    612*pow(m0,2.0/3)*pow(m1,46.0/9)*pow(m2,2.0/3) + 568*pow(m0,1.0/3)*pow(m1,52.0/9)*pow(m2,1.0/3) + \
    145*pow(m0,26.0/9)*pow(m1,4.0/3)*pow(m2,20.0/9) + 260*pow(m0,1.0/9)*pow(m1,56.0/9)*pow(m2,1.0/9) + \
    127*pow(m0,4.0/9)*pow(m1,50.0/9)*pow(m2,4.0/9) + 332*pow(m0,7.0/9)*pow(m1,44.0/9)*pow(m2,7.0/9) + \
    163*pow(m0,10.0/9)*pow(m1,38.0/9)*pow(m2,10.0/9) + 322*pow(m0,10.0/9)*pow(m1,44.0/9)*pow(m2,4.0/9) + \ 
    688*pow(m0,13.0/9)*pow(m1,38.0/9)*pow(m2,7.0/9) + 115*pow(m0,16.0/9)*pow(m1,26.0/9)*pow(m2,16.0/9) + \
    732*pow(m0,16.0/9)*pow(m1,32.0/9)*pow(m2,10.0/9) + 304*pow(m0,19.0/9)*pow(m1,26.0/9)*pow(m2,13.0/9) + \
    151*pow(m0,22.0/9)*pow(m1,14.0/9)*pow(m2,22.0/9) + 376*pow(m0,25.0/9)*pow(m1,14.0/9)*pow(m2,19.0/9) + \
    314*pow(m0,28.0/9)*pow(m1,8.0/9)*pow(m2,22.0/9))) /(72*pow(m1,17.0/18)*(pow(m0,7.0/4)*pow(m2,3.0/4) + \
    pow(m0,23.0/36)*pow(m1,14.0/9)*pow(m2,11.0/36))^(5.0/8)*(pow(m1,4.0/3) + pow(m0,2.0/3)*pow(m2,2.0/3) + \
    2*pow(m0,2.0/9)*pow(m1,8.0/9)*pow(m2,2.0/9))^(1.0/4)*(2*pow(m1,32.0/9) + pow(m0,16.0/9)*pow(m2,16.0/9) + \ 
    2*m0*pow(m1,14.0/9)*m2 + 2*pow(m0,1.0/9)*pow(m1,10.0/3)*pow(m2,1.0/9) + pow(m0,4.0/9)*pow(m1,8.0/3)*pow(m2,4.0/9))^(9.0/8));

    return diff;
}

real CoagDiffM2(m0,m1,m2,Kf)
{
    real diff = (Kf* ( \
    64*pow(m1,100.0/9) + 134*pow(m0,127.0/18)*pow(m2,73.0/18) - 8*m0*pow(m1,82.0/9)*m2 + \
    320*pow(m0,2.0/9)*pow(m1,32.0/3)*pow(m2,2.0/9) + 176*pow(m0,11.0/9)*pow(m1,26.0/3)*pow(m2,11.0/9) - \
    4*pow(m0,7.0/2)*pow(m1,46.0/9)*pow(m2,5.0/2) + 346*pow(m0,35.0/9)*pow(m1,16.0/3)*pow(m2,17.0/9) + \
    167*pow(m0,38.0/9)*pow(m1,14.0/3)*pow(m2,20.0/9) + 41*pow(m0,37.0/18)*m1^8*pow(m2,19.0/18) + \
    8*pow(m0,44.0/9)*pow(m1,10.0/3)*pow(m2,26.0/9) + 496*pow(m0,19.0/6)*pow(m1,52.0/9)*pow(m2,13.0/6) + \
    128*pow(m0,2.0/3)*pow(m1,88.0/9)*pow(m2,2.0/3) + 220*pow(m0,17.0/6)*pow(m1,58.0/9)*pow(m2,11.0/6) + \
    94*pow(m0,31.0/18)*pow(m1,26.0/3)*pow(m2,13.0/18) + 304*pow(m0,1.0/3)*pow(m1,94.0/9)*pow(m2,1.0/3) + \
    568*pow(m0,13.0/6)*pow(m1,70.0/9)*pow(m2,7.0/6) + 292*pow(m0,11.0/6)*pow(m1,76.0/9)*pow(m2,5.0/6) + \
    56*pow(m0,1.0/9)*pow(m1,98.0/9)*pow(m2,1.0/9) + 16*pow(m0,4.0/9)*pow(m1,92.0/9)*pow(m2,4.0/9) + \
    350*pow(m0,34.0/9)*pow(m1,50.0/9)*pow(m2,16.0/9) + 314*pow(m0,43.0/9)*pow(m1,32.0/9)*pow(m2,25.0/9) + \
    64*pow(m0,46.0/9)*pow(m1,26.0/9)*pow(m2,28.0/9) + 202*pow(m0,49.0/9)*pow(m1,20.0/9)*pow(m2,31.0/9) + \
    98*pow(m0,29.0/18)*pow(m1,80.0/9)*pow(m2,11.0/18) + 284*pow(m0,35.0/18)*pow(m1,74.0/9)*pow(m2,17.0/18) + \
    690*pow(m0,41.0/18)*pow(m1,68.0/9)*pow(m2,23.0/18) + 330*pow(m0,47.0/18)*pow(m1,62.0/9)*pow(m2,29.0/18) - \
    50*pow(m0,59.0/18)*pow(m1,50.0/9)*pow(m2,41.0/18))) / (144*pow(m0,1.0/18)*pow(m1,25.0/36)*pow(m2,8.0/9)*(2*pow(m1,20.0/9) + \
    pow(m0,29.0/18)*pow(m2,11.0/18))^(5.0/8)*(pow(m1,10.0/3) + pow(m0,13.0/6)*pow(m2,7.0/6) + \
    2*pow(m0,2.0/9)*pow(m1,26.0/9)*pow(m2,2.0/9))^(1.0/4)*(2*pow(m1,50.0/9) + pow(m0,59.0/18)*pow(m2,41.0/18) \
    + 2*m0*pow(m1,32.0/9)*m2 + 2*pow(m0,1.0/9)*pow(m1,16.0/3)*pow(m2,1.0/9) + \
    pow(m0,4.0/9)*pow(m1,14.0/3)*pow(m2,4.0/9))^(9.0/8)) ;

    return diff;
}

