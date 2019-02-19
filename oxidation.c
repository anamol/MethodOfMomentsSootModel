#include "udf.h"
#include "mom.h"

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

    real alpha_ox = C_UDMI(c,t,32);

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

    sourceO2 = -2*k5*O2Mc*alpha_ox*chiSoot*pi*pow(Cs,2)*MtwoThird;
    /*sourceOH = -gammaOH*OHMc*sqrt(pi*kB*cellTemp/(2*mOH))*pi*pow(Cs,2)*MtwoThird*avogad;*/
    sourceOH = -2*gammaOH*OHMc*sqrt(pi*kB*cellTemp/(2*mOH))*alpha_ox*chiSoot*pi*pow(Cs,2)*MtwoThird;
    source = sourceO2 + sourceOH;

    C_UDMI(c,t,5) = sourceO2;
    C_UDMI(c,t,6) = sourceOH;

    C_UDMI(c,t,7) = source*fudgeFac1;

    dS[eqn] = -2*k5*O2Mc*alpha_ox*chiSoot*pi*pow(Cs,2)*pow(cellM0,twoThird/3)*pow(cellM2,-oneThird/3)*pow(cellM1,-oneThird/3)*8.0/9.0;
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
    real O2multConstant = -2 * A5 * cellPressure * pi * pow(Cs, 2) * muTwoThird * cellM0 / (O2MW * RGas);

    real tempRatio = cellTempRMS/(cellTempAvg);
    real OHMfRatio = OHMfRMS/(OHMfAvg + 1e-8);
    real TARatio = TARMS/TAAvg;
    real O2MfRatio = O2MfRMS / (O2MfAvg + 1e-8);

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
        integrandOH = OHMfAvg * (f_low + f_high + 2 * f_xi) * tempInt / (2.0 * pdfIntegrandOH);
        
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

        pdfIntegrandOH = (P_high + P_low + 2 * P_yi) * OHInt / 2.0;
        integrandOH = 1 / sqrt(cellTempAvg) * (f_low + f_high + 2 * f_yi) * OHInt / (2.0 * pdfIntegrandOH);

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


        real P_low_low = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);
        real f_low_low = OHLower / sqrt(tempLower) * P_low_low;

        real P_high_low = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);
        real f_high_low = OHLower / sqrt(tempUpper) * P_high_low;

        real P_low_high = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);
        real f_low_high = OHUpper / sqrt(tempLower) * P_low_high;

        real P_high_high = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);
        real f_high_high = OHUpper / sqrt(tempUpper) * P_high_high;

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

            f_xi_low = f_xi_low + OHLower / sqrt(temp) * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * \
            gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);

            P_xi_low = P_xi_low + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);

            f_xi_high = f_xi_high + OHUpper / sqrt(temp) * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * \
            gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);

            P_xi_high = P_xi_high + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);

            f_low_yi = f_low_yi + OH / sqrt(tempLower) * gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * \
            gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            P_low_yi = P_low_yi + gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            f_high_yi = f_high_yi + OH / sqrt(tempUpper) * gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * \
            gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            P_high_yi = P_high_yi + gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            for (ctr1 = 1; ctr1 < 13; ctr1++)
            {
                OH_new = OHLower + ctr1 * OHInt;
                f_xi_yi = f_xi_yi + OH / sqrt(temp) * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * \
                gaussianMonoPDFCalc(OH_new, OHMfAvg, OHMfRMS);

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
        /*real alpha = HACAalphaCalc(cellTempAvg, cellM0, cellM1); */
        real alpha = alphaOxCalculator(ta_max, TAAvg);
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

        /*real alpha_low = HACAalphaCalc(tempLower, cellM0, cellM1);
        real alpha_high = HACAalphaCalc(tempUpper, cellM0, cellM1); */

        real alpha_low = alphaOxCalculator(ta_max, TAAvg);
        real alpha_high = alphaOxCalculator(ta_max, TAAvg);

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
            /*real alpha = HACAalphaCalc(temp, cellM0, cellM1); */
            real alpha = alphaOxCalculator(ta_max, TAAvg);
            real ChiSoot = pdfChiSootCalc(temp, cellPressure,RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);

            f_xi = f_xi + 1 / temp * exp(-Ea5/(R*temp)) * alpha * ChiSoot * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
            P_xi = P_xi + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
        }

        pdfIntegrandO2 = (P_low + P_high + 2 * P_xi) * tempInt / 2.0;
        integrandO2 = O2MfAvg * (f_low + f_high + 2 * f_xi) * tempInt / (2.0 * pdfIntegrandO2);
        
    }

    else if (tempRatio < 0.02 && O2MfRatio >= 0.02)
    {
        /*real alpha = HACAalphaCalc(cellTempAvg, cellM0, cellM1); */
        real alpha = alphaOxCalculator(ta_max, TAAvg);
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

        pdfIntegrandO2 = (P_high + P_low + 2 * P_yi) * O2Int / 2.0;
        integrandO2 = 1 / cellTempAvg * exp(-Ea5/(R*cellTempAvg)) * alpha * chiSoot * \
        (f_low + f_high + 2 * f_yi) * O2Int / (2.0 * pdfIntegrandO2);
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

        /*real alpha_low = HACAalphaCalc(tempLower, cellM0, cellM1);
        real alpha_high = HACAalphaCalc(tempUpper, cellM0, cellM1);*/

        real alpha_low = alphaOxCalculator(ta_max, TAAvg);
        real alpha_high = alphaOxCalculator(ta_max, TAAvg);

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
            /*real alpha = HACAalphaCalc(temp, cellM0, cellM1); */
            real alpha = alphaOxCalculator(ta_max, TAAvg);
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

    real OHdiff = sourceOH / muTwoThird * muTwoThirdDiff;
    real O2diff = sourceO2 / muTwoThird * muTwoThirdDiff;

    dS[eqn] = OHdiff + O2diff;
    real source = sourceOH + sourceO2;

    C_UDMI(c,t,5) = sourceOH;
    C_UDMI(c,t,6) = sourceO2;
    C_UDMI(c,t,7) = source;

    return source;
}

DEFINE_SOURCE(pdf_m_1_OxSource_alphaOxPDF,c,t,dS,eqn)
{
	/*real x[ND_ND];
	C_CENTROID(x,c,t);
	Message ("Cell Number: %d, X: %.5g, Y: %.5g \n", c, x[0], x[1]);*/
    real sourceOH = 0.0;
    real integrandOH = 0.0;
    real pdfIntegrandOH = 0.0;

    real sourceO2 = 0.0;
    real integrandO2 = 0.0;
    real pdfIntegrandO2 = 0.0;

    real temp;
    real TA;
    real TA_new;
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

    real ta_max = 6.381;

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

    /*real OHmultConstant =  -gammaOH * cellPressure * pow(Cs, 2) * sqrt(pi * kB/(2 * mOH)) * cellM0 * muTwoThird /(OHMW * RGas);*/

    real OHmultConstant =  -gammaOH * cellPressure * pow(Cs, 2) * sqrt(pi * kB/(2 * mOH)) * cellM0 * muTwoThird /(OHMW * RGas);
    real O2multConstant = -2 * A5 * cellPressure * pi * pow(Cs, 2) * muTwoThird * cellM0 * O2MfAvg / (O2MW * RGas);

    real tempRatio = cellTempRMS/(cellTempAvg);
    real OHMfRatio = OHMfRMS/(OHMfAvg + 1e-8);
    real TARatio = TARMS/(TAAvg + 1e-8);
    real O2MfRatio = O2MfRMS / (O2MfAvg + 1e-8);

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
        integrandOH = OHMfAvg * (f_low + f_high + 2 * f_xi) * tempInt / (2.0 * pdfIntegrandOH);
        
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

        pdfIntegrandOH = (P_high + P_low + 2 * P_yi) * OHInt / 2.0;
        integrandOH = 1 / sqrt(cellTempAvg) * (f_low + f_high + 2 * f_yi) * OHInt / (2.0 * pdfIntegrandOH);

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


        real P_low_low = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);
        real f_low_low = OHLower / sqrt(tempLower) * P_low_low;

        real P_high_low = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);
        real f_high_low = OHLower / sqrt(tempUpper) * P_high_low;

        real P_low_high = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);
        real f_low_high = OHUpper / sqrt(tempLower) * P_low_high;

        real P_high_high = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);
        real f_high_high = OHUpper / sqrt(tempUpper) * P_high_high;

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

            f_xi_low = f_xi_low + OHLower / sqrt(temp) * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * \
            gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);

            P_xi_low = P_xi_low + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);

            f_xi_high = f_xi_high + OHUpper / sqrt(temp) * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * \
            gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);

            P_xi_high = P_xi_high + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);

            f_low_yi = f_low_yi + OH / sqrt(tempLower) * gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * \
            gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            P_low_yi = P_low_yi + gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            f_high_yi = f_high_yi + OH / sqrt(tempUpper) * gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * \
            gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            P_high_yi = P_high_yi + gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            for (ctr1 = 1; ctr1 < 13; ctr1++)
            {
                OH_new = OHLower + ctr1 * OHInt;
                f_xi_yi = f_xi_yi + OH / sqrt(temp) * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * \
                gaussianMonoPDFCalc(OH_new, OHMfAvg, OHMfRMS);

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

    if (tempRatio < 0.02 && TARatio < 0.02)
    {
        real alpha = alphaOxCalculator(ta_max, TAAvg);
        /*real alpha = 0.35;*/
        real chiSoot = pdfChiSootCalc(cellTempAvg, cellPressure, RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
        integrandO2 = 1 / cellTempAvg * exp(-Ea5/(R*cellTempAvg)) * alpha * chiSoot;
        pdfIntegrandO2 = 1.0;

        C_UDMI(c,t,50) = 0.0;
        C_UDMI(c,t,51) = 0.0;
        C_UDMI(c,t,52) = 0.0;
        C_UDMI(c,t,53) = 0.0;
        C_UDMI(c,t,54) = 0.0;
        C_UDMI(c,t,55) = 0.0;
        C_UDMI(c,t,56) = 0.0;
        C_UDMI(c,t,57) = 0.0;
        C_UDMI(c,t,58) = 0.0;
        C_UDMI(c,t,59) = 0.0;
        C_UDMI(c,t,60) = 0.0;

    }


    else if (tempRatio >= 0.02 && TARatio < 0.02)
    {

        real tempLower;

        real tempUpper = cellTempAvg + 3.5 * cellTempRMS;

        if (cellTempAvg - 3.5 * cellTempRMS > 300) { tempLower = cellTempAvg - 3.5 * cellTempRMS; }
        else { tempLower = 300; }

        real tempInt = 0.5 * cellTempRMS;

        /*real alpha_low = HACAalphaCalc(tempLower, cellM0, cellM1);
        real alpha_high = HACAalphaCalc(tempUpper, cellM0, cellM1); */

        real alpha_low = alphaOxCalculator(ta_max, TAAvg);
        real alpha_high = alphaOxCalculator(ta_max, TAAvg);

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
            /*real alpha = HACAalphaCalc(temp, cellM0, cellM1); */
            real alpha = alphaOxCalculator(ta_max, TAAvg);
            real ChiSoot = pdfChiSootCalc(temp, cellPressure,RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);

            f_xi = f_xi + 1 / temp * exp(-Ea5/(R*temp)) * alpha * ChiSoot * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
            P_xi = P_xi + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
        }

        pdfIntegrandO2 = (P_low + P_high + 2 * P_xi) * tempInt / 2.0;
        integrandO2 = (f_low + f_high + 2 * f_xi) * tempInt / (2.0 * pdfIntegrandO2);
        
    }

    else if (tempRatio < 0.02 && TARatio >= 0.02)
    {

        real chiSoot = pdfChiSootCalc(cellTempAvg, cellPressure, RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
        real TALower; real TAUpper;

        if (TAAvg - 3.5 * TARMS > 0) { TALower = TAAvg - 3.5 * TARMS; }
        else { TALower = 0.0; }

        TAUpper = TAAvg + 3.5 * TARMS;

        real TAInt = 0.5 * TARMS;

        int ctr;
        real alpha_low = alphaOxCalculator(ta_max, TALower);
        real alpha_high = alphaOxCalculator(ta_max, TAUpper);

        real f_low = alpha_low * gaussianMonoPDFCalc(TALower, TAAvg, TARMS);
        real P_low = gaussianMonoPDFCalc(TALower, TAAvg, TARMS);

        real f_high = alpha_high * gaussianMonoPDFCalc(TAUpper, TAAvg, TARMS);
        real P_high = gaussianMonoPDFCalc(TAUpper, TAAvg, TARMS);

        real f_yi = 0.0; real P_yi = 0.0;

        for (ctr = 1; ctr < 13; ctr ++)
        {
            TA = TALower + TAInt * ctr;
            real alpha = alphaOxCalculator(ta_max, TA);
            f_yi = f_yi + alpha * gaussianMonoPDFCalc(TA, TAAvg, TARMS);
            P_yi = P_yi + gaussianMonoPDFCalc(TA, TAAvg, TARMS);
        }

        pdfIntegrandO2 = (P_high + P_low + 2 * P_yi) * TAInt / 2.0;
        integrandO2 = 1 / cellTempAvg * exp(-Ea5/(R*cellTempAvg)) * chiSoot * \
        (f_low + f_high + 2 * f_yi) * TAInt / (2.0 * pdfIntegrandO2);
    }

    else 
    {
    	/*real alpha = alphaOxCalculator(ta_max, TAAvg);*/

        real tempUpper; real tempLower;
        real TAUpper; real TALower;

        tempUpper = cellTempAvg + 3.5 * cellTempRMS;
        if (cellTempAvg - 3.5 * cellTempRMS > 300) { tempLower = cellTempAvg - 3.5 * cellTempRMS; }
        else { tempLower = 300; }

        if (TAAvg - 3.5 * TARMS > 0) { TALower = TAAvg - 3.5 * TARMS; }
        else { TALower = 0.0; }

        TAUpper = TAAvg + 3.5 * TARMS;

        real tempInt = 0.5 * cellTempRMS; real TAInt = 0.5 * TARMS;

        real ChiSoot_low = pdfChiSootCalc(tempLower, cellPressure,RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
        real ChiSoot_high = pdfChiSootCalc(tempUpper, cellPressure,RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);

        real f_low_low = O2IntegrandCalc_TApdf(tempLower, cellTempAvg, cellTempRMS, TALower, TAAvg, TARMS, ta_max, ChiSoot_low);
        real P_low_low = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(TALower, TAAvg, TARMS);

        real f_high_low = O2IntegrandCalc_TApdf(tempUpper, cellTempAvg, cellTempRMS, TALower, TAAvg, TARMS, ta_max, ChiSoot_high);
        real P_high_low = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(TALower, TAAvg, TARMS);

        real f_low_high = O2IntegrandCalc_TApdf(tempLower, cellTempAvg, cellTempRMS, TAUpper, TAAvg, TARMS, ta_max, ChiSoot_low);
        real P_low_high = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(TALower, TAAvg, TARMS);

        real f_high_high = O2IntegrandCalc_TApdf(tempUpper, cellTempAvg, cellTempRMS, TAUpper, TAAvg, TARMS, ta_max, ChiSoot_high);
        real P_high_high = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(TALower, TAAvg, TARMS);

        real f_xi_low = 0.0; real P_xi_low = 0.0;
        real f_xi_high = 0.0; real P_xi_high = 0.0;
        real f_low_yi = 0.0; real P_low_yi = 0.0;
        real f_high_yi = 0.0; real P_high_yi = 0.0;
        real f_xi_yi = 0.0; real P_xi_yi = 0.0;

        int ctr;
        int ctr1;

        for (ctr = 1; ctr < 13; ctr++)
        {
            temp = tempLower + ctr * tempInt; TA = TALower + ctr * TAInt;
    
            real ChiSoot = pdfChiSootCalc(temp, cellPressure,RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
            f_xi_low = f_xi_low + O2IntegrandCalc_TApdf(temp, cellTempAvg, cellTempRMS, TALower, TAAvg, TARMS, ta_max, ChiSoot);
            P_xi_low = P_xi_low + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(TALower, TAAvg, TARMS);

            f_xi_high = f_xi_high + O2IntegrandCalc_TApdf(temp, cellTempAvg, cellTempRMS, TAUpper, TAAvg, TARMS, ta_max, ChiSoot);
            P_xi_high = P_xi_high + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(TAUpper, TAAvg, TARMS);

            f_low_yi = f_low_yi + O2IntegrandCalc_TApdf(tempLower, cellTempAvg, cellTempRMS, TA, TAAvg, TARMS, ta_max, ChiSoot_low);
            P_low_yi = P_low_yi + gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(TA, TAAvg, TARMS);

            f_high_yi = f_high_yi + O2IntegrandCalc_TApdf(tempUpper, cellTempAvg, cellTempRMS, TA, TAAvg, TARMS, ta_max, ChiSoot_high);
            P_high_yi = P_high_yi + gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(TA, TAAvg, TARMS);

            for (ctr1 = 1; ctr1 < 13; ctr1++)
            {
                TA_new = TALower + ctr1 * TAInt;
                f_xi_yi = f_xi_yi + O2IntegrandCalc_TApdf(temp, cellTempAvg, cellTempRMS, TA_new, TAAvg, TARMS, ta_max, ChiSoot);
                P_xi_yi = P_xi_yi + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(TA_new, TAAvg, TARMS);
            }
        }
        C_UDMI(c,t,50) = f_xi_yi;
        C_UDMI(c,t,51) = f_high_low;
        C_UDMI(c,t,52) = f_low_high;
        C_UDMI(c,t,53) = f_high_high;
        C_UDMI(c,t,54) = f_low_low;
        C_UDMI(c,t,55) = f_low_yi;
        C_UDMI(c,t,56) = f_high_yi;
        C_UDMI(c,t,57) = f_xi_low;
        C_UDMI(c,t,58) = f_xi_high;
        C_UDMI(c,t,59) = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS); 
        C_UDMI(c,t,60) = gaussianMonoPDFCalc(TALower, TAAvg, TARMS);



        pdfIntegrandO2 = 0.25 * tempInt * TAInt * (P_low_high + P_high_high + P_high_low + P_low_low + \
            2 * (P_low_yi + P_high_yi + P_xi_low + P_xi_high) + 4 * P_xi_yi);

        integrandO2 = 0.25 * tempInt * TAInt * (f_low_high + f_high_high + f_high_low + f_low_low + \
            2 * (f_low_yi + f_high_yi + f_xi_low + f_xi_high) + 4 * f_xi_yi);

        integrandO2 = integrandO2 / pdfIntegrandO2; 

        /*real alpha = alphaOxCalculator(ta_max, TAAvg);
        real chiSoot = pdfChiSootCalc(cellTempAvg, cellPressure, RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
        integrandO2 = 1 / cellTempAvg * exp(-Ea5/(R*cellTempAvg)) * alpha * chiSoot +f_low_low + f_high_low;
        pdfIntegrandO2 = 1.0 + P_low_low + P_high_low; */

    }

    sourceOH = OHmultConstant * integrandOH;
    sourceO2 = O2multConstant * integrandO2;

    real muTwoThirdDiff = fracMomFirstDiffM1(twoThird, cellM0, cellM1, cellM2);

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

    real alpha_ox = C_UDMI(c,t,32);

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

    sourceO2 = k5*O2Mc*alpha_ox*chiSoot*pi*pow(Cs,2)*(4*MtwoThird - 4*MfiveThird);
    /*sourceOH = gammaOH*OHMc*sqrt(pi*kB*cellTemp/(2*mOH))*pi*pow(Cs,2)*(4*MtwoThird - 4*MfiveThird);*/
    sourceOH = gammaOH*OHMc*sqrt(pi*kB*cellTemp/(2*mOH))*alpha_ox*chiSoot*pi*pow(Cs,2)*(4*MtwoThird - 4*MfiveThird);
    source = sourceO2 + sourceOH;

    C_UDMI(c,t,8) = sourceO2;
    C_UDMI(c,t,9) = sourceOH;

    dS[eqn] = k5*O2Mc*alpha_ox*chiSoot*pi*pow(Cs,2)*(-4*pow(cellM0,twoThird/3)*pow(cellM2,-11*oneThird/3)*\
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

    real OHmultConstant =  gammaOH * cellPressure * pow(Cs, 2) * sqrt(pi * kB/(2 * mOH)) * avogad * cellM0 * \
    (4 * muTwoThird - 4 * muFiveThird) /(OHMW * RGas);

    real O2multConstant = A5 * cellPressure * pi * pow(Cs, 2) * (4 * muTwoThird - 4 * muFiveThird) * cellM0 / (O2MW * RGas);

    real tempRatio = cellTempRMS/(cellTempAvg);
    real OHMfRatio = OHMfRMS/(OHMfAvg + 1e-8);
    real TARatio = TARMS/TAAvg;
    real O2MfRatio = O2MfRMS / (O2MfAvg + 1e-8);

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
        integrandOH = OHMfAvg * (f_low + f_high + 2 * f_xi) * tempInt / (2.0 * pdfIntegrandOH);
        
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
        integrandOH = 1 / sqrt(cellTempAvg) * (f_low + f_high + 2 * f_yi) * OHInt / (2.0 * pdfIntegrandOH);

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


        real P_low_low = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);
        real f_low_low = OHLower / sqrt(tempLower) * P_low_low;

        real P_high_low = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);
        real f_high_low = OHLower / sqrt(tempUpper) * P_high_low;

        real P_low_high = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);
        real f_low_high = OHUpper / sqrt(tempLower) * P_low_high;

        real P_high_high = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);
        real f_high_high = OHUpper / sqrt(tempUpper) * P_high_high;

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

            f_xi_low = f_xi_low + OHLower / sqrt(temp) * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * \
            gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);

            P_xi_low = P_xi_low + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);

            f_xi_high = f_xi_high + OHUpper / sqrt(temp) * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * \
            gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);

            P_xi_high = P_xi_high + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);

            f_low_yi = f_low_yi + OH / sqrt(tempLower) * gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * \
            gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            P_low_yi = P_low_yi + gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            f_high_yi = f_high_yi + OH / sqrt(tempUpper) * gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * \
            gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            P_high_yi = P_high_yi + gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            for (ctr1 = 1; ctr1 < 13; ctr1++)
            {
                OH_new = OHLower + ctr1 * OHInt;
                f_xi_yi = f_xi_yi + OH / sqrt(temp) * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * \
                gaussianMonoPDFCalc(OH_new, OHMfAvg, OHMfRMS);
                
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
        /*real alpha = HACAalphaCalc(cellTempAvg, cellM0, cellM1); */
        real alpha = alphaOxCalculator(ta_max, TAAvg);
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

        /*real alpha_low = HACAalphaCalc(tempLower, cellM0, cellM1);
        real alpha_high = HACAalphaCalc(tempUpper, cellM0, cellM1); */

        real alpha_low = alphaOxCalculator(ta_max, TAAvg);
        real alpha_high = alphaOxCalculator(ta_max, TAAvg);

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
            /*real alpha = HACAalphaCalc(temp, cellM0, cellM1); */
            real alpha = alphaOxCalculator(ta_max, TAAvg);  
            real ChiSoot = pdfChiSootCalc(temp, cellPressure,RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);

            f_xi = f_xi + pow(temp, - 1.0) * exp(-Ea5/(R*temp)) * alpha * ChiSoot * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
            P_xi = P_xi + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
        }

        pdfIntegrandO2 = (P_low + P_high + 2 * P_xi) * tempInt / 2.0;
        integrandO2 = O2MfAvg * (f_low + f_high + 2 * f_xi) * tempInt / (2.0 * pdfIntegrandO2);
        
    }

    else if (tempRatio < 0.02 && O2MfRatio >= 0.02)
    {
        /*real alpha = HACAalphaCalc(cellTempAvg, cellM0, cellM1); */
        real alpha = alphaOxCalculator(ta_max, TAAvg);      
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

        pdfIntegrandO2 = (P_high + P_low + 2 * P_yi) * O2Int / 2.0;
        integrandO2 = 1 / cellTempAvg * exp(-Ea5/(R*cellTempAvg)) * alpha * chiSoot * \
        (f_low + f_high + 2 * f_yi) * O2Int / (2.0 * pdfIntegrandO2);
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

        /*real alpha_low = HACAalphaCalc(tempLower, cellM0, cellM1);
        real alpha_high = HACAalphaCalc(tempUpper, cellM0, cellM1); */

        real alpha_low = alphaOxCalculator(ta_max, TAAvg);
        real alpha_high = alphaOxCalculator(ta_max, TAAvg);

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
            /*real alpha = HACAalphaCalc(temp, cellM0, cellM1); */
            real alpha = alphaOxCalculator(ta_max, TAAvg);        
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

    real OHdiff = sourceOH / (4 * muTwoThird - 4 * muFiveThird) * (4 * muTwoThirdDiff - 4 * muFiveThirdDiff);
    real O2diff = sourceO2 / (4 * muTwoThird - 4 * muFiveThird) * (4 * muTwoThirdDiff - 4 * muFiveThirdDiff);

    dS[eqn] = OHdiff + O2diff;
    real source = sourceOH + sourceO2;

    C_UDMI(c,t,8) = sourceOH;
    C_UDMI(c,t,9) = sourceO2;
    C_UDMI(c,t,10) = source;

    return source;
}

DEFINE_SOURCE(pdf_m_2_OxSource_alphaOxPDF,c,t,dS,eqn)
{
    real sourceOH = 0.0;
    real integrandOH = 0.0;
    real pdfIntegrandOH = 0.0;

    real sourceO2 = 0.0;
    real integrandO2 = 0.0;
    real pdfIntegrandO2 = 0.0;

    real temp;
    real TA;
    real TA_new;
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

    real ta_max = 13.5;

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
    (4 * muTwoThird - 4 * muFiveThird) /(OHMW * RGas);

    real O2multConstant = A5 * cellPressure * pi * pow(Cs, 2) * (4 * muTwoThird - 4 * muFiveThird) * cellM0 * O2MfAvg / (O2MW * RGas);

    real tempRatio = cellTempRMS/(cellTempAvg);
    real OHMfRatio = OHMfRMS/(OHMfAvg + 1e-8);
    real TARatio = TARMS/(TAAvg + 1e-8);
    real O2MfRatio = O2MfRMS / (O2MfAvg + 1e-8);

    C_UDMI(c,t,33) = tempRatio;
    C_UDMI(c,t,34) = OHMfRatio;
    C_UDMI(c,t,35) = TARatio;

    if (tempRatio < 0.02 && OHMfRatio < 0.02)
    {

        integrandOH = OHMfAvg /sqrt(cellTempAvg);
        pdfIntegrandOH = 1.0;

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


        real P_low_low = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);
        real f_low_low = OHLower / sqrt(tempLower) * P_low_low;

        real P_high_low = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);
        real f_high_low = OHLower / sqrt(tempUpper) * P_high_low;

        real P_low_high = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);
        real f_low_high = OHUpper / sqrt(tempLower) * P_low_high;

        real P_high_high = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);
        real f_high_high = OHUpper / sqrt(tempUpper) * P_high_high;

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

            f_xi_low = f_xi_low + OHLower / sqrt(temp) * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * \
            gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);

            P_xi_low = P_xi_low + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHLower, OHMfAvg, OHMfRMS);

            f_xi_high = f_xi_high + OHUpper / sqrt(temp) * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * \
            gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);

            P_xi_high = P_xi_high + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OHUpper, OHMfAvg, OHMfRMS);

            f_low_yi = f_low_yi + OH / sqrt(tempLower) * gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * \
            gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            P_low_yi = P_low_yi + gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            f_high_yi = f_high_yi + OH / sqrt(tempUpper) * gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * \
            gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            P_high_yi = P_high_yi + gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(OH, OHMfAvg, OHMfRMS);

            for (ctr1 = 1; ctr1 < 13; ctr1++)
            {
                OH_new = OHLower + ctr1 * OHInt;
                f_xi_yi = f_xi_yi + OH / sqrt(temp) * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * \
                gaussianMonoPDFCalc(OH_new, OHMfAvg, OHMfRMS);
                
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


    if (tempRatio < 0.02 && TARatio < 0.02)
    {
        real alpha = alphaOxCalculator(ta_max, TAAvg);
        real chiSoot = pdfChiSootCalc(cellTempAvg, cellPressure, RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
        integrandO2 = 1 / cellTempAvg * exp(-Ea5/(R*cellTempAvg)) * alpha * chiSoot;
        pdfIntegrandO2 = 1.0;

    }

    else if (tempRatio >= 0.02 && TARatio < 0.02)
    {

        real tempLower;

        real tempUpper = cellTempAvg + 3.5 * cellTempRMS;

        if (cellTempAvg - 3.5 * cellTempRMS > 300) { tempLower = cellTempAvg - 3.5 * cellTempRMS; }
        else { tempLower = 300; }

        real tempInt = 0.5 * cellTempRMS;

        /*real alpha_low = HACAalphaCalc(tempLower, cellM0, cellM1);
        real alpha_high = HACAalphaCalc(tempUpper, cellM0, cellM1); */

        real alpha_low = alphaOxCalculator(ta_max, TAAvg);
        real alpha_high = alphaOxCalculator(ta_max, TAAvg);

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
            /*real alpha = HACAalphaCalc(temp, cellM0, cellM1); */
            real alpha = alphaOxCalculator(ta_max, TAAvg);
            real ChiSoot = pdfChiSootCalc(temp, cellPressure,RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);

            f_xi = f_xi + 1 / temp * exp(-Ea5/(R*temp)) * alpha * ChiSoot * gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
            P_xi = P_xi + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS);
        }

        pdfIntegrandO2 = (P_low + P_high + 2 * P_xi) * tempInt / 2.0;
        integrandO2 = (f_low + f_high + 2 * f_xi) * tempInt / (2.0 * pdfIntegrandO2);
        
    }

    else if (tempRatio < 0.02 && TARatio >= 0.02)
    {

        real chiSoot = pdfChiSootCalc(cellTempAvg, cellPressure, RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
        real TALower; real TAUpper;

        if (TAAvg - 3.5 * TARMS > 0) { TALower = TAAvg - 3.5 * TARMS; }
        else { TALower = 0.0; }

        TAUpper = TAAvg + 3.5 * TARMS;

        real TAInt = 0.5 * TARMS;

        int ctr;
        real alpha_low = alphaOxCalculator(ta_max, TALower);
        real alpha_high = alphaOxCalculator(ta_max, TAUpper);

        real f_low = alpha_low * gaussianMonoPDFCalc(TALower, TAAvg, TARMS);
        real P_low = gaussianMonoPDFCalc(TALower, TAAvg, TARMS);

        real f_high = alpha_high * gaussianMonoPDFCalc(TAUpper, TAAvg, TARMS);
        real P_high = gaussianMonoPDFCalc(TAUpper, TAAvg, TARMS);

        real f_yi = 0.0; real P_yi = 0.0;

        for (ctr = 1; ctr < 13; ctr ++)
        {
            TA = TALower + TAInt * ctr;
            real alpha = alphaOxCalculator(ta_max, TA);
            f_yi = f_yi + alpha * gaussianMonoPDFCalc(TA, TAAvg, TARMS);
            P_yi = P_yi + gaussianMonoPDFCalc(TA, TAAvg, TARMS);
        }

        pdfIntegrandO2 = (P_high + P_low + 2 * P_yi) * TAInt / 2.0;
        integrandO2 = 1 / cellTempAvg * exp(-Ea5/(R*cellTempAvg)) * chiSoot * \
        (f_low + f_high + 2 * f_yi) * TAInt / (2.0 * pdfIntegrandO2);
    }

    else 
    {
        real tempUpper; real tempLower;
        real TAUpper; real TALower;

        tempUpper = cellTempAvg + 3.5 * cellTempRMS;
        if (cellTempAvg - 3.5 * cellTempRMS > 300) { tempLower = cellTempAvg - 3.5 * cellTempRMS; }
        else { tempLower = 300; }

        if (TAAvg - 3.5 * TARMS > 0) { TALower = TAAvg - 3.5 * TARMS; }
        else { TALower = 0.0; }

        TAUpper = TAAvg + 3.5 * TARMS;

        real tempInt = 0.5 * cellTempRMS; real TAInt = 0.5 * TARMS;


        real ChiSoot_low = pdfChiSootCalc(tempLower, cellPressure,RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
        real ChiSoot_high = pdfChiSootCalc(tempUpper, cellPressure,RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);

        real f_low_low = O2IntegrandCalc_TApdf(tempLower, cellTempAvg, cellTempRMS, TALower, TAAvg, TARMS, ta_max, ChiSoot_low);
        real P_low_low = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(TALower, TAAvg, TARMS);

        real f_high_low = O2IntegrandCalc_TApdf(tempUpper, cellTempAvg, cellTempRMS, TALower, TAAvg, TARMS, ta_max, ChiSoot_high);
        real P_high_low = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(TALower, TAAvg, TARMS);

        real f_low_high = O2IntegrandCalc_TApdf(tempLower, cellTempAvg, cellTempRMS, TAUpper, TAAvg, TARMS, ta_max, ChiSoot_low);
        real P_low_high = gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(TALower, TAAvg, TARMS);

        real f_high_high = O2IntegrandCalc_TApdf(tempUpper, cellTempAvg, cellTempRMS, TAUpper, TAAvg, TARMS, ta_max, ChiSoot_high);
        real P_high_high = gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(TALower, TAAvg, TARMS);

        real f_xi_low = 0.0; real P_xi_low = 0.0;
        real f_xi_high = 0.0; real P_xi_high = 0.0;
        real f_low_yi = 0.0; real P_low_yi = 0.0;
        real f_high_yi = 0.0; real P_high_yi = 0.0;
        real f_xi_yi = 0.0; real P_xi_yi = 0.0;

        int ctr;
        int ctr1;

        for (ctr = 1; ctr < 13; ctr++)
        {
            temp = tempLower + ctr * tempInt; TA = TALower + ctr * TAInt;
    
            real ChiSoot = pdfChiSootCalc(temp, cellPressure,RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
            f_xi_low = f_xi_low + O2IntegrandCalc_TApdf(temp, cellTempAvg, cellTempRMS, TALower, TAAvg, TARMS, ta_max, ChiSoot);
            P_xi_low = P_xi_low + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(TALower, TAAvg, TARMS);

            f_xi_high = f_xi_high + O2IntegrandCalc_TApdf(temp, cellTempAvg, cellTempRMS, TAUpper, TAAvg, TARMS, ta_max, ChiSoot);
            P_xi_high = P_xi_high + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(TAUpper, TAAvg, TARMS);

            f_low_yi = f_low_yi + O2IntegrandCalc_TApdf(tempLower, cellTempAvg, cellTempRMS, TA, TAAvg, TARMS, ta_max, ChiSoot_low);
            P_low_yi = P_low_yi + gaussianMonoPDFCalc(tempLower, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(TA, TAAvg, TARMS);

            f_high_yi = f_high_yi + O2IntegrandCalc_TApdf(tempUpper, cellTempAvg, cellTempRMS, TA, TAAvg, TARMS, ta_max, ChiSoot_high);
            P_high_yi = P_high_yi + gaussianMonoPDFCalc(tempUpper, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(TA, TAAvg, TARMS);

            for (ctr1 = 1; ctr1 < 13; ctr1++)
            {
                TA_new = TALower + ctr1 * TAInt;
                f_xi_yi = f_xi_yi + O2IntegrandCalc_TApdf(temp, cellTempAvg, cellTempRMS, TA_new, TAAvg, TARMS, ta_max, ChiSoot);
                P_xi_yi = P_xi_yi + gaussianMonoPDFCalc(temp, cellTempAvg, cellTempRMS) * gaussianMonoPDFCalc(TA_new, TAAvg, TARMS);
            }
        }
        
        pdfIntegrandO2 = 0.25 * tempInt * TAInt * (P_low_high + P_high_high + P_high_low + P_low_low + \
            2 * (P_low_yi + P_high_yi + P_xi_low + P_xi_high) + 4 * P_xi_yi);

        integrandO2 = 0.25 * tempInt * TAInt * (f_low_high + f_high_high + f_high_low + f_low_low + \
            2 * (f_low_yi + f_high_yi + f_xi_low + f_xi_high) + 4 * f_xi_yi);

        integrandO2 = integrandO2 / pdfIntegrandO2;

        /*real alpha = alphaOxCalculator(ta_max, TAAvg);
        real chiSoot = pdfChiSootCalc(cellTempAvg, cellPressure, RGas, HMf, OHMfAvg, H2Mf, H2OMf, C2H2Mf, O2MfAvg);
        integrandO2 = 1 / cellTempAvg * exp(-Ea5/(R*cellTempAvg)) * alpha * chiSoot +f_low_low + f_high_low;
        pdfIntegrandO2 = 1.0 + P_low_low + P_high_low; */

    }
    

    sourceOH = OHmultConstant * integrandOH;
    sourceO2 = O2multConstant * integrandO2;

    real muTwoThirdDiff = fracMomFirstDiffM2(twoThird, cellM0, cellM1, cellM2);
    real muFiveThirdDiff = fracMomFirstDiffM2(5.0 * oneThird, cellM0, cellM1, cellM2);

    real OHdiff = sourceOH / (4 * muTwoThird - 4 * muFiveThird) * (4 * muTwoThirdDiff - 4 * muFiveThirdDiff);
    real O2diff = sourceO2 / (4 * muTwoThird - 4 * muFiveThird) * (4 * muTwoThirdDiff - 4 * muFiveThirdDiff);

    dS[eqn] = OHdiff + O2diff;
    real source = sourceOH + sourceO2;

    C_UDMI(c,t,8) = sourceOH;
    C_UDMI(c,t,9) = sourceO2;
    C_UDMI(c,t,10) = source;

    return source;
}
