#include "udf.h"
#include "mom.h"

DEFINE_SOURCE (m_1_C2H2Source,c,t,dS,eqn)
{
	/*

	First moment HACA surface growth source term. Assumes no turbulence interaction.

	*/

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
	/* 
	
	First moment HACA surface growth source term.Assumes turbulence interaction by considering 
	gaussian PDFs of temperature and nucleating species. PDFs are assumed to be independent.
	
	*/

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

        real f_high = pow(tempUpper, n4) / tempUpper * exp(-Ea4/(R*tempUpper)) * alpha_high * ChiSoot_high * \
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
        integrand = C2H2MfAvg * (f_low + f_high + 2 * f_xi) * tempInt / (2.0 * pdfIntegrand);
        
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

        pdfIntegrand = (P_high + P_low + 2 * P_yi) * C2H2Int / 2.0;
        integrand = pow(cellTempAvg, n4 - 1.0) * exp(-Ea4/(R*cellTempAvg)) * alpha * chiSoot * \
        (f_low + f_high + 2 * f_yi) * C2H2Int / (2.0 * pdfIntegrand);
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

    }

    C_UDMI(c,t,33) = pdfIntegrand;
    

    source = multConstant * integrand;

    real muTwoThirdDiff = fracMomFirstDiffM1(twoThird, cellM0, cellM1, cellM2);

    dS[eqn] = source / muTwoThird * muTwoThirdDiff ;

    /* dS[eqn] = 2 * A4 * cellPressure/(RGas*C2H2MW)*pow(Cs, 2) * pi * muTwoThirdDiff * integrand; */

    C_UDMI(c,t,3) = source;

    return source;
}

DEFINE_SOURCE (temppdf_m_1_C2H2Source,c,t,dS,eqn)
{
    /* 
    
    First moment HACA surface growth source term.Assumes turbulence interaction by considering 
    gaussian PDF of temperature.
    
    */

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
    

    if (tempRatio < 0.02)
    {
        real alpha = HACAalphaCalc(cellTempAvg, cellM0, cellM1);
        real chiSoot = pdfChiSootCalc(cellTempAvg, cellPressure, RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);
        integrand = pow(cellTempAvg, n4) / cellTempAvg * exp(-Ea4/(R*cellTempAvg)) * alpha * chiSoot * C2H2MfAvg;
        pdfIntegrand = 1.0;

    }

    else
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

        real f_high = pow(tempUpper, n4) / tempUpper * exp(-Ea4/(R*tempUpper)) * alpha_high * ChiSoot_high * \
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
        integrand = C2H2MfAvg * (f_low + f_high + 2 * f_xi) * tempInt / (2.0 * pdfIntegrand);
        
    }

    

    source = multConstant * integrand;

    real muTwoThirdDiff = fracMomFirstDiffM1(twoThird, cellM0, cellM1, cellM2);

    dS[eqn] = source / muTwoThird * muTwoThirdDiff ;

    /* dS[eqn] = 2 * A4 * cellPressure/(RGas*C2H2MW)*pow(Cs, 2) * pi * muTwoThirdDiff * integrand; */

    C_UDMI(c,t,3) = source;

    return source;
}

DEFINE_SOURCE (m_2_C2H2Source,c,t,dS,eqn)
{
    /*

	Second moment HACA surface growth source term. Assumes no turbulence interaction.

	*/

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
    
    /* 
	
	Second moment HACA surface growth source term.Assumes turbulence interaction by considering 
	gaussian PDFs of temperature and nucleating species. PDFs are assumed to be independent.
	
	*/

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
        integrand = C2H2MfAvg * (f_low + f_high + 2 * f_xi) * tempInt / (2.0 * pdfIntegrand);
        
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

        pdfIntegrand = (P_high + P_low + 2 * P_yi) * C2H2Int / 2.0;
        integrand = pow(cellTempAvg, n4 - 1.0) * exp(-Ea4/(R*cellTempAvg)) * alpha * chiSoot * \
        (f_low + f_high + 2 * f_yi) * C2H2Int / (2.0 * pdfIntegrand);
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

    }

    C_UDMI(c,t,33) = pdfIntegrand;
    

    source = multConstant * integrand;

    real muTwoThirdDiff = fracMomFirstDiffM2(twoThird, cellM0, cellM1, cellM2);
    real muFiveThirdDiff = fracMomFirstDiffM2(5 * oneThird, cellM0, cellM1, cellM2);

    dS[eqn] = source / (4 * muTwoThird + 4 * muFiveThird) * (4 * muTwoThirdDiff + 4 * muFiveThirdDiff);

    /* dS[eqn] = A4 * cellPressure/(RGas*C2H2MW)*pow(Cs, 2) * pi * (4 * muTwoThirdDiff + 4 * muFiveThirdDiff) * integrand; */

    C_UDMI(c,t,4) = source;

    return source;
}

DEFINE_SOURCE (temppdf_m_2_C2H2Source,c,t,dS,eqn)
{
    
    /* 
    
    Second moment HACA surface growth source term.Assumes turbulence interaction by considering 
    gaussian PDF of temperature.
    
    */

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
    

    if (tempRatio < 0.02)
    {
        real alpha = HACAalphaCalc(cellTempAvg, cellM0, cellM1);
        real chiSoot = pdfChiSootCalc(cellTempAvg, cellPressure, RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);
        integrand = pow(cellTempAvg, n4) / cellTempAvg * exp(-Ea4/(R*cellTempAvg)) * alpha * chiSoot * C2H2MfAvg;
        pdfIntegrand = 1.0;

    }

    else
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

        real f_high = pow(tempUpper, n4) / tempUpper * exp(-Ea4/(R*tempUpper)) * alpha_high * ChiSoot_high * \
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
        integrand = C2H2MfAvg * (f_low + f_high + 2 * f_xi) * tempInt / (2.0 * pdfIntegrand);
        
    }

    real chiSoot = pdfChiSootCalc(cellTempAvg, cellPressure, RGas, HMf, OHMf, H2Mf, H2OMf, C2H2MfAvg, O2Mf);
    C_UDMI(c,t,11) = chiSoot;
    

    source = multConstant * integrand;

    real muTwoThirdDiff = fracMomFirstDiffM2(twoThird, cellM0, cellM1, cellM2);
    real muFiveThirdDiff = fracMomFirstDiffM2(5 * oneThird, cellM0, cellM1, cellM2);

    dS[eqn] = source / (4 * muTwoThird + 4 * muFiveThird) * (4 * muTwoThirdDiff + 4 * muFiveThirdDiff);

    /* dS[eqn] = A4 * cellPressure/(RGas*C2H2MW)*pow(Cs, 2) * pi * (4 * muTwoThirdDiff + 4 * muFiveThirdDiff) * integrand; */

    C_UDMI(c,t,4) = source;

    return source;
}