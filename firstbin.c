#include "udf.h"
#include "mom.h"

DEFINE_SOURCE(first_bin_coag,c,t,dS,eqn)
{
	real source;
    real cellRho = C_R(c,t);
    real cellTemp = C_T(c,t);
    real cellPressure = C_P(c,t) + Patm;
    real lamVisc = C_MU_L(c,t);
    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2);
    real cellN1 = cellRho*C_UDSI(c,t,3);
    real m1 = 2 * NPyr;

    real muoneThird = fractionalMoments(oneThird, cellM0, cellM1, cellM2); 
    real meanFreePath = kB * cellTemp/(sqrt(2) * pi * pow(dAir,2) * cellPressure);
    real avgDia = pow(6 * mC/(pi * rhoSoot), oneThird) * muoneThird; 
    real Kn = 2 * meanFreePath/avgDia;

    real Kf = vDW * sqrt(6 * kB * cellTemp/rhoSoot) * pow(3 * mC/(4 * pi * rhoSoot), oneSixth);
    real Kc = 2 * kB * cellTemp/(3 * lamVisc);
    real KcPrime = 2.514 * meanFreePath * pow(pi * rhoSoot/6, oneThird);

    real MminusOneHalf = fractionalMoments(-oneHalf, cellM0, cellM1, cellM2) * cellM0;
    real MoneThird = fractionalMoments(oneThird, cellM0, cellM1, cellM2) * cellM0;
    real MminusOneThird = fractionalMoments(-oneThird, cellM0, cellM1, cellM2) * cellM0;
    real MminusTwoThird = fractionalMoments(-twoThird, cellM0, cellM1, cellM2) * cellM0;
    real MoneHalf = fractionalMoments(oneHalf, cellM0, cellM1, cellM2) * cellM0;
    real MoneSixth = fractionalMoments(oneSixth, cellM0, cellM1, cellM2) * cellM0;
    real MminusOneSixth = fractionalMoments(-oneSixth, cellM0, cellM1, cellM2) * cellM0;
    real MsevenSixth = fractionalMoments(7 * oneSixth, cellM0, cellM1, cellM2) * cellM0;
    real MfiveSixth = fractionalMoments(5 * oneSixth, cellM0, cellM1, cellM2) * cellM0;
    real MthreeHalf = fractionalMoments(3 * oneHalf, cellM0, cellM1, cellM2) * cellM0;
    real MelevenSixth = fractionalMoments(11 * oneSixth, cellM0, cellM1, cellM2) * cellM0;
    real MthirteenSixth = fractionalMoments(13 * oneSixth, cellM0, cellM1, cellM2) * cellM0;

    real f_0 = pow(m1, twoThird) * MminusOneHalf + MoneSixth + 2 * pow(m1, oneThird) * MminusOneSixth;

    real f_1 = pow(m1, 5 * oneThird) * MminusOneHalf + m1 * MoneSixth + 2 * pow(m1, 4 * oneThird) * MminusOneSixth + \
    	pow(m1, twoThird) * MoneHalf + MsevenSixth + 2 * pow(m1, oneThird) * MfiveSixth;

    real f_2 = pow(m1, 8 * oneThird) * MminusOneHalf + 2 * pow(m1, 7 * oneThird) * MminusOneSixth + \
    	pow(m1, 2) * MoneSixth + 2 * pow(m1, 5 * oneThird) * MoneHalf + 4 * pow(m1, 4 * oneThird) * MfiveSixth + \
    	2 * m1 * MsevenSixth + pow(m1, twoThird) * MthreeHalf + 2 * pow(m1, oneThird) * MelevenSixth + MthirteenSixth;

    real f_oneHalf = pow(f_0, 3.0/8.0) * pow(f_1, 3.0/4.0) * pow(f_2, -1.0/8.0);

    real Gc = cellN1 * Kc * (2 * cellM0 + pow(m1, oneThird) * MminusOneThird + pow(m1, -oneThird) * MoneThird + \ 
    	KcPrime * (2 * pow(m1, -oneThird) * cellM0 + MminusOneThird + pow(m1, oneThird) * MminusTwoThird)) * normParameter
    real Gf = cellN1 * pow(m1, -oneHalf) * Kf * normParameter;


    if (Kn <= 0.1)
    {
    	source = Gf;
    }
    else if (Kn >= 8.0)
    {
    	source = Gc;
    }
    else
    {
    	source = Gc * Gf/(Gc + Gf);
    }

    dS[eqn] = 0.0;

    C_UDMI(c,t,52) = source;

    return source;
}

DEFINE_SOURCE(first_bin_oxidation,c,t,dS,eqn)
{
	real source;
    real cellRho = C_R(c,t);
    real cellTemp = C_T(c,t);
    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2);
    real cellN1 = cellRho*C_UDSI(c,t,3);

    int idO2 = mixture_specie_index(mat, "o2");
    int idOH = mixture_specie_index(mat, "oh");

    real O2Mf = Pdf_Yi(c,t,idO2);
    real OHMf = Pdf_Yi(c,t,idOH);

    real O2Mc = O2Mf*cellRho/O2MW;
    real OHMc = OHMf*cellRho/OHMW;
    
    real k5 = A5*exp(-Ea5/(R*cellTemp)) * O2Mc;
    real active_sites = 0.35 * (1 - 2 * 0.3) * siteDensity; 

    real diameter = pow(6 * 32 * mC/(rhoSoot * pi),oneThird);
    real surface_area = pi * pow(diameter, 2);
    
    real delM = 16 * mC;

    source = (-2 * k5 * active_sites * cellN1 * surface_area / 16);

    C_UDMI(c,t,51) = source;

    dS[eqn] = source / cellN1;

    return source;
}

DEFINE_SOURCE(first_bin_growth,c,t,dS,eqn)
{
	real source;
    real cellRho = C_R(c,t);
    real cellTemp = C_T(c,t);
    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2);
    real cellN1 = cellRho*C_UDSI(c,t,3);

    int idC2H2 = mixture_specie_index(mat1, "c2h2");
    real C2H2Mf = Pdf_Yi(c,t,idC2H2);
    real C2H2Mc = C2H2Mf*cellRho/C2H2MW;

    
    real k4 = A4*pow(cellTemp, n4)*exp(-Ea4/(R*cellTemp)) * C2H2Mc;
    real active_sites = 0.35 * (1 - 2 * 0.3) * siteDensity; 

    real diameter = pow(6 * 32 * mC/(rhoSoot * pi),oneThird);
    real surface_area = pi * pow(diameter, 2);
    
    real delM = 2 * mC;

    source = (-2 * k4 * active_sites * cellN1 * surface_area / 2);

    C_UDMI(c,t,50) = source;

    dS[eqn] = source / cellN1;

    return source;

}