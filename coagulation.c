#include "udf.h"
#include "mom.h"

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