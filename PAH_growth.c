#include "udf.h"
#include "mom.h"

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