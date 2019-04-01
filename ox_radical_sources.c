#include "udf.h"
#include "mom.h"

/* Surface oxidation model
   Units mol/cm^3, s, K  */
#define S1Af 6.09e7
#define S1nf 1.85
#define S1Eaf 7448.0
#define S1Ar 1.07e4
#define S1nr 2.65
#define S1Ear 2796.0
#define S2A 3.26e13
#define S2n 0.17
#define S3A 1.29e14
#define S3Ea 1812.0
#define S4A 6.42e10
#define S4n 4.0
#define S4Ea 29183.0
#define S5A 3.89e3
#define S5n 2.683
#define S5Ea 369.0
#define S6A 1e14
#define S7A 2.4e12
#define S7Ea 2328.0
#define S8A 3.54e11
#define S8n 0.505
#define S8Ea 306.0
#define S10Af 9.8e7
#define S10nf 2.65
#define S10Eaf 8039.0
#define S10Ar 1.6e4
#define S10nr 2.63
#define S10Ear 2137.0
#define S11A 4.86e13
#define S11n 0.13
#define S12A 3.89e3
#define S12n 2.683
#define S12Ea 369.0
#define S13A 9.4e13
#define S13Ea 1771.0
#define S14A 4.34e14
#define S14Ea 984.0
#define S15A 1.47e14
#define S15Ea 632.0
#define S16A 2.0e14
#define S16Ea 2670.0

#define alpha_o 0.35

real h_calc(real temp, real tempmid, real ah1, real ah2, real ah3, real ah4, real ah5, real ah6, 
	real ah7, real al1, real al2, real al3, real al4, real al5, real al6, real al7);

real s_calc(real temp, real tempmid, real ah1, real ah2, real ah3, real ah4, real ah5, real ah6, 
	real ah7, real al1, real al2, real al3, real al4, real al5, real al6, real al7);

real gibbs_calc(int rxn, real temp);

DEFINE_SOURCE(CAH,c,t,dS,eqn)
{
	real cellTemp = C_T(c,t);
	real cellRho = C_R(c,t) * 1e-3; /* in g/cc */

	real CAHMc = UDSI(c,t,3) * C_R(c,t) * 1e-6;
	real CAradMc = UDSI(c,t,4) * C_R(c,t) * 1e-6;
	real CAOMc = UDSI(c,t,5) * C_R(c,t) * 1e-6;
	real CR5Mc = UDSI(c,t,6) * C_R(c,t) * 1e-6;
	real CZHMc = UDSI(c,t,7) * C_R(c,t) * 1e-6;
	real CZradMc = UDSI(c,t,8) * C_R(c,t) * 1e-6;

	Material *mat1 = THREAD_MATERIAL(t);

    int idOH = mixture_specie_index(mat1, "oh");
    real OHMf = Pdf_Yi(c,t,idOH);
    real OHMc = OHMf*cellRho/OHMW;

    int idO = mixture_specie_index(mat1, "o");
    real OMf = Pdf_Yi(c,t,idO);
    real OMc = OMf*cellRho/OMW;

    int idO2 = mixture_specie_index(mat1, "o2");
    real O2Mf = Pdf_Yi(c,t,idO2);
    real O2Mc = O2Mf*cellRho/O2MW;

    int idH = mixture_specie_index(mat1, "h");
    real HMf = Pdf_Yi(c,t,idH);
    real HMc = HMf*cellRho/HMW;

    int idH2 = mixture_specie_index(mat1, "h2");
    real H2Mf = Pdf_Yi(c,t,idH2);
    real H2Mc = H2Mf*cellRho/H2MW;

    real M0 = C_UDSI(c,t,0) * C_R(c,t) / 1e6;   /* Conv				*/
    real M1 = C_UDSI(c,t,1) * C_R(c,t) / 1e6;	/*     to 			*/
    real M2 = C_UDSI(c,t,2) * C_R(c,t) / 1e6;	/*		 /cm^3  	*/
    real Mtwothirds = fractionalMoments(twoThird, M0, M1, M2) * M0;

    real A = CAHMc + CAradMc +CAOMc + CR5Mc;
	real Z = CZHMc + CZradMc;

	real surfArea = pi * pow(6 * mC * 1e3/(pi * rhoSoot * 1e-3), twoThird) * Mtwothirds;
	real chi_nominal = siteDensity * 1e-4 / (avogad * 1e-3) * surfArea;

	real alpha = (A + Z)/chi_nominal;

	real M0_dot = C_UDMI(c,t,50) / 1e6;
	real M1_dot = C_UDMI(c,t,51) / 1e6;
	real M2_dot = C_UDMI(c,t,52) / 1e6;
	real Mtwothirds_dot = 2.0/9.0 * pow(M0, -7.0/9.0) * pow(M1, 8.0/9.0) * pow(M2, -1.0/9.0) * M0_dot
		+ 8.0/9.0 * pow(M0, 2.0/9.0) * pow(M1, -1.0/9.0) * pow(M2, -1.0/9.0) * M1_dot
		- 1.0/9.0 * pow(M0, 2.0/9.0) * pow(M1, 8.0/9.0) * pow(M2, -10.0/9.0) * M2_dot;

	real k1f = S1Af * pow(cellTemp, S1nf) * exp(-S1Eaf/cellTemp) * HMc;
	real k1r = S1Ar * pow(cellTemp, S1nr) * exp(-S1Ear/cellTemp) * H2Mc;
	real k2f = S2A * pow(cellTemp, S2n) * HMc;
	real Keq2 = exp(-gibbs_calc(2, cellTemp)/(R * cellTemp));
	real k2r = k2f/Keq2;
	real k4 = S4A * pow(cellTemp, S4n) * exp(-S4Ea/cellTemp);
	real k5f = S5A * pow(cellTemp, S5n) * exp(-S5Ea/cellTemp) * OHMc;
	real Keq5 = exp(-gibbs_calc(5, cellTemp)/(R * cellTemp));
	real k5r = k5f/Keq5;
	real k6 = S6A * OHMc;
	real k7f = S7A * exp(-S7A/cellTemp) * OMc;
	real Keq7 = exp(-gibbs_calc(7, cellTemp)/(R * cellTemp));
	real k7r = k7f/Keq7;
	real k8 = S8A * pow(cellTemp, S8n) * exp(-S8Ea/cellTemp) * OMc;

	real source_reac = -k1f * CAHMc + k1r * CAradMc + k2f * CAradMc - k2r * CAHMc - k5f * CAHMc + k5r * CAradMc 
		- k7f * CAHMc + k7r * CAOMc + 2 * (k8 * CR5Mc + 1/2 * (k4 * CAOMc + k6 * CAradMc + k8 * CR5Mc)) 
		- 2 * pow(alpha, 2) *  chi_nominal * Mtwothirds_dot / Mtwothirds;

	real nuc_term = C_UDMI(c,t,0) * 1e-6 * alpha_o * (2 * NPyr);

	real source = source_reac + nuc_term;

	dS[eqn] = -k1f - k2r - k5f - k7f ;

	return source;

}

DEFINE_SOURCE(CArad,c,t,dS,eqn)
{
	real cellTemp = C_T(c,t);
	real cellRho = C_R(c,t) * 1e-3; /* in g/cc */

	real CAHMc = UDSI(c,t,3) * C_R(c,t) * 1e-6;
	real CAradMc = UDSI(c,t,4) * C_R(c,t) * 1e-6;
	real CAOMc = UDSI(c,t,5) * C_R(c,t) * 1e-6;
	real CR5Mc = UDSI(c,t,6) * C_R(c,t) * 1e-6;
	real CZHMc = UDSI(c,t,7) * C_R(c,t) * 1e-6;
	real CZradMc = UDSI(c,t,8) * C_R(c,t) * 1e-6;

	Material *mat1 = THREAD_MATERIAL(t);

    int idOH = mixture_specie_index(mat1, "oh");
    real OHMf = Pdf_Yi(c,t,idOH);
    real OHMc = OHMf*cellRho/OHMW;

    int idO = mixture_specie_index(mat1, "o");
    real OMf = Pdf_Yi(c,t,idO);
    real OMc = OMf*cellRho/OMW;

    int idO2 = mixture_specie_index(mat1, "o2");
    real O2Mf = Pdf_Yi(c,t,idO2);
    real O2Mc = O2Mf*cellRho/O2MW;

    int idH = mixture_specie_index(mat1, "h");
    real HMf = Pdf_Yi(c,t,idH);
    real HMc = HMf*cellRho/HMW;

    int idH2 = mixture_specie_index(mat1, "h2");
    real H2Mf = Pdf_Yi(c,t,idH2);
    real H2Mc = H2Mf*cellRho/H2MW;

    real M0 = C_UDSI(c,t,0) * C_R(c,t) / 1e6;   /* Conv				*/
    real M1 = C_UDSI(c,t,1) * C_R(c,t) / 1e6;	/*     to 			*/
    real M2 = C_UDSI(c,t,2) * C_R(c,t) / 1e6;	/*		 /cm^3  	*/
    real Mtwothirds = fractionalMoments(twoThird, M0, M1, M2) * M0;

	real k1f = S1Af * pow(cellTemp, S1nf) * exp(-S1Eaf/cellTemp) * HMc;
	real k1r = S1Ar * pow(cellTemp, S1nr) * exp(-S1Ear/cellTemp) * H2Mc;
	real k2f = S2A * pow(cellTemp, S2n) * HMc;
	real Keq2 = exp(-gibbs_calc(2, cellTemp)/(R * cellTemp));
	real k2r = k2f/Keq2;
	real k3f = S3A * exp(-S3Ea/cellTemp) * O2Mc;
	real Keq3 = exp(-gibbs_calc(3, cellTemp)/(R * cellTemp));
	real k3r = k3f/Keq3;
	real k5f = S5A * pow(cellTemp, S5n) * exp(-S5Ea/cellTemp) * OHMc;
	real Keq5 = exp(-gibbs_calc(5, cellTemp)/(R * cellTemp));
	real k5r = k5f/Keq5;
	real k6 = S6A * OHMc;

	real source_reac = k1f * CAHMc - k1r * CAradMc - k2f * CAradMc + k2r * CAHMc - k3f * CAradMc 
		+ k3r * CAOMc + k5f * CAHMc - k5r * CAradMc - k6 * CAradMc;

	dS[eqn] = - k1r - k2f - k3f - k5r - k6;

	return source_reac;
}

DEFINE_SOURCE(CAO,c,t,dS,eqn)
{
	real cellTemp = C_T(c,t);
	real cellRho = C_R(c,t) * 1e-3; /* in g/cc */

	real CAHMc = UDSI(c,t,3) * C_R(c,t) * 1e-6;
	real CAradMc = UDSI(c,t,4) * C_R(c,t) * 1e-6;
	real CAOMc = UDSI(c,t,5) * C_R(c,t) * 1e-6;
	real CR5Mc = UDSI(c,t,6) * C_R(c,t) * 1e-6;
	real CZHMc = UDSI(c,t,7) * C_R(c,t) * 1e-6;
	real CZradMc = UDSI(c,t,8) * C_R(c,t) * 1e-6;

	Material *mat1 = THREAD_MATERIAL(t);

    int idO = mixture_specie_index(mat1, "o");
    real OMf = Pdf_Yi(c,t,idO);
    real OMc = OMf*cellRho/OMW;

    int idH = mixture_specie_index(mat1, "h");
    real HMf = Pdf_Yi(c,t,idH);
    real HMc = HMf*cellRho/HMW;

    real M0 = C_UDSI(c,t,0) * C_R(c,t) / 1e6;   /* Conv				*/
    real M1 = C_UDSI(c,t,1) * C_R(c,t) / 1e6;	/*     to 			*/
    real M2 = C_UDSI(c,t,2) * C_R(c,t) / 1e6;	/*		 /cm^3  	*/
    real Mtwothirds = fractionalMoments(twoThird, M0, M1, M2) * M0;

	real k3f = S3A * exp(-S3Ea/cellTemp) * HMc;
	real Keq3 = exp(-gibbs_calc(3, cellTemp)/(R * cellTemp));
	real k3r = k3f/Keq3;
	real k4 = S4A * pow(cellTemp, S4n) * exp(-S4Ea/cellTemp);
	real k7f = S7A * exp(-S7A/cellTemp) * OMc;
	real Keq7 = exp(-gibbs_calc(7, cellTemp)/(R*cellTemp));
	real k7r = k7f/Keq7;

	real source_reac = k3f * CAradMc - k3r * CAOMc - k4 * CAOMc + k7f * CAHMc - k7r * CAOMc;

	dS[eqn] = -k3r - k4 - k7r;

	return source_reac;
}

DEFINE_SOURCE(CR5,c,t,dS,eqn)
{
	real cellTemp = C_T(c,t);
	real cellRho = C_R(c,t) * 1e-6; /* in g/cc */

	real CAHMc = UDSI(c,t,3) * C_R(c,t) * 1e-6;
	real CAradMc = UDSI(c,t,4) * C_R(c,t) * 1e-6;
	real CAOMc = UDSI(c,t,5) * C_R(c,t) * 1e-6;
	real CR5Mc = UDSI(c,t,6) * C_R(c,t) * 1e-6;
	real CZHMc = UDSI(c,t,7) * C_R(c,t) * 1e-6;
	real CZradMc = UDSI(c,t,8) * C_R(c,t) * 1e-6;

	Material *mat1 = THREAD_MATERIAL(t);

    int idOH = mixture_specie_index(mat1, "oh");
    real OHMf = Pdf_Yi(c,t,idOH);
    real OHMc = OHMf*cellRho/OHMW;

    int idO = mixture_specie_index(mat1, "o");
    real OMf = Pdf_Yi(c,t,idO);
    real OMc = OMf*cellRho/OMW;

    real M0 = C_UDSI(c,t,0) * C_R(c,t) / 1e6;   /* Conv				*/
    real M1 = C_UDSI(c,t,1) * C_R(c,t) / 1e6;	/*     to 			*/
    real M2 = C_UDSI(c,t,2) * C_R(c,t) / 1e6;	/*		 /cm^3  	*/

	real k4 = S4A * pow(cellTemp, S4n) * exp(-S4Ea/cellTemp);
	real k6 = S6A * OHMc;
	real k8 = S8A * pow(cellTemp, S8n) * exp(-S8Ea/cellTemp) * OMc;

	real source_reac = k4 * CAOMc + k6 * CAradMc - k8 * CR5Mc -  1/2 * (k4 * CAOMc + k6 * CAradMc + k8 * CR5Mc);

	dS[eqn] = - k8 - 1/2 * k8;

	return source_reac;
}

DEFINE_SOURCE(CZH,c,t,dS,eqn)
{
	real cellTemp = C_T(c,t);
	real cellRho = C_R(c,t) * 1e-6; /* in g/cc */

	real CAHMc = UDSI(c,t,3) * C_R(c,t) * 1e-6;
	real CAradMc = UDSI(c,t,4) * C_R(c,t) * 1e-6;
	real CAOMc = UDSI(c,t,5) * C_R(c,t) * 1e-6;
	real CR5Mc = UDSI(c,t,6) * C_R(c,t) * 1e-6;
	real CZHMc = UDSI(c,t,7) * C_R(c,t) * 1e-6;
	real CZradMc = UDSI(c,t,8) * C_R(c,t) * 1e-6;

	Material *mat1 = THREAD_MATERIAL(t);

    int idH = mixture_specie_index(mat1, "h");
    real HMf = Pdf_Yi(c,t,idH);
    real HMc = HMf*cellRho/HMW;

    int idH2 = mixture_specie_index(mat1, "h2");
    real H2Mf = Pdf_Yi(c,t,idH2);
    real H2Mc = H2Mf*cellRho/H2MW;

    int idO = mixture_specie_index(mat1, "o");
    real OMf = Pdf_Yi(c,t,idO);
    real OMc = OMf*cellRho/OMW;

    real M0 = C_UDSI(c,t,0) * C_R(c,t) / 1e6;   /* Conv				*/
    real M1 = C_UDSI(c,t,1) * C_R(c,t) / 1e6;	/*     to 			*/
    real M2 = C_UDSI(c,t,2) * C_R(c,t) / 1e6;	/*		 /cm^3  	*/
    real Mtwothirds = fractionalMoments(twoThird, M0, M1, M2) * M0; 

    real M0_dot = C_UDMI(c,t,50) / 1e6;
	real M1_dot = C_UDMI(c,t,51) / 1e6;
	real M2_dot = C_UDMI(c,t,52) / 1e6;

    real Mtwothirds_dot = 2.0/9.0 * pow(M0, -7.0/9.0) * pow(M1, 8.0/9.0) * pow(M2, -1.0/9.0) * M0_dot
		+ 8.0/9.0 * pow(M0, 2.0/9.0) * pow(M1, -1.0/9.0) * pow(M2, -1.0/9.0) * M1_dot
		- 1.0/9.0 * pow(M0, 2.0/9.0) * pow(M1, 8.0/9.0) * pow(M2, -10.0/9.0) * M2_dot;

	real A = CAHMc + CAradMc +CAOMc + CR5Mc;
	real Z = CZHMc + CZradMc;

	real surfArea = pi * pow(6 * mC * 1e3/(pi * rhoSoot * 1e-3), twoThird) * Mtwothirds;
	real chi_nominal = siteDensity * 1e-4 / (avogad * 1e-3) * surfArea;

	real alpha = (A + Z)/chi_nominal;

	real k8 = S8A * pow(cellTemp, S8n) * exp(-S8Ea/cellTemp) * OMc;
	real k10f = S10Af * pow(cellTemp, S10nf) * exp(-S10Eaf/cellTemp) * HMc;
	real k10r = S10Ar * pow(cellTemp, S10nr) * exp(-S10Ear/cellTemp) * H2Mc;
	real k11f = S11A * pow(cellTemp, S11n) * HMc;
	real Keq11 = exp(-gibbs_calc(11, cellTemp)/(R*cellTemp));
	real k11r = k11f/Keq11;

	real source_reac = 2 * pow(alpha, 2) *  chi_nominal * Mtwothirds_dot / Mtwothirds
		- k10f * CZHMc + k10r * CZradMc + k11f * CZradMc - k11r * CZHMc;

	real nuc_term = C_UDMI(c,t,0) * 1e-6 * (1 - 2 * alpha_o) * (2 * NPyr); 

	real source = source_reac + nuc_term;

	dS[eqn] = -k10f - k11r;

	return source;
}

DEFINE_SOURCE(CZrad,c,t,dS,eqn)
{
	real cellTemp = C_T(c,t);
	real cellRho = C_R(c,t) * 1e-6; /* in g/cc */

	real CAHMc = UDSI(c,t,3) * C_R(c,t) * 1e-6;
	real CAradMc = UDSI(c,t,4) * C_R(c,t) * 1e-6;
	real CAOMc = UDSI(c,t,5) * C_R(c,t) * 1e-6;
	real CR5Mc = UDSI(c,t,6) * C_R(c,t) * 1e-6;
	real CZHMc = UDSI(c,t,7) * C_R(c,t) * 1e-6;
	real CZradMc = UDSI(c,t,8) * C_R(c,t) * 1e-6;

	Material *mat1 = THREAD_MATERIAL(t);

    int idH = mixture_specie_index(mat1, "h");
    real HMf = Pdf_Yi(c,t,idH);
    real HMc = HMf*cellRho/HMW;

    int idH2 = mixture_specie_index(mat1, "h2");
    real H2Mf = Pdf_Yi(c,t,idH2);
    real H2Mc = H2Mf*cellRho/H2MW;

    real M0 = C_UDSI(c,t,0) * C_R(c,t) / 1e6;   /* Conv				*/
    real M1 = C_UDSI(c,t,1) * C_R(c,t) / 1e6;	/*     to 			*/
    real M2 = C_UDSI(c,t,2) * C_R(c,t) / 1e6;	/*		 /cm^3  	*/
    real Mtwothirds = fractionalMoments(twoThird, M0, M1, M2) * M0; 

	real k10f = S10Af * pow(cellTemp, S10nf) * exp(-S10Eaf/cellTemp) * HMc;
	real k10r = S10Ar * pow(cellTemp, S10nr) * exp(-S10Ear/cellTemp) * H2Mc;
	real k11f = S11A * pow(cellTemp, S11n) * HMc;
	real Keq11 = exp(-gibbs_calc(11, cellTemp)/(R*cellTemp));
	real k11r = k11f/Keq11;

	real source_reac = k10f * CZHMc - k10r * CZradMc - k11f * CZradMc + k11r * CZHMc;

	dS[eqn] = -k10r - k11f;

	return source_reac;
}

DEFINE_EXECUTE_AT_END(add_source_terms)
{
	Domain *d;
	Thread *t;

	cell_t c;

	thread_loop_c(t,d)
	{
		begin_c_loop(c,t)
		{
			C_UDMI(c,t,50) = C_UDMI(c,t,0) + C_UDMI(c,t,15) + C_UDMI(c,t,43);
			/*					nuc 	+		coag		+	thermoph 	*/			
			C_UDMI(c,t,51) = C_UDMI(c,t,1) + C_UDMI(c,t,3) + C_UDMI(c,t,7) + C_UDMI(c,t,44);
			/*					nuc 	+		c2h2		+	oxid 		+	thermoph 	*/			
			C_UDMI(c,t,52) = C_UDMI(c,t,2) + C_UDMI(c,t,4) + C_UDMI(c,t,10) + C_UDMI(c,t,18) + C_UDMI(c,t,45);
			/*					nuc 	+		c2h2		+	oxid 		+	coag		+	thermoph 	*/	

			real CAHMc = UDSI(c,t,3) * C_R(c,t) * 1e-6;
			real CAradMc = UDSI(c,t,4) * C_R(c,t) * 1e-6;
			real CAOMc = UDSI(c,t,5) * C_R(c,t) * 1e-6;
			real CR5Mc = UDSI(c,t,6) * C_R(c,t) * 1e-6;
			real CZHMc = UDSI(c,t,7) * C_R(c,t) * 1e-6;
			real CZradMc = UDSI(c,t,8) * C_R(c,t) * 1e-6;

			real M0 = C_UDSI(c,t,0) * C_R(c,t) / 1e6;   /* Conv				*/
    		real M1 = C_UDSI(c,t,1) * C_R(c,t) / 1e6;	/*     to 			*/
    		real M2 = C_UDSI(c,t,2) * C_R(c,t) / 1e6;	/*		 /cm^3  	*/
    		real Mtwothirds = fractionalMoments(twoThird, M0, M1, M2) * M0;

    		real A = CAHMc + CAradMc +CAOMc + CR5Mc;
			real Z = CZHMc + CZradMc;

			real surfArea = pi * pow(6 * mC * 1e3/(pi * rhoSoot * 1e-3), twoThird) * Mtwothirds;
			real chi_nominal = siteDensity * 1e-4 / (avogad * 1e-3) * surfArea;

			real alpha = (A + Z)/chi_nominal;
			C_UDMI(c,t,53) = alpha;
		}
	}
}

real gibbs_calc(int rxn, real temp)
/* returns gibbs energy change in J/mol */
{
	real gibbs = 0.0;

	switch(rxn) {
		case 2: /* CA* + H <=> CAH  
				    a  + b <=>  c    */
		{
			real ah1 = 3.29950506E+00, real ah2 = 6.30133365E-02, real ah3 = -3.79760083E-05, real ah4 = 1.08180756E-08, real ah5 = -1.18007697E-12;
			real ah6 = 4.47516414E+04, real ah7 = 5.41223132E+00, real al1 = -8.00768796E+00, real al2 = 1.03041289E-01, real al3 = -8.38190998E-05;
			real al4 = 2.76491726E-08, real al5 = -8.88842209E-13, real al6 = 4.70598674E+04, real al7 = 6.07299723E+01, real al8 = 2.10265738E+04;
			real a_tempmid = 1000;
			real aH = h_calc(temp, a_tempmid, ah1, ah2, ah3, ah4, ah5, ah6, ah7, al1, al2, al3, al4, al5, al6, al7);
			real aS = s_calc(temp, a_tempmid, ah1, ah2, ah3, ah4, ah5, ah6, ah7, al1, al2, al3, al4, al5, al6, al7);

			real bh1 = 0.25000000E+01, real bh2 = 0.00000000E+00, real bh3 = 0.00000000E+00, real bh4 = 0.00000000E+00, real bh5 = 0.00000000E+00;
	 		real bh6 = 0.25473660E+05, real bh7 = -0.44668285E+00, real bl1 = 0.25000000E+01, real bl2 = 0.00000000E+00, real bl3 = 0.00000000E+00;
	 		real bl4 = 0.00000000E+00, real bl5 = 0.00000000E+00, real bl6 = 0.25473660E+05, real bl7 = -0.44668285E+00, real bl8 = 0.26219035E+05;
	 		real b_tempmid = 1000;
	 		real bH = h_calc(temp, b_tempmid, bh1, bh2, bh3, bh4, bh5, bh6, bh7, bl1, bl2, bl3, bl4, bl5, bl6, bl7);
			real bS = s_calc(temp, b_tempmid, bh1, bh2, bh3, bh4, bh5, bh6, bh7, bl1, bl2, bl3, bl4, bl5, bl6, bl7);

	 		real ch1 = 1.76826275E+00, real ch2 = 6.89143506E-02, real ch3 = -4.14322176E-05, real ch4 = 1.17914309E-08, real ch5 = -1.28597061E-12;
			real ch6 = 1.45412795E+04, real ch7 = 1.06257927E+01, real cl1 = -8.72434585E+00, real cl2 =  1.05376008E-01, real cl3 = -8.01710690E-05;
			real cl4 = 2.18545974E-08, real cl5 =  1.42066606E-12, real cl6 =  1.66588912E+04, real cl7 = 6.19828860E+01, real cl8 = 2.10522087E+04;
			real c_tempmid = 1000;
			real cH = h_calc(temp, c_tempmid, ch1, ch2, ch3, ch4, ch5, ch6, ch7, cl1, cl2, cl3, cl4, cl5, cl6, cl7);
			real cS = s_calc(temp, c_tempmid, ch1, ch2, ch3, ch4, ch5, ch6, ch7, cl1, cl2, cl3, cl4, cl5, cl6, cl7);

			gibbs = cH - (bH + aH) - temp * (cS - (bS + aS));
			break;
		}

		case 3: /* CA* + O2 <=> CAO + O 
				    a  + b  <=>  c  + d   */
		{
			real ah1 = 3.29950506E+00, real ah2 = 6.30133365E-02, real ah3 = -3.79760083E-05, real ah4 = 1.08180756E-08, real ah5 = -1.18007697E-12;
			real ah6 = 4.47516414E+04, real ah7 = 5.41223132E+00, real al1 = -8.00768796E+00, real al2 = 1.03041289E-01, real al3 = -8.38190998E-05;
			real al4 = 2.76491726E-08, real al5 = -8.88842209E-13, real al6 = 4.70598674E+04, real al7 = 6.07299723E+01, real al8 = 2.10265738E+04;
			real a_tempmid = 1000;
			real aH = h_calc(temp, a_tempmid, ah1, ah2, ah3, ah4, ah5, ah6, ah7, al1, al2, al3, al4, al5, al6, al7);
			real aS = s_calc(temp, a_tempmid, ah1, ah2, ah3, ah4, ah5, ah6, ah7, al1, al2, al3, al4, al5, al6, al7);

			real bh1 = 3.66096065E+00, real bh2 = 6.56365811E-04, real bh3 = -1.41149627E-07, real bh4 = 2.05797935E-11, real bh5 =-1.29913436E-15;
			real bh6 = -1.21597718E+03, real bh7 = 3.41536279E+00, real bl1 = 3.78245636E+00, real bl2 = -2.99673416E-03, real bl3 = 9.84730201E-06;
			real bl4 = -9.68129509E-09, real bl5 = 3.24372837E-12, real bl6 = -1.06394356E+03, real bl7 = 3.65767573E+00; 
			real b_tempmid = 1000.0;
			real bH = h_calc(temp, b_tempmid, bh1, bh2, bh3, bh4, bh5, bh6, bh7, bl1, bl2, bl3, bl4, bl5, bl6, bl7);
			real bS = s_calc(temp, b_tempmid, bh1, bh2, bh3, bh4, bh5, bh6, bh7, bl1, bl2, bl3, bl4, bl5, bl6, bl7);

			real ch1 = 2.10591364E+01, real ch2 = 2.82563070E-02, real ch3 = -1.03328686E-05, real ch4 = 1.68867034E-09, real ch5 = -1.01974767E-13;
			real ch6 = 4.09143507E+03, real ch7 = -8.84963398E+01, real cl1 = -1.15176448E+00, real cl2 = 6.11354512E-02, real cl3 = 3.20151083E-05; 
			real cl4 = -9.94285290E-08, real cl5 = 4.79990043E-11, real cl6 = 1.14058756E+04, real cl7 = 3.25584836E+01, real cl8 = 1.38887800E+04;
			real c_tempmid = 1000.0;
			real cH = h_calc(temp, c_tempmid, ch1, ch2, ch3, ch4, ch5, ch6, ch7, cl1, cl2, cl3, cl4, cl5, cl6, cl7);
			real cS = s_calc(temp, c_tempmid, ch1, ch2, ch3, ch4, ch5, ch6, ch7, cl1, cl2, cl3, cl4, cl5, cl6, cl7);

	 		real dh1 = 2.54363697E+00, real dh2 = -2.73162486E-05, real dh3 = -4.19029520E-09, real dh4 = 4.95481845E-12, real dh5 = -4.79553694E-16;
	 		real dh6 = 2.92260120E+04, real dh7 = 4.92229457E+00, real dl1 = 3.16826710E+00, real dl2 = -3.27931884E-03, real dl3 = 6.64306396E-06;
			real dl4 = -6.12806624E-09, real dl5 = 2.11265971E-12, real dl6 = 2.91222592E+04, real dl7 = 2.05193346E+00, real dl8 = 2.99687009E+04;
			real d_tempmid = 1000.0;
			real dH = h_calc(temp, d_tempmid, dh1, dh2, dh3, dh4, dh5, dh6, dh7, dl1, dl2, dl3, dl4, dl5, dl6, dl7);
			real dS = s_calc(temp, d_tempmid, dh1, dh2, dh3, dh4, dh5, dh6, dh7, dl1, dl2, dl3, dl4, dl5, dl6, dl7);

			gibbs = dH + cH - (aH + bH) - temp * (dS + cS - (aS + bS));
			break;
		}

		case 5: /* CAH + OH <=> CA* + H2O
				    a  + b  <=>  c  +  d    */
		{
			real ah1 = 1.76826275E+00, real ah2 = 6.89143506E-02, real ah3 = -4.14322176E-05, real ah4 = 1.17914309E-08, real ah5 = -1.28597061E-12;
			real ah6 = 1.45412795E+04, real ah7 = 1.06257927E+01, real al1 = -8.72434585E+00, real al2 =  1.05376008E-01, real al3 = -8.01710690E-05;
			real al4 = 2.18545974E-08, real al5 =  1.42066606E-12, real al6 =  1.66588912E+04, real al7 = 6.19828860E+01, real al8 = 2.10522087E+04;
			real a_tempmid = 1000.0;
			real aH = h_calc(temp, a_tempmid, ah1, ah2, ah3, ah4, ah5, ah6, ah7, al1, al2, al3, al4, al5, al6, al7);
			real aS = s_calc(temp, a_tempmid, ah1, ah2, ah3, ah4, ah5, ah6, ah7, al1, al2, al3, al4, al5, al6, al7);

	 		real bh1 = 2.83853033E+00, real bh2 = 1.10741289E-03, real bh3 = -2.94000209E-07, real bh4 = 4.20698729E-11, real bh5 = -2.42289890E-15;
	 		real bh6 = 3.69780808E+03, real bh7 = 5.84494652E+00, real bl1 = 3.99198424E+00, real bl2 = -2.40106655E-03, real bl3 = 4.61664033E-06;
			real bl4 = -3.87916306E-09, real bl5 = 1.36319502E-12, real bl6 = 3.36889836E+03, real bl7 = -1.03998477E-01, real bl8 = 4.48615380E+03;
			real b_tempmid = 1000.0;
			real bH = h_calc(temp, b_tempmid, bh1, bh2, bh3, bh4, bh5, bh6, bh7, bl1, bl2, bl3, bl4, bl5, bl6, bl7);
			real bS = s_calc(temp, b_tempmid, bh1, bh2, bh3, bh4, bh5, bh6, bh7, bl1, bl2, bl3, bl4, bl5, bl6, bl7);

			real ch1 = 3.29950506E+00, real ch2 = 6.30133365E-02, real ch3 = -3.79760083E-05, real ch4 = 1.08180756E-08, real ch5 = -1.18007697E-12;
			real ch6 = 4.47516414E+04, real ch7 = 5.41223132E+00, real cl1 = -8.00768796E+00, real cl2 = 1.03041289E-01, real cl3 = -8.38190998E-05;
			real cl4 = 2.76491726E-08, real cl5 = -8.88842209E-13, real cl6 = 4.70598674E+04, real cl7 = 6.07299723E+01, real cl8 = 2.10265738E+04;
			real c_tempmid = 1000.0;
			real cH = h_calc(temp, c_tempmid, ch1, ch2, ch3, ch4, ch5, ch6, ch7, cl1, cl2, cl3, cl4, cl5, cl6, cl7);
			real cS = s_calc(temp, c_tempmid, ch1, ch2, ch3, ch4, ch5, ch6, ch7, cl1, cl2, cl3, cl4, cl5, cl6, cl7);

			real dh1 = 0.26770389E+01, real dh2 = 0.29731816E-02, real dh3 = -0.77376889E-06, real dh4 = 0.94433514E-10, real dh5 = -0.42689991E-14;
			real dh6 = -0.29885894E+05, real dh7 = 0.68825500E+01, real dl1 = 0.41986352E+01, real dl2 = -0.20364017E-02, real dl3 = 0.65203416E-05;
			real dl4 = -0.54879269E-08, real dl5 = 0.17719680E-11, real dl6 = -0.30293726E+05, real dl7 = -0.84900901E+00, real dl8 =-0.29084817E+05;
			real d_tempmid = 1000.0;
			real dH = h_calc(temp, d_tempmid, dh1, dh2, dh3, dh4, dh5, dh6, dh7, dl1, dl2, dl3, dl4, dl5, dl6, dl7);
			real dS = s_calc(temp, d_tempmid, dh1, dh2, dh3, dh4, dh5, dh6, dh7, dl1, dl2, dl3, dl4, dl5, dl6, dl7);

			gibbs = dH + cH - (aH + bH) - temp * (dS + cS - (aS + bS));
			break;
		}

		case 7: /* CAH + O <=> CAO + OH
				    a  + b <=>  c  + d    */
		{
			real ah1 = 1.76826275E+00, real ah2 = 6.89143506E-02, real ah3 = -4.14322176E-05, real ah4 = 1.17914309E-08, real ah5 = -1.28597061E-12;
			real ah6 = 1.45412795E+04, real ah7 = 1.06257927E+01, real al1 = -8.72434585E+00, real al2 =  1.05376008E-01, real al3 = -8.01710690E-05;
			real al4 = 2.18545974E-08, real al5 =  1.42066606E-12, real al6 =  1.66588912E+04, real al7 = 6.19828860E+01, real al8 = 2.10522087E+04;
			real a_tempmid = 1000.0;
			real aH = h_calc(temp, a_tempmid, ah1, ah2, ah3, ah4, ah5, ah6, ah7, al1, al2, al3, al4, al5, al6, al7);
			real aS = s_calc(temp, a_tempmid, ah1, ah2, ah3, ah4, ah5, ah6, ah7, al1, al2, al3, al4, al5, al6, al7);

			real bh1 = 2.54363697E+00, real bh2 = -2.73162486E-05, real bh3 = -4.19029520E-09, real bh4 = 4.95481845E-12, real bh5 = -4.79553694E-16;
	 		real bh6 = 2.92260120E+04, real bh7 = 4.92229457E+00, real bl1 = 3.16826710E+00, real bl2 = -3.27931884E-03, real bl3 = 6.64306396E-06;
			real bl4 = -6.12806624E-09, real bl5 = 2.11265971E-12, real bl6 = 2.91222592E+04, real bl7 = 2.05193346E+00, real bl8 = 2.99687009E+04;
			real b_tempmid = 1000.0;
			real bH = h_calc(temp, b_tempmid, bh1, bh2, bh3, bh4, bh5, bh6, bh7, bl1, bl2, bl3, bl4, bl5, bl6, bl7);
			real bS = s_calc(temp, b_tempmid, bh1, bh2, bh3, bh4, bh5, bh6, bh7, bl1, bl2, bl3, bl4, bl5, bl6, bl7);

			real ch1 = 2.10591364E+01, real ch2 = 2.82563070E-02, real ch3 = -1.03328686E-05, real ch4 = 1.68867034E-09, real ch5 = -1.01974767E-13;
			real ch6 = 4.09143507E+03, real ch7 = -8.84963398E+01, real cl1 = -1.15176448E+00, real cl2 = 6.11354512E-02, real cl3 = 3.20151083E-05; 
			real cl4 = -9.94285290E-08, real cl5 = 4.79990043E-11, real cl6 = 1.14058756E+04, real cl7 = 3.25584836E+01, real cl8 = 1.38887800E+04;
			real c_tempmid = 1000.0;
			real cH = h_calc(temp, c_tempmid, ch1, ch2, ch3, ch4, ch5, ch6, ch7, cl1, cl2, cl3, cl4, cl5, cl6, cl7);
			real cS = s_calc(temp, c_tempmid, ch1, ch2, ch3, ch4, ch5, ch6, ch7, cl1, cl2, cl3, cl4, cl5, cl6, cl7);

			real dh1 = 2.83853033E+00, real dh2 = 1.10741289E-03, real dh3 = -2.94000209E-07, real dh4 = 4.20698729E-11, real dh5 = -2.42289890E-15;
	 		real dh6 = 3.69780808E+03, real dh7 = 5.84494652E+00, real dl1 = 3.99198424E+00, real dl2 = -2.40106655E-03, real dl3 = 4.61664033E-06;
			real dl4 = -3.87916306E-09, real dl5 = 1.36319502E-12, real dl6 = 3.36889836E+03, real dl7 = -1.03998477E-01, real dl8 = 4.48615380E+03;
			real d_tempmid = 1000.0;
			real dH = h_calc(temp, d_tempmid, dh1, dh2, dh3, dh4, dh5, dh6, dh7, dl1, dl2, dl3, dl4, dl5, dl6, dl7);
			real dS = s_calc(temp, d_tempmid, dh1, dh2, dh3, dh4, dh5, dh6, dh7, dl1, dl2, dl3, dl4, dl5, dl6, dl7);

			gibbs = dH + cH - (aH + bH) - temp * (dS + cS - (aS + bS));
			break;
		}

		case 11: /* CZ* + H <=> CZH  
				    a  + b <=>  c    */
		{
			real ah1 = 3.29950506E+00, real ah2 = 6.30133365E-02, real ah3 = -3.79760083E-05, real ah4 = 1.08180756E-08, real ah5 = -1.18007697E-12;
			real ah6 = 4.47516414E+04, real ah7 = 5.41223132E+00, real al1 = -8.00768796E+00, real al2 = 1.03041289E-01, real al3 = -8.38190998E-05;
			real al4 = 2.76491726E-08, real al5 = -8.88842209E-13, real al6 = 4.70598674E+04, real al7 = 6.07299723E+01, real al8 = 2.10265738E+04;
			real a_tempmid = 1000;
			real aH = h_calc(temp, a_tempmid, ah1, ah2, ah3, ah4, ah5, ah6, ah7, al1, al2, al3, al4, al5, al6, al7);
			real aS = s_calc(temp, a_tempmid, ah1, ah2, ah3, ah4, ah5, ah6, ah7, al1, al2, al3, al4, al5, al6, al7);

			real bh1 = 0.25000000E+01, real bh2 = 0.00000000E+00, real bh3 = 0.00000000E+00, real bh4 = 0.00000000E+00, real bh5 = 0.00000000E+00;
	 		real bh6 = 0.25473660E+05, real bh7 = -0.44668285E+00, real bl1 = 0.25000000E+01, real bl2 = 0.00000000E+00, real bl3 = 0.00000000E+00;
	 		real bl4 = 0.00000000E+00, real bl5 = 0.00000000E+00, real bl6 = 0.25473660E+05, real bl7 = -0.44668285E+00, real bl8 = 0.26219035E+05;
	 		real b_tempmid = 1000;
	 		real bH = h_calc(temp, b_tempmid, bh1, bh2, bh3, bh4, bh5, bh6, bh7, bl1, bl2, bl3, bl4, bl5, bl6, bl7);
			real bS = s_calc(temp, b_tempmid, bh1, bh2, bh3, bh4, bh5, bh6, bh7, bl1, bl2, bl3, bl4, bl5, bl6, bl7);

	 		real ch1 = 1.76826275E+00, real ch2 = 6.89143506E-02, real ch3 = -4.14322176E-05, real ch4 = 1.17914309E-08, real ch5 = -1.28597061E-12;
			real ch6 = 1.45412795E+04, real ch7 = 1.06257927E+01, real cl1 = -8.72434585E+00, real cl2 =  1.05376008E-01, real cl3 = -8.01710690E-05;
			real cl4 = 2.18545974E-08, real cl5 =  1.42066606E-12, real cl6 =  1.66588912E+04, real cl7 = 6.19828860E+01, real cl8 = 2.10522087E+04;
			real c_tempmid = 1000;
			real cH = h_calc(temp, c_tempmid, ch1, ch2, ch3, ch4, ch5, ch6, ch7, cl1, cl2, cl3, cl4, cl5, cl6, cl7);
			real cS = s_calc(temp, c_tempmid, ch1, ch2, ch3, ch4, ch5, ch6, ch7, cl1, cl2, cl3, cl4, cl5, cl6, cl7);

			gibbs = cH - (bH + aH) - temp * (cS - (bS + aS));
			break;
		}
	}

	return gibbs;

}

real h_calc(real temp, real tempmid, real ah1, real ah2, real ah3, real ah4, real ah5, real ah6, 
	real ah7, real al1, real al2, real al3, real al4, real al5, real al6, real al7)
{
	real h = 0.0;
	if (temp > tempmid)
	{
		h = (ah1 * temp + ah2/2 * pow(temp, 2) + ah3/3 * pow(temp, 3) + ah4/4 * pow(temp, 4) + ah5/5 * pow(temp, 5) + ah6) * R;
	}

	else
	{
		h = (al1 * temp + al2/2 * pow(temp, 2) + al3/3 * pow(temp, 3) + al4/4 * pow(temp, 4) + al5/5 * pow(temp, 5) + al6) * R;
	}

	return h;
}

real s_calc(real temp, real tempmid, real ah1, real ah2, real ah3, real ah4, real ah5, real ah6, 
	real ah7, real al1, real al2, real al3, real al4, real al5, real al6, real al7)
{
	real s = 0.0;
	if (temp > tempmid)
	{
		s = (ah1 * log(temp) + ah2 * temp + ah3/2 * pow(temp,2) + ah4/3 * pow(temp, 3) + ah5/4 * pow(temp, 4) + ah7) * R;
	}

	else
	{
		s = (al1 * log(temp) + al2 * temp + al3/2 * pow(temp,2) + al4/3 * pow(temp, 3) + al5/4 * pow(temp, 4) + al7) * R;
	}

	return s;
}
