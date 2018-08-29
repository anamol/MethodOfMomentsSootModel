#include "udf.h"
#include "mom.h"

DEFINE_SOURCE (m_0_NucSourcePyr,c,t,dS,eqn)
{
	/* 
	
	Zeroth moment source term for nucleation assuming pyrene as the nucleating species.
	Assumes no turbulence interaction.

	*/

    real sourcePyr = pyrNucSource(c,t); 

    dS[eqn] = 0.0;
    C_UDMI(c,t,0) = sourcePyr;
    return sourcePyr;
}

DEFINE_SOURCE (pdf_m_0_NucSourcePyr,c,t,dS,eqn)
{
    /* 
	
	Zeroth moment source term for nucleation assuming pyrene as the nucleating species.
	Assumes turbulence interaction by considering gaussian PDFs of temperature and 
	nucleating species. PDFs are assumed to be independent.
	
	*/

    real sourcePyr = pdfpyrNucSource(c,t); 

    dS[eqn] = 0.0;
    C_UDMI(c,t,0) = sourcePyr;
    return sourcePyr;
}

DEFINE_SOURCE (m_0_NucSourceBenz,c,t,dS,eqn)
{
	/* 
	
	Zeroth moment source term for nucleation assuming benzene as the nucleating species.
	Assumes no turbulence interaction.

	*/

    real sourceBenz = benzNucSource(c,t); 

    dS[eqn] = 0.0;
    C_UDMI(c,t,0) = sourceBenz;
    return sourceBenz;
}

DEFINE_SOURCE (m_1_NucSourcePyr,c,t,dS,eqn)
{
	/* 
	
	First moment source term for nucleation assuming pyrene as the nucleating species.
	Assumes no turbulence interaction.

	*/

    real sourcePyr = (2*NPyr)*pyrNucSource(c,t); 
    real source = sourcePyr; 


    dS[eqn] = 0.0;
    C_UDMI(c,t,1) = source;
    return source;
}

DEFINE_SOURCE (pdf_m_1_NucSourcePyr,c,t,dS,eqn)
{
	/* 
	
	First moment source term for nucleation assuming pyrene as the nucleating species.
	Assumes turbulence interaction by considering gaussian PDFs of temperature and 
	nucleating species. PDFs are assumed to be independent.
	
	*/

    real sourcePyr = (2*NPyr)*pdfpyrNucSource(c,t); 
    real source = sourcePyr; 


    dS[eqn] = 0.0;
    C_UDMI(c,t,1) = source;
    return source;
}

DEFINE_SOURCE (m_1_NucSourceBenz,c,t,dS,eqn)
{
	/* 
	
	First moment source term for nucleation assuming benzene as the nucleating species.
	Assumes no turbulence interaction.

	*/

    real sourceBenz = (2*NBenzene)*benzNucSource(c,t); 
    real source = sourceBenz;

    dS[eqn] = 0.0;
    C_UDMI(c,t,1) = source;
    return source;
}

DEFINE_SOURCE (m_2_NucSourcePyr,c,t,dS,eqn)
{
	/* 
	
	Second moment source term for nucleation assuming pyrene as the nucleating species.
	Assumes no turbulence interaction.

	*/

    real sourcePyr = (2*NPyr)*(2*NPyr)*pyrNucSource(c,t);
    real source = sourcePyr; 

    dS[eqn] = 0.0;
    C_UDMI(c,t,2) = source;
    return source;
}

DEFINE_SOURCE (pdf_m_2_NucSourcePyr,c,t,dS,eqn)
{
	/* 
	
	First moment source term for nucleation assuming pyrene as the nucleating species.
	Assumes turbulence interaction by considering gaussian PDFs of temperature and 
	nucleating species. PDFs are assumed to be independent.
	
	*/
	
    real sourcePyr = (2*NPyr)*(2*NPyr)*pdfpyrNucSource(c,t);
    real source = sourcePyr; 

    dS[eqn] = 0.0;
    C_UDMI(c,t,2) = source;
    return source;
}

DEFINE_SOURCE (m_2_NucSourceBenz,c,t,dS,eqn)
{
	/* 
	
	First moment source term for nucleation assuming benzene as the nucleating species.
	Assumes no turbulence interaction.

	*/

    real sourceBenz = (2*NBenzene)*(2*NBenzene)*benzNucSource(c,t); 
    real source = sourceBenz;

    dS[eqn] = 0.0;
    C_UDMI(c,t,2) = source;
    return source;
}