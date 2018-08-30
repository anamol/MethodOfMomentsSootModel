#include "udf.h"
#include "mom.h"

DEFINE_DIFFUSIVITY(uds_diff,c,t,i)
{
    real diff = C_MU_EFF(c,t)/(ScT*PrT); 

    return diff;
}

DEFINE_PROPERTY(abs_coeff_Widmann,c,t)
{
    real cellTemp = C_T(c,t);
    real cellRho = C_R(c,t);
    real cellM1 = C_UDSI(c,t,1);
    real sootfv = cellM1*mC*cellRho/rhoSoot*normParameter;

    real absCoeff = 2370*cellTemp*sootfv;
    return absCoeff;
}

DEFINE_PROPERTY(abs_coeff_Sazhin,c,t)
{
    real cellTemp = C_T(c,t);
    real cellRho = C_R(c,t);
    real cellM1 = C_UDSI(c,t,1);
    real sootfv = cellM1*mC*cellRho/rhoSoot*normParameter;

    real absCoeff = 1232.4*rhoSoot*sootfv*(1+4.8e-4*(cellTemp-2000));
    return absCoeff;
}

DEFINE_WSGGM_ABS_COEFF(wsggm_abs_coeff_Sazhin, c, t, xi, p_t, s, soot_conc, Tcell, nb, ab_wsggm, ab_soot)
{
	real cellTemp = C_T(c,t);
    real cellRho = C_R(c,t);
    real cellM1 = C_UDSI(c,t,1);
    real sootfv = cellM1*mC*cellRho/rhoSoot*normParameter;

    real absCoeff = 1232.4*rhoSoot*sootfv*(1+4.8e-4*(cellTemp-2000));

    
    *ab_wsggm = *ab_wsggm + absCoeff;
    
}

DEFINE_WSGGM_ABS_COEFF(wsggm_abs_coeff_Widmann, c, t, xi, p_t, s, soot_conc, Tcell, nb, ab_wsggm, ab_soot)
{
	real cellTemp = C_T(c,t);
    real cellRho = C_R(c,t);
    real cellM1 = C_UDSI(c,t,1);
    real sootfv = cellM1*mC*cellRho/rhoSoot*normParameter;

    real absCoeff = 2370*cellTemp*sootfv;

    *ab_wsggm = *ab_wsggm + absCoeff;
}

DEFINE_PROPERTY(thermal_conductivity,c,t)
{
    real cellTemp = C_T(c,t);
    real Tr = cellTemp/132.52;


    real thermCond = 4.358e-3*(33.972*pow(Tr,-1) - 164.702*pow(Tr,-twoThird) + 262.108*pow(Tr,-oneThird) 
        - 21.534 - 443.455*pow(Tr,oneThird) + 607.339*pow(Tr, twoThird) - 368.79*Tr + 111.296*pow(Tr, 2*twoThird) 
        - 13.412*pow(Tr, 5*oneThird));

    return thermCond;
}