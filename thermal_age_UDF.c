#include "udf.h"

DEFINE_INIT(thermal_age_setup,domain)
{
	if(NULLP(user_particle_vars)) Init_User_Particle_Vars();
	strcpy(user_particle_vars[0].name, "thermal-age");
	strcpy(user_particle_vars[0].label, "Thermal Age");
	strcpy(user_particle_vars[1].name, "thermal-age-0");
	strcpy(user_particle_vars[1].label, "Thermal Age 0");
}

DEFINE_DPM_SCALAR_UPDATE(thermal_age, cell, thread, initialize, p)
{
	cellTemp = C_T(cell, thread);
	if(initialize)
	{
		P_USER_REAL(p,0) = 0.0;
		P_USER_REAL(p,1) = cellTemp;
	}
	else
	{
		P_USER_REAL(p,0) = P_USER_REAL(p,0) + P_DT(p) * 0.5 * cellTemp;
		P_USER_CELL(p,1) = cellTemp;
	}
}