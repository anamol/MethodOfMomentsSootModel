#include "udf.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"

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
	real particleTemp = P_T(p);
	if(initialize)
	{
		P_USER_REAL(p,0) = 0.0;
		P_USER_REAL(p,1) = particleTemp;
	}
	else
	{
		P_USER_REAL(p,0) = P_USER_REAL(p,0) + P_DT(p) * 0.5 * P_USER_CELL(p,1) * particleTemp;
		P_USER_CELL(p,1) = particleTemp;
	}
}

DEFINE_INIT(init_alpha,d)
{
	FILE* file = fopen("cell_values_alpha.csv", "r");

	char line[256];

	int cell_id[62400]; 
	float thermal_age[62400];
	float thermal_std[62400];
	int i = 0;
	const char s[2] = ",";

	while (fgets(line, sizeof(line), file))
	{
		
		char *token;

		token = strtok(line,s);
		cell_id[i] = atoi(token);

		token = strtok(NULL,s);
		thermal_age[i] = atof(token);

		token = strtok(NULL,s);
		thermal_std[i] = atof(token);

		token = strtok(NULL,s);
		alpha[i] = atof(token);

		i = i + 1;
	}

	fclose(file);

	cell_t c;
	Thread* t;
	real xc[ND_ND];

	thread_loop_c(t,d)
	{
		begin_c_loop_all(c,t)
		{
			int cell_number = c;
			C_UDMI(c,t,30) = thermal_age[cell_number - 1];
			C_UDMI(c,t,31) = thermal_std[cell_number - 1];
			C_UDMI(c,t,32) = alpha[cell_number - 1];

		}
	}
}