#include "udf.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"

#define mesh_size 62400

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
		P_USER_REAL(p,1) = particleTemp;
	}
}

DEFINE_ON_DEMAND(alpha)
{
	FILE* file = fopen("//home/anamol/Fluent/SmallFlameFinerMesh_newPDF/cell_values_alpha.csv", "r");

	char line[256];

	int cell_id[mesh_size]; 
	float thermal_age[mesh_size];
	float thermal_std[mesh_size];
	float alpha[mesh_size];
	float x_coord[mesh_size];
	float y_coord[mesh_size];
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

		token = strtok(NULL,s);
		x_coord[i] = atof(token);

		token = strtok(NULL,s);
		y_coord[i] = atof(token);


		i = i + 1;
	}

	fclose(file);

	Domain* d;
	cell_t c;
	Thread* t;
	real xc[ND_ND];

	d = Get_Domain(1);

	thread_loop_c(t,d)
	{
		begin_c_loop(c,t)
		{
			C_CENTROID(xc,c,t);
			float distance[mesh_size];
			float min_distance = 500.00;
			int index;
			int i;

			for(i = 0; i < mesh_size; i++)
			{
				distance[i] = pow((xc[0] - x_coord[i]),2) + pow((xc[1] - y_coord[i]),2);
				
				if (distance[i] < min_distance) 
				{
					min_distance = distance[i];
					index = i;
				}

			}

			C_UDMI(c,t,30) = thermal_age[index];
			C_UDMI(c,t,31) = thermal_std[index];
			C_UDMI(c,t,32) = alpha[index]; 
			

		}
		end_c_loop(c,t)
	}
}
