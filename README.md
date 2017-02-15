 Soot Method of Moments with Interpolative Closure Model For ANSYS FLUENT 17.2

=====================
*   Anamol Pundle                             *
*   Department of Mechanical Engineering      *
*   University of Washington, Seattle         *


This code implements the Method of Moments with Interpolative Closure by Frenkalch et. al 
(as given in M.Frenklach, Method of Moments with Interpolative Closure, Chem. Eng. Sci.
Volume 57, Issue 12, June 2002, Pages 2229–2239) for use with ANSYS FLUENT 17.2. The 
implementation is for three moments through FLUENT UDF's. The aggregation model and the 
UDF for the absorption coefficient is not yet implemented.  

This code can currently only be used for turbulent non-premixed flames with the 
flamelet model. This file contains:  

1. Source term UDF's for each moment  
    Moment-0 : nucleation, coagulation  
    Moment-1 : nucleation, HACA surface growth, oxidation via O2 and OH  
    Moment-2 : nucleation, coagulation, HACA surface growth, oxidation via O2 and OH  

2. UDF for calculating effective turbulent diffusivity.  
 

RUNNING THE CODE FOR STEADY STATE SOLUTION:  
1. Declare three user defined scalars(UDS) and x user defined memories in Fluent.  
2. Compile the UDF's in this file.  
3. Select only the nucleation source terms for each UDS.  
4. Set the effective turbulent diffusivity of each UDS equal to the eff. turb. diff. UDF.  
5. Solve only the three UDS equations until convergence (switch other equations off).  
6. Add the HACA source term to the first and second moment UDS's and the coagulation source   
   term to the zeroth moment UDS. Solve until convergence, starting with low URF's.  
7. Add the oxidation source terms to the first and second moment UDS's. Solve until
   convergence.  
8. Add the coagulation term to the second moment UDS. Solve until convergence.  



USER DEFINED MEMORY ALLOCATION  
UDM 0 : Nucleation source term for zeroth moment    
UDM 1 : Nucleation source term for first moment  
UDM 2 : Nucleation source term for second moment   
UDM 3 : C2H2 surface growth term for first moment  
UDM 4 : C2H2 surface growth term for second moment  
UDM 5 : O2 oxidation term for first moment  
UDM 6 : OH oxidation term for first moment  
UDM 7 : O2 oxidation term for second moment  
UDM 8 : OH oxidation term for second moment  