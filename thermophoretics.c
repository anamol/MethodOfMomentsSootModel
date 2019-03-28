#include "udf.h"
#include "mom.h"
#include "stdbool.h"

DEFINE_SOURCE(m_0_WallThermophoretic,c,t,dS,eqn)
{
    real cellRho = C_R(c,t);
    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2); 
    real A[ND_ND];
    real faceArea;

    int n;
    int ID1 = 12;
    int ID2 = 13; 
    bool bdyCell = false;
    bool x_grad = false;
    face_t f;
    Thread *tf;
    real source;
    real Vt;

    c_face_loop(c, t, n)
    {
        tf = C_FACE_THREAD(c,t,n);
        f = C_FACE(c,t,n);

        if (BOUNDARY_FACE_THREAD_P(tf))
        {
            F_AREA(A, f, tf);
            faceArea = pow(A[0] * A[0] + A[1] * A[1], 0.5);

            if (THREAD_ID(tf) == ID1)
            {
                bdyCell = true;
                x_grad = false;
            }
            else if (THREAD_ID(tf) == ID2)
            {
                bdyCell = true;
                x_grad = true;
            }

            else 
            {
                bdyCell = false;
            }
            
        }

    }

    if (bdyCell)
    {
        real Kn = C_UDMI(c,t,22);
        real A = 1.2;
        real B = 0.41;
        real C = 0.88;
        real Cs = 1.17;
        real Ct = 2.18;
        real Cm = 1.14;
        real mu = C_MU_L(c,t);
        real nu = mu / cellRho;
        real cellTemp = C_T(c,t);
        real TempGradX = C_T_G(c,t)[0];
        real TempGradY = C_T_G(c,t)[1];
        real ksoot = 0.26;
        real kg = C_K_L(c,t);
        
        real Cc = 1 + Kn * (A + B * exp(-C/Kn));
        real Kt = 2 * Cs * Cc * (kg/ksoot + Ct * Kn) / ((1 + 3 * Cm * Kn) * (1 + 2 * kg/ksoot + Ct * Kn));

        if (x_grad)
        {
            Vt = fabs(mu * Kt  * TempGradX) / cellTemp;
        }
        else 
        {
            Vt = fabs(mu * Kt * TempGradY) / cellTemp;
        }
            
    }

    else
    {
        Vt = 0.0;
    }

    source = -cellM0 * faceArea * Vt / (C_VOLUME(c,t) * 2 * pi);
    C_UDMI(c,t,43) = source;
    return source;
}

DEFINE_SOURCE(m_1_WallThermophoretic,c,t,dS,eqn)
{
    real cellRho = C_R(c,t);
    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2); 
    real A[ND_ND];
    real faceArea;

    int n;
    int ID1 = 12;
    int ID2 = 13; 
    bool bdyCell = false;
    bool x_grad = false;
    face_t f;
    Thread *tf;
    real source;
    real Vt;

    c_face_loop(c, t, n)
    {
        tf = C_FACE_THREAD(c,t,n);
        f = C_FACE(c,t,n);

        if (BOUNDARY_FACE_THREAD_P(tf))
        {
            F_AREA(A, f, tf);
            faceArea = pow(A[0] * A[0] + A[1] * A[1], 0.5);

            if (THREAD_ID(tf) == ID1)
            {
                bdyCell = true;
                x_grad = false;
            }
            else if (THREAD_ID(tf) == ID2)
            {
                bdyCell = true;
                x_grad = true;
            }

            else 
            {
                bdyCell = false;
            }
            
        }

    }

    if (bdyCell)
    {
        real Kn = C_UDMI(c,t,22);
        real A = 1.2;
        real B = 0.41;
        real C = 0.88;
        real Cs = 1.17;
        real Ct = 2.18;
        real Cm = 1.14;
        real mu = C_MU_L(c,t);
        real nu = mu / cellRho;
        real cellTemp = C_T(c,t);
        real TempGradX = C_T_G(c,t)[0];
        real TempGradY = C_T_G(c,t)[1];
        real ksoot = 0.26;
        real kg = C_K_L(c,t);
        
        real Cc = 1 + Kn * (A + B * exp(-C/Kn));
        real Kt = 2 * Cs * Cc * (kg/ksoot + Ct * Kn) / ((1 + 3 * Cm * Kn) * (1 + 2 * kg/ksoot + Ct * Kn));

        if (x_grad)
        {
            Vt = fabs(mu * Kt * TempGradX) / cellTemp;
        }
        else 
        {
            Vt = fabs(mu * Kt * TempGradY) / cellTemp;
        }
            
    }

    else
    {
        Vt = 0.0;
    }

    source = -cellM1 * faceArea * Vt / (C_VOLUME(c,t) * 2 * pi);
    C_UDMI(c,t,44) = source;
    return source;
}

DEFINE_SOURCE(m_2_WallThermophoretic,c,t,dS,eqn)
{
    real cellRho = C_R(c,t);
    real cellM0 = cellRho*C_UDSI(c,t,0);
    real cellM1 = cellRho*C_UDSI(c,t,1);
    real cellM2 = cellRho*C_UDSI(c,t,2); 
    real A[ND_ND];
    real faceArea;

    int n;
    int ID1 = 12;
    int ID2 = 13; 
    bool bdyCell = false;
    bool x_grad = false;
    face_t f;
    Thread *tf;
    real source;
    real Vt;

    c_face_loop(c, t, n)
    {
        tf = C_FACE_THREAD(c,t,n);
        f = C_FACE(c,t,n);

        if (BOUNDARY_FACE_THREAD_P(tf))
        {
            F_AREA(A, f, tf);
            faceArea = pow(A[0] * A[0] + A[1] * A[1], 0.5);

            if (THREAD_ID(tf) == ID1)
            {
                bdyCell = true;
                x_grad = false;
            }
            else if (THREAD_ID(tf) == ID2)
            {
                bdyCell = true;
                x_grad = true;
            }

            else 
            {
                bdyCell = false;
            }
            
        }

    }

    if (bdyCell)
    {
        real Kn = C_UDMI(c,t,22);
        real A = 1.2;
        real B = 0.41;
        real C = 0.88;
        real Cs = 1.17;
        real Ct = 2.18;
        real Cm = 1.14;
        real mu = C_MU_L(c,t);
        real nu = mu / cellRho;
        real cellTemp = C_T(c,t);
        real TempGradX = C_T_G(c,t)[0];
        real TempGradY = C_T_G(c,t)[1];
        real ksoot = 0.26;
        real kg = C_K_L(c,t);
        
        real Cc = 1 + Kn * (A + B * exp(-C/Kn));
        real Kt = 2 * Cs * Cc * (kg/ksoot + Ct * Kn) / ((1 + 3 * Cm * Kn) * (1 + 2 * kg/ksoot + Ct * Kn));

        if (x_grad)
        {
            Vt = fabs(mu * Kt * TempGradX) / cellTemp;
        }
        else 
        {
            Vt = fabs(mu * Kt * TempGradY) / cellTemp;
        }
            
    }

    else
    {
        Vt = 0.0;
    }

    source = -cellM2 * faceArea * Vt / (C_VOLUME(c,t) * 2 * pi);

    if (source > 1.0e20)
    {
        source = 1.0e19;
    }
    C_UDMI(c,t,45) = source;
    return source;
}