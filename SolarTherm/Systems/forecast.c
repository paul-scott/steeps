// poly.c

#include <math.h>

double st_mip(
    double x,
    double * P_curtail,
    double * P_direct,
    double * P_EES_in,
    double * P_EES_out,
    double * m_h2_stor_in,
    double * m_h2_stor_out
)
{
    P_curtail[0] = 0;
    P_direct[0] = 0;
    P_EES_in[0] = 0;
    P_EES_out[0] = 0;
    m_h2_stor_in[0] = 0;
    m_h2_stor_out[0] = 0;
}

