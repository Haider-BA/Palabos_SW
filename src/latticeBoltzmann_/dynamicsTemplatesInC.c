

#include "latticeBoltzmann_/dynamicsTemplatesInC.h"

double bgk_ma2_collision_InC(struct UniformSequence cellPopulation, double rhoBar, struct UniformSequence_const j_, double omega, struct DescriptorInC* descriptor_)
{
    double one_m_omega = (double)1 - omega;
    double t0_omega = (*descriptor_).t[0] * omega;
    double t1_omega = (*descriptor_).t[1] * omega;
    double t4_omega = (*descriptor_).t[4] * omega;

    double* f = (double*) (cellPopulation.sequence);
    const double* j= (const double*) (j_.sequence);
    double invRho = (*descriptor_).invRho(rhoBar);


    double jSqr   = j[0]*j[0] + j[1]*j[1] + j[2]*j[2];
    double kx     = (double)3 * j[0];
    double ky     = (double)3 * j[1];
    double kz     = (double)3 * j[2];
    double kxSqr_ = invRho / (double)2 * kx*kx;
    double kySqr_ = invRho / (double)2 * ky*ky;
    double kzSqr_ = invRho / (double)2 * kz*kz;
    double kxky_  = invRho * kx*ky;
    double kxkz_  = invRho * kx*kz;
    double kykz_  = invRho * ky*kz;

    double C1 = rhoBar + invRho*(double)3*jSqr;
    double C2, C3;

    // i=0
    C3 = -kxSqr_ - kySqr_ - kzSqr_;
    f[0] *= one_m_omega; f[0] += t0_omega * (C1+C3);

    // i=1 and i=10
    C2 = -kx;
    C3 = -kySqr_ - kzSqr_;
    f[1]  *= one_m_omega; f[1]  += t1_omega * (C1+C2+C3);
    f[10] *= one_m_omega; f[10] += t1_omega * (C1-C2+C3);

    // i=2 and i=11
    C2 = -ky;
    C3 = -kxSqr_ - kzSqr_;
    f[2]  *= one_m_omega; f[2]  += t1_omega * (C1+C2+C3);
    f[11] *= one_m_omega; f[11] += t1_omega * (C1-C2+C3);

    // i=3 and i=12
    C2 = -kz;
    C3 = -kxSqr_ - kySqr_;
    f[3]  *= one_m_omega; f[3]  += t1_omega * (C1+C2+C3);
    f[12] *= one_m_omega; f[12] += t1_omega * (C1-C2+C3);

    // i=4 and i=13
    C2 = -kx - ky;
    C3 = kxky_ - kzSqr_;
    f[4]  *= one_m_omega; f[4]  += t4_omega * (C1+C2+C3);
    f[13] *= one_m_omega; f[13] += t4_omega * (C1-C2+C3);

    // i=5 and i=14
    C2 = -kx + ky;
    C3 = -kxky_ - kzSqr_;
    f[5]  *= one_m_omega; f[5]  += t4_omega * (C1+C2+C3);
    f[14] *= one_m_omega; f[14] += t4_omega * (C1-C2+C3);

    // i=6 and i=15
    C2 = -kx - kz;
    C3 = kxkz_ - kySqr_;
    f[6]  *= one_m_omega; f[6]  += t4_omega * (C1+C2+C3);
    f[15] *= one_m_omega; f[15] += t4_omega * (C1-C2+C3);

    // i=7 and i=16
    C2 = -kx + kz;
    C3 = -kxkz_ - kySqr_;
    f[7]  *= one_m_omega; f[7]  += t4_omega * (C1+C2+C3);
    f[16] *= one_m_omega; f[16] += t4_omega * (C1-C2+C3);

    // i=8 and i=17
    C2 = -ky - kz;
    C3 = kykz_ - kxSqr_;
    f[8]  *= one_m_omega; f[8]  += t4_omega * (C1+C2+C3);
    f[17] *= one_m_omega; f[17] += t4_omega * (C1-C2+C3);

    // i=9 and i=18
    C2 = -ky + kz;
    C3 = -kykz_ - kxSqr_;
    f[9]  *= one_m_omega; f[9]  += t4_omega * (C1+C2+C3);
    f[18] *= one_m_omega; f[18] += t4_omega * (C1-C2+C3);

    return invRho*invRho*jSqr;
}

double bgk_ma2_equilibrium_InC (int iPop, double rhoBar, double invRho, struct UniformSequence_const j, double jSqr, struct DescriptorInC* descriptor_)
{
    ///
    double* jPop= (double*) (j.sequence);
    double c_j = (*descriptor_).c[iPop][0]*(*jPop);
    for (int iD=1; iD < (*descriptor_).d; ++iD) {
    c_j +=(*descriptor_).c[iPop][iD]*(*jPop+iD);
    }
    return (*descriptor_).t[iPop] * (
    rhoBar + (*descriptor_).invCs2 * c_j +
    (*descriptor_).invCs2/(double)2 * invRho * (
    (*descriptor_).invCs2 * c_j*c_j - jSqr )
    );
}
