#include "math.h"
#include "K_Golden.h"
//
void
DIFF_ps_t_K(double t, double & ps, double & dpsdt, double & d2psdt2)
{ // units: R, atm
  static const double a = 1.3408e6;
  static const double b = 0.53299;
  static const double c = -18717.;
  //
    ps = (a / pow(t, b)) * exp(c / t);
    dpsdt = -ps * (c / (t * t) + b / t);
    d2psdt2 = -dpsdt * (c / (t * t) + b / t) - ps * (-2. * c / (t * t * t) - b / (t * t));
}
//
void
DIFF_ps_t_K(double t, double & ps, double & dpsdt, double & d2psdt2, double & d3psdt3)
{ // units: R, atm
  static const double a = 1.3408e6;
  static const double b = 0.53299;
  static const double c = -18717.;
  //
    ps = (a / pow(t, b)) * exp(c / t);
    dpsdt = -ps * (c / (t * t) + b / t);
    d2psdt2 = -dpsdt * (c / (t * t) + b / t) - ps * (-2. * c / (t * t * t) - b / (t * t));
    d3psdt3 = -d2psdt2 * (c / (t * t) + b / t) - dpsdt * (-2. * c / (t * t * t) - b / (t * t)) -
              dpsdt * (-2. * c / (t * t * t) - b / (t * t)) -
              ps * (6. * c / (t * t * t * t) + 2. * b / (t * t * t));
}
//
int
DIFF_ts_p_K(double p, double & ts, double & dtsdp, double & d2tsdp2)
{ // units: R, atm
  double ps, dpsdt, d2psdt2;

  static const double tol_p = 1.e-6;
  static const double tnb = 1029.73 * 9. / 5.;
  // initial guess from normal boiling point in R
  ts = tnb;

  int it = 0;
  double dp = 1.;
  while (fabs(dp) / p > tol_p)
  {
    DIFF_ps_t_K(ts, ps, dpsdt, d2psdt2);
    dp = ps - p;
    ts -= dp / dpsdt;
    if (++it > 20)
    {
      return -1;
    }
  }
  DIFF_ps_t_K(ts, ps, dpsdt, d2psdt2);
  dtsdp = 1. / dpsdt;
  d2tsdp2 = -d2psdt2 * dtsdp * dtsdp * dtsdp;
  return 0;
}
//
void
DIFF_d1_t_K(double t, double & d1, double & dd1dt, double & d2d1dt2)
{ // units: R, lb/ft3
  double tt;
  //
  static const double a = 52.768;
  static const double b = -7.4975e-3;
  static const double c = -0.5255e-6;
  static const double d = 0.0498e-9;
  //
  tt = t - 459.7;
  d1 = a + tt * (b + tt * (c + tt * d));
  dd1dt = b + tt * (2. * c + tt * 3. * d);
  d2d1dt2 = 2. * c + tt * 6. * d;
}
//
void
DIFF_d1_t_K(double t, double & d1, double & dd1dt, double & d2d1dt2, double & d3d1dt3)
{ // units: R, lb/ft3
  double tt;
  //
  static const double a = 52.768;
  static const double b = -7.4975e-3;
  static const double c = -0.5255e-6;
  static const double d = 0.0498e-9;
  //
  tt = t - 459.7;
  d1 = a + tt * (b + tt * (c + tt * d));
  dd1dt = b + tt * (2. * c + tt * 3. * d);
  d2d1dt2 = 2. * c + tt * 6. * d;
  d3d1dt3 = 6. * d;
}
//
void
DIFF_v1_t_K(double t, double & v1, double & dv1dt, double & d2v1dt2)
{ // units: R, ft3/lb
  double d1, dd1dt, d2d1dt2;
  //
  DIFF_d1_t_K(t, d1, dd1dt, d2d1dt2);
  //
  v1 = 1. / d1;
  dv1dt = -dd1dt * v1 * v1;
  d2v1dt2 = v1 * v1 * (dd1dt * dd1dt * 2. * v1 - d2d1dt2);
}
//
//-----------------------------------------------------------------------------
// isothermal compressibility in the liquid phase
// The calculation relies on a corresponding states law established by
// Pasternak. The isothermal compressibility of sodium is calculated for a
// transformed temperature and the result is transformed again to obtain the
// isothermal compressibility of potassium. The correlation for sodium is
// a polynomial fit to the values derived by Fink and Leibowitz. The
// temperature range of the data is larger than that of Pasternaks evaluation.
// Both sets of values are in fairly good agreement.
//-----------------------------------------------------------------------------
//
void
betat_t_K(double t, double & betat, double & dbetatdt, double & d2betatdt2)
{ // R, 1/atm
  // betat equation in 1/MPa as a function of T*T (T in K)
  double dbetatdt_i, d2betatdt2_i;
  double dbetatdt_int, d2betatdt2_int;
  //Na coefficients
  static const double a6 =  3.10913961e-40;
  static const double a5 =- 3.02303007e-33;
  static const double a4 =  3.16640457e-27;
  static const double a3 =  7.18050234e-20;
  static const double a2 =- 3.27404587e-13;
  static const double a1 =  1.04346122e-06;
  static const double a0 =- 8.70307432;
  //atomic diameters (A) and energy parameters (K)
  static const double d_Na = 3.84;
  static const double d_K  = 4.76;
  static const double e_Na = 1970.;
  static const double e_K  = 1760.;
  //derivatives of reduced parameters
  static const double dbetaK_dbetaNa=(e_Na/e_K)*(d_K/d_Na)*(d_K/d_Na)*(d_K/d_Na);
  static const double dTNa_dTK=e_Na/e_K;
  //transformation of temperature (to K and for sodium)
  double tK = t * (5. / 9.) * dTNa_dTK;
  double tK2 = tK * tK;
  double dtK2dt = 2. * tK * (5. / 9.) * dTNa_dTK;
  double d2tK2dt2 = 2. * (5. / 9.) * (5. / 9.) * dTNa_dTK * dTNa_dTK;
  //
  betat =
      (((((a6 * tK2 + a5) * tK2 + a4) * tK2 + a3) *
            tK2 +
        a2) *
           tK2 +
       a1) *
          tK2 +
      a0;
  dbetatdt_i =
      ((((6. * a6 * tK2 + 5. * a5) * tK2 + 4. * a4) * tK2 +
        3. * a3) *
           tK2 +
       2. * a2) *
          tK2 +
      a1;
  d2betatdt2_i =
      (((30. * a6 * tK2 + 20. * a5) * tK2 + 12. * a4) * tK2 +
       6. * a3) *
          tK2 +
      2. * a2;
  betat = exp(betat);
  dbetatdt_int = betat * dbetatdt_i;
  d2betatdt2_int = dbetatdt_int * dbetatdt_i + betat * d2betatdt2_i;
  betat *= 0.101325;
  dbetatdt = 0.101325 * dbetatdt_int * dtK2dt;
  d2betatdt2 = 0.101325 * (d2betatdt2_int * dtK2dt * dtK2dt + dbetatdt_int * d2tK2dt2);
  //corresponding states law (temperature transformation already considered)
  betat *= dbetaK_dbetaNa;
  dbetatdt *= dbetaK_dbetaNa;
  d2betatdt2 *= dbetaK_dbetaNa;
}
//
void
betat_t_K(double t, double & betat, double & dbetatdt, double & d2betatdt2, double & d3betatdt3)
{ // R, 1/atm
  // betat equation in 1/MPa as a function of T*T (T in K)
  double dbetatdt_i, d2betatdt2_i, d3betatdt3_i;
  double dbetatdt_int, d2betatdt2_int, d3betatdt3_int;
  //Na coefficients
  static const double a6 =  3.10913961e-40;
  static const double a5 =- 3.02303007e-33;
  static const double a4 =  3.16640457e-27;
  static const double a3 =  7.18050234e-20;
  static const double a2 =- 3.27404587e-13;
  static const double a1 =  1.04346122e-06;
  static const double a0 =- 8.70307432;
  //atomic diameters (A) and energy parameters (K)
  static const double d_Na = 3.84;
  static const double d_K  = 4.76;
  static const double e_Na = 1970.;
  static const double e_K  = 1760.;
  //derivatives of reduced parameters
  static const double dbetaK_dbetaNa=(e_Na/e_K)*(d_K/d_Na)*(d_K/d_Na)*(d_K/d_Na);
  static const double dTNa_dTK=e_Na/e_K;
  //transformation of temperature (to K and for sodium)
  double tK = t * (5. / 9.) * dTNa_dTK;
  double tK2 = tK * tK;
  double dtK2dt = 2. * tK * (5. / 9.) * dTNa_dTK;
  double d2tK2dt2 = 2. * (5. / 9.) * (5. / 9.) * dTNa_dTK * dTNa_dTK;
  //
  betat =
      (((((a6 * tK2 + a5) * tK2 + a4) * tK2 + a3) *
            tK2 +
        a2) *
           tK2 +
       a1) *
          tK2 +
      a0;
  dbetatdt_i =
      ((((6. * a6 * tK2 + 5. * a5) * tK2 + 4. * a4) * tK2 +
        3. * a3) *
           tK2 +
       2. * a2) *
          tK2 +
      a1;
  d2betatdt2_i =
      (((30. * a6 * tK2 + 20. * a5) * tK2 + 12. * a4) * tK2 +
       6. * a3) *
          tK2 +
      2. * a2;
  d3betatdt3_i =
      ((120. * a6 * tK2 + 60. * a5) * tK2 + 24. * a4) * tK2 +
      6. * a3;
  betat = exp(betat);
  dbetatdt_int = betat * dbetatdt_i;
  d2betatdt2_int = dbetatdt_int * dbetatdt_i + betat * d2betatdt2_i;
  d3betatdt3_int = d2betatdt2_int * dbetatdt_i + dbetatdt_int * d2betatdt2_i +
                   dbetatdt_int * d2betatdt2_i + betat * d3betatdt3_i;
  betat *= 0.101325;
  dbetatdt = 0.101325 * dbetatdt_int * dtK2dt;
  d2betatdt2 = 0.101325 * (d2betatdt2_int * dtK2dt * dtK2dt + dbetatdt_int * d2tK2dt2);
  d3betatdt3 =
      0.101325 * (d3betatdt3_int * dtK2dt * dtK2dt * dtK2dt +
                  d2betatdt2_int * 2. * dtK2dt * d2tK2dt2 + d2betatdt2_int * dtK2dt * d2tK2dt2);
  //corresponding states law (temperature transformation already considered)
  betat *= dbetaK_dbetaNa;
  dbetatdt *= dbetaK_dbetaNa;
  d2betatdt2 *= dbetaK_dbetaNa;
  d3betatdt3 *= dbetaK_dbetaNa;
}
//
void
DIFF_d_tp_L_K(double t,
               double p,
               double & dl,
               double & ddldt,
               double & d2dldt2,
               double & ddldp,
               double & d2dldp2,
               double & d2dldtdp)
{ // units: R, atm, lb/ft3
  double d1, dd1dt, d2d1dt2, ps, dpsdt, d2psdt2, betat, dbetatdt, d2betatdt2;
  //
  DIFF_d1_t_K(t, d1, dd1dt, d2d1dt2);
  DIFF_ps_t_K(t, ps, dpsdt, d2psdt2);
  betat_t_K(t, betat, dbetatdt, d2betatdt2);
  //
  dl = d1 + d1 * betat * (p - ps);
  ddldt = dd1dt + dd1dt * betat * (p - ps) + d1 * dbetatdt * (p - ps) - d1 * betat * dpsdt;
  d2dldt2 = d2d1dt2 + d2d1dt2 * betat * (p - ps) + dd1dt * dbetatdt * (p - ps) -
            dd1dt * betat * dpsdt + dd1dt * dbetatdt * (p - ps) + d1 * d2betatdt2 * (p - ps) -
            d1 * dbetatdt * dpsdt - dd1dt * betat * dpsdt - d1 * dbetatdt * dpsdt -
            d1 * betat * d2psdt2;
  ddldp = d1 * betat;
  d2dldp2 = 0.;
  d2dldtdp = dd1dt * betat + d1 * dbetatdt;
}
//
void
DIFF_d_tp_L_K(double t,
               double p,
               double & dl,
               double & ddldt,
               double & d2dldt2,
               double & d3dldt3,
               double & ddldp,
               double & d2dldp2,
               double & d2dldtdp,
               double & d3dldt2dp,
               double & d3dldtdp2)
{ // units: R, atm, lb/ft3
  double d1, dd1dt, d2d1dt2, d3d1dt3, ps, dpsdt, d2psdt2, d3psdt3, betat, dbetatdt, d2betatdt2,
      d3betatdt3;
  //
  DIFF_d1_t_K(t, d1, dd1dt, d2d1dt2, d3d1dt3);
  DIFF_ps_t_K(t, ps, dpsdt, d2psdt2, d3psdt3);
  betat_t_K(t, betat, dbetatdt, d2betatdt2, d3betatdt3);
  //
  dl = d1 + d1 * betat * (p - ps);
  ddldt = dd1dt + dd1dt * betat * (p - ps) + d1 * dbetatdt * (p - ps) - d1 * betat * dpsdt;
  d2dldt2 = d2d1dt2 + d2d1dt2 * betat * (p - ps) + 2. * dd1dt * dbetatdt * (p - ps) -
            2. * dd1dt * betat * dpsdt + d1 * d2betatdt2 * (p - ps) - 2. * d1 * dbetatdt * dpsdt -
            d1 * betat * d2psdt2;
  d3dldt3 = d3d1dt3 + d3d1dt3 * betat * (p - ps) + d2d1dt2 * dbetatdt * (p - ps) -
            d2d1dt2 * betat * dpsdt + 2. * d2d1dt2 * dbetatdt * (p - ps) +
            2. * dd1dt * d2betatdt2 * (p - ps) - 2. * dd1dt * dbetatdt * dpsdt -
            2. * d2d1dt2 * betat * dpsdt - 2. * dd1dt * dbetatdt * dpsdt -
            2. * dd1dt * betat * d2psdt2 + dd1dt * d2betatdt2 * (p - ps) +
            d1 * d3betatdt3 * (p - ps) - d1 * d2betatdt2 * dpsdt - 2. * dd1dt * dbetatdt * dpsdt -
            2. * d1 * d2betatdt2 * dpsdt - 2. * d1 * dbetatdt * d2psdt2 - dd1dt * betat * d2psdt2 -
            d1 * dbetatdt * d2psdt2 - d1 * betat * d3psdt3;
  ddldp = d1 * betat;
  d2dldp2 = 0.;
  d2dldtdp = dd1dt * betat + d1 * dbetatdt;
  d3dldtdp2 = 0.0;
  d3dldt2dp = d2d1dt2 * betat + 2. * dd1dt * dbetatdt + d1 * d2betatdt2;
}
//
void
DIFF_v_tp_L_K(double t,
               double p,
               double & vl,
               double & dvldt,
               double & d2vldt2,
               double & dvldp,
               double & d2vldp2,
               double & d2vldtdp)
{ // units: R, atm, lb/ft3
  double dl, ddldt, d2dldt2, ddldp, d2dldp2, d2dldtdp;
  //
  DIFF_d_tp_L_K(t, p, dl, ddldt, d2dldt2, ddldp, d2dldp2, d2dldtdp);
  //
  vl = 1. / dl;
  dvldt = -ddldt * vl * vl;
  d2vldt2 = -d2dldt2 * vl * vl - 2. * ddldt * vl * dvldt;
  dvldp = -ddldp * vl * vl;
  d2vldp2 = -d2dldp2 * vl * vl - 2. * ddldp * vl * dvldp;
  d2vldtdp = -d2dldtdp * vl * vl - 2. * ddldt * vl * dvldp;
}
//
void
DIFF_h1_t_K(double t, double & h1, double & dh1dt, double & d2h1dt2)
{ // units: R, Btu/lb
  static const double a = 87.8783;
  static const double b = 0.2022;
  static const double c = -0.2177e-4;
  static const double d = 0.07741e-7;
  //
  h1 = a + t * (b + t * (c + t * d));
  dh1dt = b + t * (2. * c + t * 3. * d);
  d2h1dt2 = 2. * c + t * 6. * d;
}
//
void
DIFF_vu_tp_L_K(double t,
                double p,
                double & vl,
                double & dvldt,
                double & d2vldt2,
                double & dvldp,
                double & d2vldp2,
                double & d2vldtdp,
                double & ul,
                double & duldt,
                double & d2uldt2,
                double & duldp,
                double & d2uldp2,
                double & d2uldtdp)
{ // units: R, atm, Btu/lb
  static const double atmft3_toBtu = (101325. * 0.3048 * 0.3048 * 0.3048) / 1055.05585262;

  double h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  DIFF_v_tp_L_K(t, p, vl, dvldt, d2vldt2, dvldp, d2vldp2, d2vldtdp);
  DIFF_h_tp_L_K(t, p, h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);

  ul = h - p * vl * atmft3_toBtu;
  duldt = dhdt - p * dvldt * atmft3_toBtu;
  d2uldt2 = d2hdt2 - p * d2vldt2 * atmft3_toBtu;
  duldp = dhdp - (vl + p * dvldp) * atmft3_toBtu;
  d2uldp2 = d2hdp2 - (dvldp + dvldp + p * d2vldp2) * atmft3_toBtu;
  d2uldtdp = d2hdtdp - (dvldt + p * d2vldtdp) * atmft3_toBtu;
}
//
void
DIFF_h_tp_L_K(double t,
               double p,
               double & hl,
               double & dhldt,
               double & d2hldt2,
               double & dhldp,
               double & d2hldp2,
               double & d2hldtdp)
{ // units: R, atm, Btu/lb
  double ps, dpsdt, d2psdt2, h1, dh1dt, d2h1dt2, dlp, ddlpdt, d2dlpdt2, ddlpdp, d2dlpdp2, d2dlpdtdp;
  double dl, ddldt, d2dldt2, d3dldt3, ddldp, d2dldp2, d2dldtdp, d3dldt2dp, d3dldtdp2;
  double term0, dterm0dt, d2term0dt2, dterm0dp, d2term0dp2, d2term0dtdp, term1, dterm1dt,
      d2term1dt2, dterm1dp, d2term1dp2, d2term1dtdp, term2;
  static const double b = 2.721308; //this conversion factor is taken from Golden et al., a more precise value would be 2.7194811 (remains unchanged to be consistent with Golden's paper)
  //
  DIFF_ps_t_K(t, ps, dpsdt, d2psdt2);
  DIFF_h1_t_K(t, h1, dh1dt, d2h1dt2);
  DIFF_d_tp_L_K(t, p, dl, ddldt, d2dldt2, d3dldt3, ddldp, d2dldp2, d2dldtdp, d3dldt2dp, d3dldtdp2);
  //
  dlp = ddldt;
  ddlpdt = d2dldt2;
  d2dlpdt2 = d3dldt3;
  ddlpdp = d2dldtdp;
  d2dlpdp2 = d3dldtdp2;
  d2dlpdtdp = d3dldt2dp;
  term0 = 1. / dl;
  dterm0dt = -term0 * term0 * ddldt;
  d2term0dt2 = -2. * term0 * dterm0dt * ddldt - term0 * term0 * d2dldt2;
  dterm0dp = -term0 * term0 * ddldp;
  d2term0dp2 = -2. * term0 * dterm0dp * ddldp - term0 * term0 * d2dldp2;
  d2term0dtdp = -2. * term0 * dterm0dp * ddldt - term0 * term0 * d2dldtdp;
  term1 = 1. + t * term0 * dlp;
  dterm1dt = dlp * (term0 + t * dterm0dt) + t * term0 * ddlpdt;
  d2term1dt2 = ddlpdt * (term0 + t * dterm0dt) + dlp * (2. * dterm0dt + t * d2term0dt2) +
               ddlpdt * (term0 + t * dterm0dt) + t * term0 * d2dlpdt2;
  dterm1dp = t * (dterm0dp * dlp + term0 * ddlpdp);
  d2term1dp2 = t * (d2term0dp2 * dlp + 2. * dterm0dp * ddlpdp + term0 * d2dlpdp2);
  d2term1dtdp = ddlpdp * (term0 + t * dterm0dt) + dlp * (dterm0dp + t * d2term0dtdp) +
                t * (dterm0dp * ddlpdt + term0 * d2dlpdtdp);
  term2 = p - ps;
  hl = h1 + term0 * term1 * term2 * b;
  dhldt = dh1dt + (dterm0dt * term1 * term2 + term0 * (dterm1dt * term2 - term1 * dpsdt)) * b;
  d2hldt2 = d2h1dt2 +
            (d2term0dt2 * term1 * term2 + dterm0dt * dterm1dt * term2 - dterm0dt * term1 * dpsdt +
             dterm0dt * (dterm1dt * term2 - term1 * dpsdt) +
             term0 * (d2term1dt2 * term2 - dterm1dt * dpsdt - dterm1dt * dpsdt - term1 * d2psdt2)) *
                b;
  dhldp = (dterm0dp * term1 * term2 + term0 * dterm1dp * term2 + term0 * term1) * b;
  d2hldp2 = (d2term0dp2 * term1 * term2 + dterm0dp * dterm1dp * term2 + dterm0dp * term1 +
             dterm0dp * dterm1dp * term2 + term0 * d2term1dp2 * term2 + term0 * dterm1dp +
             dterm0dp * term1 + term0 * dterm1dp) *
            b;
  d2hldtdp = (d2term0dtdp * term1 * term2 + dterm0dt * dterm1dp * term2 + dterm0dt * term1 +
              dterm0dp * (dterm1dt * term2 - term1 * dpsdt) +
              term0 * (d2term1dtdp * term2 + dterm1dt - dterm1dp * dpsdt)) *
             b;
}
//
void
DIFF_s1_t_K(double t, double & s1, double & ds1dt, double & d2s1dt2)
{ // units: R, Btu/(lb*F)
  static const double a = 0.227126;
  static const double b = -0.64848e-4;
  static const double c = 0.11589e-7;
  static const double d = -0.9646;
  //
  s1 = a * log(t) + t * (b + c * t) + d;
  ds1dt = a / t + b + 2. * c * t;
  d2s1dt2 = -a / (t * t) + 2. * c;
}
//
void
DIFF_s_tp_L_K(double t,
               double p,
               double & sl,
               double & dsldt,
               double & d2sldt2,
               double & dsldp,
               double & d2sldp2,
               double & d2sldtdp)
{ // units: R, atm, Btu/lb
  double ps, dpsdt, d2psdt2, s1, ds1dt, d2s1dt2, dlp, ddlpdt, d2dlpdt2, ddlpdp, d2dlpdp2, d2dlpdtdp;
  double dl, ddldt, d2dldt2, d3dldt3, ddldp, d2dldp2, d2dldtdp, d3dldt2dp, d3dldtdp2;
  double dl_inv, term0, dterm0dt, d2term0dt2, dterm0dp, d2term0dp2, d2term0dtdp, term2;
  static const double b = 2.721308; //this conversion factor is taken from Golden et al., a more precise value would be 2.7194811 (remains unchanged to be consistent with Golden's paper)
  //
  DIFF_ps_t_K(t, ps, dpsdt, d2psdt2);
  DIFF_s1_t_K(t, s1, ds1dt, d2s1dt2);
  DIFF_d_tp_L_K(t, p, dl, ddldt, d2dldt2, d3dldt3, ddldp, d2dldp2, d2dldtdp, d3dldt2dp, d3dldtdp2);
  //
  dlp = ddldt;
  ddlpdt = d2dldt2;
  d2dlpdt2 = d3dldt3;
  ddlpdp = d2dldtdp;
  d2dlpdp2 = d3dldtdp2;
  d2dlpdtdp = d3dldt2dp;

  dl_inv = 1. / dl;
  term0 = dl_inv * dl_inv;
  dterm0dt = -2. * term0 * dl_inv * ddldt;
  d2term0dt2 = 6. * term0 * term0 * term0 * ddldt * ddldt - 2. * term0 * dl_inv * d2dldt2;
  dterm0dp = -2. * term0 * dl_inv * ddldp;
  d2term0dp2 = 6. * term0 * term0 * term0 * ddldp * ddldp - 2. * term0 * dl_inv * d2dldp2;
  d2term0dtdp = 6. * term0 * term0 * term0 * ddldp * ddldt - 2. * term0 * dl_inv * d2dldtdp;
  term2 = p - ps;
  sl = s1 + term0 * dlp * term2 * b;
  dsldt = ds1dt + (dterm0dt * dlp * term2 + term0 * (ddlpdt * term2 - dlp * dpsdt)) * b;
  d2sldt2 =
      d2s1dt2 + (d2term0dt2 * dlp * term2 + dterm0dt * ddlpdt * term2 - dterm0dt * dlp * dpsdt +
                 dterm0dt * (ddlpdt * term2 - dlp * dpsdt) +
                 term0 * (d2dlpdt2 * term2 - ddlpdt * dpsdt - ddlpdt * dpsdt - dlp * d2psdt2)) *
                    b;
  dsldp = (dterm0dp * dlp * term2 + term0 * (ddlpdp * term2 + dlp)) * b;
  d2sldp2 = (d2term0dp2 * dlp * term2 + dterm0dp * ddlpdp * term2 + dterm0dp * dlp +
             dterm0dp * (ddlpdp * term2 + dlp) + term0 * (d2dlpdp2 * term2 + ddlpdp + ddlpdp)) *
            b;
  d2sldtdp = (d2term0dtdp * dlp * term2 + dterm0dt * ddlpdp * term2 + dterm0dt * dlp +
              dterm0dp * (ddlpdt * term2 - dlp * dpsdt) +
              term0 * (d2dlpdtdp * term2 + ddlpdt - ddlpdp * dpsdt)) *
             b;
}
//
//
//
void
DIFF_qc_t_K(double t,
             double ps,
             double dpsdt,
             double d2psdt2,
             double & x1,
             double & dx1dt,
             double & d2x1dt2,
             double & x2,
             double & dx2dt,
             double & d2x2dt2,
             double & x4,
             double & dx4dt,
             double & d2x4dt2,
             double & abar,
             double & dabardt,
             double & d2abardt2,
             double & bd,
             double & dbddt,
             double & d2bddt2,
             double & b2,
             double & db2dt,
             double & d2b2dt2,
             double & b4,
             double & db4dt,
             double & d2b4dt2)
{ //! units: R
  //
  double two, four, u, w, x1p, dwdt, d2wdt2;
  double dtwodt, d2twodt2, dfourdt, d2fourdt2, dudt, d2udt2;
  double num, dnumdt, d2numdt2, den, ddendt, d2dendt2;
  double x1_p2, x1_p3, x1_p4;
  //
  static const double a_2 = -8.9033;
  static const double b_2 = 12250.1;
  static const double a_4 = -23.394;
  static const double b_4 = 31694.6;
  static const double c = 39.102;
  //
  two = exp(a_2 + (b_2 / t));
  dtwodt = two * (-b_2 / (t * t));
  d2twodt2 = dtwodt * (-b_2 / (t * t)) + two * (2. * b_2 / (t * t * t));
  four = exp(a_4 + (b_4 / t));
  dfourdt = four * (-b_4 / (t * t));
  d2fourdt2 = dfourdt * (-b_4 / (t * t)) + four * (2. * b_4 / (t * t * t));
  u = ps * two;
  dudt = dpsdt * two + ps * dtwodt;
  d2udt2 = d2psdt2 * two + 2. * dpsdt * dtwodt + ps * d2twodt2;
  w = four * ps * ps * ps;
  dwdt = dfourdt * ps * ps * ps + 3. * four * ps * ps * dpsdt;
  d2wdt2 = d2fourdt2 * ps * ps * ps + dfourdt * 3. * ps * ps * dpsdt +
           3. * dfourdt * ps * ps * dpsdt + +3. * four * 2. * ps * dpsdt * dpsdt +
           3. * four * ps * ps * d2psdt2;
  x1 = 0.8;
  x1p = 0.;
  dx1dt = 0.;
  d2x1dt2 = 0.;
  while (fabs(x1 - x1p) / x1 > 1.e-6)
  {
    x1p = x1;
    x1_p2 = x1 * x1;
    x1_p3 = x1_p2 * x1;
    x1_p4 = x1_p3 * x1;
    num = w * x1_p4 + u * x1_p2 + x1 - 1.;
    dnumdt = dwdt * x1_p4 + w * 4. * x1_p3 * dx1dt + dudt * x1_p2 + u * 2. * x1 * dx1dt + dx1dt;
    d2numdt2 = d2wdt2 * x1_p4 + dwdt * 4. * x1_p3 * dx1dt + dwdt * 4. * x1_p3 * dx1dt +
               w * 12. * x1_p2 * dx1dt * dx1dt + w * 4. * x1_p3 * d2x1dt2 + d2udt2 * x1_p2 +
               dudt * 2. * x1 * dx1dt + dudt * 2. * x1 * dx1dt + u * 2. * dx1dt * dx1dt +
               u * 2. * x1 * d2x1dt2 + d2x1dt2;
    den = 4. * w * x1_p3 + 2. * u * x1 + 1.;
    ddendt = 4. * dwdt * x1_p3 + 4. * w * 3. * x1_p2 * dx1dt + 2. * dudt * x1 + 2. * u * dx1dt;
    d2dendt2 = 4. * d2wdt2 * x1_p3 + 4. * dwdt * 3. * x1_p2 * dx1dt +
               4. * dwdt * 3. * x1_p2 * dx1dt + 4. * w * 6. * x1 * dx1dt * dx1dt +
               4. * w * 3. * x1_p2 * d2x1dt2 + 2. * d2udt2 * x1 + 2. * dudt * dx1dt +
               2. * dudt * dx1dt + 2. * u * d2x1dt2;
    x1 = x1 - num / den;
    dx1dt = dx1dt - (dnumdt * den - num * ddendt) / (den * den);
    d2x1dt2 = d2x1dt2 -
              ((d2numdt2 * den + dnumdt * ddendt - dnumdt * ddendt - num * d2dendt2) * (den * den) -
               (dnumdt * den - num * ddendt) * 2. * den * ddendt) /
                  (den * den * den * den);
  }
  x2 = u * x1 * x1;
  dx2dt = dudt * x1 * x1 + u * 2. * x1 * dx1dt;
  d2x2dt2 = d2udt2 * x1 * x1 + dudt * 2. * x1 * dx1dt + dudt * 2. * x1 * dx1dt +
            u * 2. * dx1dt * dx1dt + u * 2. * x1 * d2x1dt2;
  x4 = 1. - x1 - x2;
  dx4dt = -dx1dt - dx2dt;
  d2x4dt2 = -d2x1dt2 - d2x2dt2;
  abar = c * (x1 + 2. * x2 + 4. * x4);
  dabardt = c * (dx1dt + 2. * dx2dt + 4. * dx4dt);
  d2abardt2 = c * (d2x1dt2 + 2. * d2x2dt2 + 4. * d2x4dt2);
  //
  bd = x1 + 2. * x2 + 4. * x4;
  dbddt = dx1dt + 2. * dx2dt + 4. * dx4dt;
  d2bddt2 = d2x1dt2 + 2. * d2x2dt2 + 4. * d2x4dt2;
  b2 = 2. * x2 / bd;
  num = dx2dt * bd - x2 * dbddt;
  dnumdt = d2x2dt2 * bd + dx2dt * dbddt - dx2dt * dbddt - x2 * d2bddt2;
  den = bd * bd;
  ddendt = 2. * bd * dbddt;
  db2dt = 2. * num / den;
  d2b2dt2 = 2. * (dnumdt * den - num * ddendt) / (den * den);
  b4 = 4. * x4 / bd;
  num = dx4dt * bd - x4 * dbddt;
  dnumdt = d2x4dt2 * bd + dx4dt * dbddt - dx4dt * dbddt - x4 * d2bddt2;
  den = bd * bd;
  ddendt = 2. * bd * dbddt;
  db4dt = 4. * num / den;
  d2b4dt2 = 4. * (dnumdt * den - num * ddendt) / (den * den);
}
//
void
DIFF_qc_tp_K(double t,
              double p,
              double & x1,
              double & dx1dt,
              double & d2x1dt2,
              double & dx1dp,
              double & d2x1dp2,
              double & d2x1dtdp,
              double & x2,
              double & dx2dt,
              double & d2x2dt2,
              double & dx2dp,
              double & d2x2dp2,
              double & d2x2dtdp,
              double & x4,
              double & dx4dt,
              double & d2x4dt2,
              double & dx4dp,
              double & d2x4dp2,
              double & d2x4dtdp,
              double & abar,
              double & dabardt,
              double & d2abardt2,
              double & dabardp,
              double & d2abardp2,
              double & d2abardtdp,
              double & bd,
              double & dbddt,
              double & d2bddt2,
              double & dbddp,
              double & d2bddp2,
              double & d2bddtdp,
              double & b2,
              double & db2dt,
              double & d2b2dt2,
              double & db2dp,
              double & d2b2dp2,
              double & d2b2dtdp,
              double & b4,
              double & db4dt,
              double & d2b4dt2,
              double & db4dp,
              double & d2b4dp2,
              double & d2b4dtdp)
{ // units: R
  //
  double two, four, u, w, x1p, dwdt, d2wdt2, dwdp, d2wdp2;
  double dtwodt, d2twodt2, dfourdt, d2fourdt2, dudt, d2udt2;
  double dudp, d2udp2;
  double d2udtdp, d2wdtdp;
  double num, dnumdt, d2numdt2, den, ddendt, d2dendt2;
  double dnumdp, d2numdp2, ddendp, d2dendp2, d2numdtdp, d2dendtdp;
  double tnum, dtnumdt, pnum, dpnumdp, dtnumdp;
  double p_p2, p_p3, x1_p2, x1_p3, x1_p4;
  //
  static const double a_2 = -8.9033;
  static const double b_2 = 12250.1;
  static const double a_4 = -23.394;
  static const double b_4 = 31694.6;
  static const double c = 39.102;
  //
  two = exp(a_2 + (b_2 / t));
  dtwodt = two * (-b_2 / (t * t));
  d2twodt2 = dtwodt * (-b_2 / (t * t)) + two * (2. * b_2 / (t * t * t));
  four = exp(a_4 + (b_4 / t));
  dfourdt = four * (-b_4 / (t * t));
  d2fourdt2 = dfourdt * (-b_4 / (t * t)) + four * (2. * b_4 / (t * t * t));
  u = p * two;
  dudt = p * dtwodt;
  d2udt2 = p * d2twodt2;
  dudp = two;
  d2udp2 = 0.;
  d2udtdp = dtwodt;
  p_p2 = p * p;
  p_p3 = p_p2 * p;
  w = four * p_p3;
  dwdt = dfourdt * p_p3;
  d2wdt2 = d2fourdt2 * p_p3;
  dwdp = four * 3. * p_p2;
  d2wdp2 = four * 6. * p;
  d2wdtdp = dfourdt * 3. * p_p2;
  x1 = 0.8;
  x1p = 0.;
  dx1dt = 0.;
  d2x1dt2 = 0.;
  dx1dp = 0.;
  d2x1dp2 = 0.;
  d2x1dtdp = 0.;
  while (fabs(x1 - x1p) / x1 > 1.e-6)
  {
    x1p = x1;
    x1_p2 = x1 * x1;
    x1_p3 = x1_p2 * x1;
    x1_p4 = x1_p3 * x1;
    num = w * x1_p4 + u * x1_p2 + x1 - 1.0;
    dnumdt = dwdt * x1_p4 + w * 4. * x1_p3 * dx1dt + dudt * x1_p2 + u * 2. * x1 * dx1dt + dx1dt;
    d2numdt2 = d2wdt2 * x1_p4 + dwdt * 4. * x1_p3 * dx1dt + dwdt * 4. * x1_p3 * dx1dt +
               w * 12. * x1_p2 * dx1dt * dx1dt + w * 4. * x1_p3 * d2x1dt2 + d2udt2 * x1_p2 +
               dudt * 2. * x1 * dx1dt + dudt * 2. * x1 * dx1dt + u * 2. * dx1dt * dx1dt +
               u * 2. * x1 * d2x1dt2 + d2x1dt2;
    dnumdp = dwdp * x1_p4 + w * 4. * x1_p3 * dx1dp + dudp * x1_p2 + u * 2. * x1 * dx1dp + dx1dp;
    d2numdp2 = d2wdp2 * x1_p4 + dwdp * 4. * x1_p3 * dx1dp + dwdp * 4. * x1_p3 * dx1dp +
               w * 12. * x1_p2 * dx1dp * dx1dp + w * 4. * x1_p3 * d2x1dp2 + d2udp2 * x1_p2 +
               dudp * 2. * x1 * dx1dp + dudp * 2. * x1 * dx1dp + u * 2. * dx1dp * dx1dp +
               u * 2. * x1 * d2x1dp2 + d2x1dp2;
    d2numdtdp = d2wdtdp * x1_p4 + dwdt * 4. * x1_p3 * dx1dp + dwdp * 4. * x1_p3 * dx1dt +
                w * 12. * x1_p2 * dx1dp * dx1dt + w * 4. * x1_p3 * d2x1dtdp + d2udtdp * x1_p2 +
                dudt * 2. * x1 * dx1dp + dudp * 2. * x1 * dx1dt + u * 2. * dx1dp * dx1dt +
                u * 2. * x1 * d2x1dtdp + d2x1dtdp;
    den = 4. * w * x1_p3 + 2. * u * x1 + 1.;
    ddendt = 4. * dwdt * x1_p3 + 4. * w * 3. * x1_p2 * dx1dt + 2. * dudt * x1 + 2. * u * dx1dt;
    d2dendt2 = 4. * d2wdt2 * x1_p3 + 4. * dwdt * 3. * x1_p2 * dx1dt +
               4. * dwdt * 3. * x1_p2 * dx1dt + 4. * w * 6. * x1 * dx1dt * dx1dt +
               4. * w * 3. * x1_p2 * d2x1dt2 + 2. * d2udt2 * x1 + 2. * dudt * dx1dt +
               2. * dudt * dx1dt + 2. * u * d2x1dt2;
    ddendp = 4. * dwdp * x1_p3 + 4. * w * 3. * x1_p2 * dx1dp + 2. * dudp * x1 + 2. * u * dx1dp;
    d2dendp2 = 4. * d2wdp2 * x1_p3 + 12. * dwdp * x1_p2 * dx1dp + 4. * dwdp * 3. * x1_p2 * dx1dp +
               4. * w * 6. * x1 * dx1dp * dx1dp + 4. * w * 3. * x1_p2 * d2x1dp2 + 2. * d2udp2 * x1 +
               2. * dudp * dx1dp + 2. * dudp * dx1dp + 2. * u * d2x1dp2;
    d2dendtdp = 4. * d2wdtdp * x1_p3 + 12. * dwdp * x1_p2 * dx1dt + 4. * dwdt * 3. * x1_p2 * dx1dp +
                4. * w * 6. * x1 * dx1dt * dx1dp + 4. * w * 3. * x1_p2 * d2x1dtdp +
                2. * d2udtdp * x1 + 2. * dudp * dx1dt + 2. * dudt * dx1dp + 2. * u * d2x1dtdp;
    x1 = x1 - num / den;
    dx1dt = dx1dt - (dnumdt * den - num * ddendt) / (den * den);
    d2x1dt2 = d2x1dt2 -
              ((d2numdt2 * den + dnumdt * ddendt - dnumdt * ddendt - num * d2dendt2) * (den * den) -
               (dnumdt * den - num * ddendt) * 2. * den * ddendt) /
                  (den * den * den * den);
    dx1dp = dx1dp - (dnumdp * den - num * ddendp) / (den * den);
    d2x1dp2 = d2x1dp2 -
              ((d2numdp2 * den + dnumdp * ddendp - dnumdp * ddendp - num * d2dendp2) * (den * den) -
               (dnumdp * den - num * ddendp) * 2. * den * ddendp) /
                  (den * den * den * den);
    d2x1dtdp = d2x1dtdp - ((d2numdtdp * den + dnumdt * ddendp - dnumdp * ddendt - num * d2dendtdp) *
                               (den * den) -
                           (dnumdt * den - num * ddendt) * 2. * den * ddendp) /
                              (den * den * den * den);
  }
  x2 = u * x1 * x1;
  dx2dt = dudt * x1 * x1 + u * 2. * x1 * dx1dt;
  d2x2dt2 = d2udt2 * x1 * x1 + dudt * 2. * x1 * dx1dt + dudt * 2. * x1 * dx1dt +
            u * 2. * dx1dt * dx1dt + u * 2. * x1 * d2x1dt2;
  dx2dp = dudp * x1 * x1 + u * 2. * x1 * dx1dp;
  d2x2dp2 = d2udp2 * x1 * x1 + dudp * 2. * x1 * dx1dp + dudp * 2. * x1 * dx1dp +
            u * 2. * dx1dp * dx1dp + u * 2. * x1 * d2x1dp2;
  d2x2dtdp = d2udtdp * x1 * x1 + dudt * 2. * x1 * dx1dp + dudp * 2. * x1 * dx1dt +
             u * 2. * dx1dp * dx1dt + u * 2. * x1 * d2x1dtdp;
  x4 = 1. - x1 - x2;
  dx4dt = -dx1dt - dx2dt;
  d2x4dt2 = -d2x1dt2 - d2x2dt2;
  dx4dp = -dx1dp - dx2dp;
  d2x4dp2 = -d2x1dp2 - d2x2dp2;
  d2x4dtdp = -d2x1dtdp - d2x2dtdp;
  abar = c * (x1 + 2. * x2 + 4. * x4);
  dabardt = c * (dx1dt + 2. * dx2dt + 4. * dx4dt);
  d2abardt2 = c * (d2x1dt2 + 2. * d2x2dt2 + 4. * d2x4dt2);
  dabardp = c * (dx1dp + 2. * dx2dp + 4. * dx4dp);
  d2abardp2 = c * (d2x1dp2 + 2. * d2x2dp2 + 4. * d2x4dp2);
  d2abardtdp = c * (d2x1dtdp + 2. * d2x2dtdp + 4. * d2x4dtdp);
  //
  bd = x1 + 2. * x2 + 4. * x4;
  dbddt = dx1dt + 2. * dx2dt + 4. * dx4dt;
  d2bddt2 = d2x1dt2 + 2. * d2x2dt2 + 4. * d2x4dt2;
  dbddp = dx1dp + 2. * dx2dp + 4. * dx4dp;
  d2bddp2 = d2x1dp2 + 2. * d2x2dp2 + 4. * d2x4dp2;
  d2bddtdp = d2x1dtdp + 2. * d2x2dtdp + 4. * d2x4dtdp;
  b2 = 2. * x2 / bd;
  tnum = dx2dt * bd - x2 * dbddt;
  dtnumdt = d2x2dt2 * bd + dx2dt * dbddt - dx2dt * dbddt - x2 * d2bddt2;
  den = bd * bd;
  ddendt = 2. * bd * dbddt;
  db2dt = 2. * tnum / den;
  d2b2dt2 = 2. * (dtnumdt * den - tnum * ddendt) / (den * den);
  pnum = dx2dp * bd - x2 * dbddp;
  dpnumdp = d2x2dp2 * bd + dx2dp * dbddp - dx2dp * dbddp - x2 * d2bddp2;
  ddendp = 2. * bd * dbddp;
  db2dp = 2. * pnum / den;
  d2b2dp2 = 2. * (dpnumdp * den - pnum * ddendp) / (den * den);
  dtnumdp = d2x2dtdp * bd + dx2dt * dbddp - dx2dp * dbddt - x2 * d2bddtdp;
  d2b2dtdp = 2. * (dtnumdp * den - tnum * ddendp) / (den * den);
  b4 = 4. * x4 / bd;
  tnum = dx4dt * bd - x4 * dbddt;
  dtnumdt = d2x4dt2 * bd + dx4dt * dbddt - dx4dt * dbddt - x4 * d2bddt2;
  den = bd * bd;
  ddendt = 2. * bd * dbddt;
  db4dt = 4. * tnum / den;
  d2b4dt2 = 4. * (dtnumdt * den - tnum * ddendt) / (den * den);
  pnum = dx4dp * bd - x4 * dbddp;
  dpnumdp = d2x4dp2 * bd + dx4dp * dbddp - dx4dp * dbddp - x4 * d2bddp2;
  ddendp = 2. * bd * dbddp;
  db4dp = 4. * pnum / den;
  d2b4dp2 = 4. * (dpnumdp * den - pnum * ddendp) / (den * den);
  dtnumdp = d2x4dtdp * bd + dx4dt * dbddp - dx4dp * dbddt - x4 * d2bddtdp;
  d2b4dtdp = 4. * (dtnumdp * den - tnum * ddendp) / (den * den);
}
//
void
DIFF_v2_t_K(double t, double & v2, double & dv2dt, double & d2v2dt2)
{ // units: R, ft3/lb
  double ps, dpsdt, d2psdt2, x1, dx1dt, d2x1dt2, x2, dx2dt, d2x2dt2, x4, dx4dt, d2x4dt2, abar,
      dabardt, d2abardt2, bd, dbddt, d2bddt2, b2, db2dt, d2b2dt2, b4, db4dt, d2b4dt2;
  double num, dnumdt, den, ddendt;
  //
  static const double R = 0.730229;
  //
  DIFF_ps_t_K(t, ps, dpsdt, d2psdt2);
  DIFF_qc_t_K(t,
               ps,
               dpsdt,
               d2psdt2,
               x1,
               dx1dt,
               d2x1dt2,
               x2,
               dx2dt,
               d2x2dt2,
               x4,
               dx4dt,
               d2x4dt2,
               abar,
               dabardt,
               d2abardt2,
               bd,
               dbddt,
               d2bddt2,
               b2,
               db2dt,
               d2b2dt2,
               b4,
               db4dt,
               d2b4dt2);
  //
  v2 = R * t / (abar * ps);
  num = abar * ps - t * (dabardt * ps + abar * dpsdt);
  dnumdt = dabardt * ps + abar * dpsdt - (dabardt * ps + abar * dpsdt) -
           t * (d2abardt2 * ps + dabardt * dpsdt + dabardt * dpsdt + abar * d2psdt2);
  den = (abar * ps) * (abar * ps);
  ddendt = 2. * (abar * ps) * (dabardt * ps + abar * dpsdt);
  dv2dt = R * num / den;
  d2v2dt2 = R * (dnumdt * den - num * ddendt) / (den * den);
}
//
void
DIFF_v_tp_G_K(double t,
               double p,
               double & vv,
               double & dvvdt,
               double & d2vvdt2,
               double & dvvdp,
               double & d2vvdp2,
               double & d2vvdtdp)
{ // units: R, ft3/lb
  double x1h, dx1hdt, d2x1hdt2, dx1hdp, d2x1hdp2, d2x1hdtdp, x2h, dx2hdt, d2x2hdt2, dx2hdp,
      d2x2hdp2, d2x2hdtdp, x4h, dx4hdt, d2x4hdt2, dx4hdp, d2x4hdp2, d2x4hdtdp, abarh, dabarhdt,
      d2abarhdt2, dabarhdp, d2abarhdp2, d2abarhdtdp, bdh, dbdhdt, d2bdhdt2, dbdhdp, d2bdhdp2,
      d2bdhdtdp, b2h, db2hdt, d2b2hdt2, db2hdp, d2b2hdp2, d2b2hdtdp, b4h, db4hdt, d2b4hdt2, db4hdp,
      d2b4hdp2, d2b4hdtdp;
  double den, ddendt, d2dendt2, ddendp, d2dendp2, d2dendtdp, num_dt, dnum_dtdt, den1, dden1dt;
  //
  static const double R = 0.730229;
  //
  DIFF_qc_tp_K(t,
                p,
                x1h,
                dx1hdt,
                d2x1hdt2,
                dx1hdp,
                d2x1hdp2,
                d2x1hdtdp,
                x2h,
                dx2hdt,
                d2x2hdt2,
                dx2hdp,
                d2x2hdp2,
                d2x2hdtdp,
                x4h,
                dx4hdt,
                d2x4hdt2,
                dx4hdp,
                d2x4hdp2,
                d2x4hdtdp,
                abarh,
                dabarhdt,
                d2abarhdt2,
                dabarhdp,
                d2abarhdp2,
                d2abarhdtdp,
                bdh,
                dbdhdt,
                d2bdhdt2,
                dbdhdp,
                d2bdhdp2,
                d2bdhdtdp,
                b2h,
                db2hdt,
                d2b2hdt2,
                db2hdp,
                d2b2hdp2,
                d2b2hdtdp,
                b4h,
                db4hdt,
                d2b4hdt2,
                db4hdp,
                d2b4hdp2,
                d2b4hdtdp);
  //
  den = abarh * p;
  ddendt = dabarhdt * p;
  d2dendt2 = d2abarhdt2 * p;
  ddendp = dabarhdp * p + abarh;
  d2dendp2 = d2abarhdp2 * p + 2. * dabarhdp;
  d2dendtdp = d2abarhdtdp * p + dabarhdt;
  vv = R * t / den;
  num_dt = den - t * ddendt;
  dnum_dtdt = -t * d2dendt2;
  den1 = den * den;
  dden1dt = 2. * den * ddendt;
  dvvdt = R * num_dt / den1;
  d2vvdt2 = R * (dnum_dtdt * den1 - num_dt * dden1dt) / (den1 * den1);
  dvvdp = -R * t / (den * den) * ddendp;
  d2vvdp2 = R * t * (2. / (den * den * den) * ddendp * ddendp - d2dendp2 / (den * den));
  d2vvdtdp = -R * ((den1 - t * dden1dt) / (den1 * den1) * ddendp + t / (den * den) * d2dendtdp);
}
//
void
DIFF_dhv_t_K(double t, double & dhv, double & ddhvdt, double & d2dhvdt2)
{ // units: R, Btu/lb
  double ps, dpsdt, d2psdt2, x1, dx1dt, d2x1dt2, x2, dx2dt, d2x2dt2, x4, dx4dt, d2x4dt2, abar,
      dabardt, d2abardt2, bd, dbddt, d2bddt2, b2, db2dt, d2b2dt2, b4, db4dt, d2b4dt2;
  double dhf1, ddhf1dt, d2dhf1dt2, dhf2, ddhf2dt, d2dhf2dt2;
  double dhf4, ddhf4dt, d2dhf4dt2, num, dnumdt, d2numdt2;
  //
  static const double a = 21856.5;
  static const double b = -2.1734;
  static const double c = 7.0470e-4;
  static const double d = -1.6816e-7;
  static const double e2 = -13500.;
  static const double e4 = -34920.;
  //
  DIFF_ps_t_K(t, ps, dpsdt, d2psdt2);
  DIFF_qc_t_K(t,
               ps,
               dpsdt,
               d2psdt2,
               x1,
               dx1dt,
               d2x1dt2,
               x2,
               dx2dt,
               d2x2dt2,
               x4,
               dx4dt,
               d2x4dt2,
               abar,
               dabardt,
               d2abardt2,
               bd,
               dbddt,
               d2bddt2,
               b2,
               db2dt,
               d2b2dt2,
               b4,
               db4dt,
               d2b4dt2);
  //
  dhf1 = a + t * (b + t * (c + t * d));
  ddhf1dt = b + t * (2. * c + t * 3. * d);
  d2dhf1dt2 = 2. * c + t * 6. * d;
  dhf2 = 2. * dhf1 + e2;
  ddhf2dt = 2. * ddhf1dt;
  d2dhf2dt2 = 2. * d2dhf1dt2;
  dhf4 = 4. * dhf1 + e4;
  ddhf4dt = 4. * ddhf1dt;
  d2dhf4dt2 = 4. * d2dhf1dt2;
  num = x1 * dhf1 + x2 * dhf2 + x4 * dhf4;
  dnumdt = dx1dt * dhf1 + x1 * ddhf1dt + dx2dt * dhf2 + x2 * ddhf2dt + dx4dt * dhf4 + x4 * ddhf4dt;
  d2numdt2 = d2x1dt2 * dhf1 + dx1dt * ddhf1dt + dx1dt * ddhf1dt + x1 * d2dhf1dt2 + d2x2dt2 * dhf2 +
             dx2dt * ddhf2dt + dx2dt * ddhf2dt + x2 * d2dhf2dt2 + d2x4dt2 * dhf4 + dx4dt * ddhf4dt +
             dx4dt * ddhf4dt + x4 * d2dhf4dt2;
  dhv = 1.8 * num / abar;
  ddhvdt = 1.8 * (dnumdt * abar - num * dabardt) / (abar * abar);
  d2dhvdt2 =
      1.8 *
      ((d2numdt2 * abar + dnumdt * dabardt - dnumdt * dabardt - num * d2abardt2) * (abar * abar) -
       (dnumdt * abar - num * dabardt) * 2. * abar * dabardt) /
      (abar * abar * abar * abar);
}
//
void
DIFF_h2_t_K(double t, double & h2, double & dh2dt, double & d2h2dt2)
{ // units: R, Btu/lb
  double h1, dh1dt, d2h1dt2, dhv, ddhvdt, d2dhvdt2;
  //
  DIFF_h1_t_K(t, h1, dh1dt, d2h1dt2);
  DIFF_dhv_t_K(t, dhv, ddhvdt, d2dhvdt2);
  //
  h2 = h1 + dhv;
  dh2dt = dh1dt + ddhvdt;
  d2h2dt2 = d2h1dt2 + d2dhvdt2;
}
//
void
DIFF_vu_tp_G_K(double t,
                double p,
                double & vv,
                double & dvvdt,
                double & d2vvdt2,
                double & dvvdp,
                double & d2vvdp2,
                double & d2vvdtdp,
                double & uv,
                double & duvdt,
                double & d2uvdt2,
                double & duvdp,
                double & d2uvdp2,
                double & d2uvdtdp)
{ // units: R, atm, Btu/lb

  static const double atmft3_toBtu = (101325. * 0.3048 * 0.3048 * 0.3048) / 1055.05585262;

  double h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  DIFF_v_tp_G_K(t, p, vv, dvvdt, d2vvdt2, dvvdp, d2vvdp2, d2vvdtdp);
  DIFF_h_tp_G_K(t, p, h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);

  uv = h - p * vv * atmft3_toBtu;
  duvdt = dhdt - p * dvvdt * atmft3_toBtu;
  d2uvdt2 = d2hdt2 - p * d2vvdt2 * atmft3_toBtu;
  duvdp = dhdp - (vv + p * dvvdp) * atmft3_toBtu;
  d2uvdp2 = d2hdp2 - (dvvdp + dvvdp + p * d2vvdp2) * atmft3_toBtu;
  d2uvdtdp = d2hdtdp - (dvvdt + p * d2vvdtdp) * atmft3_toBtu;
}
//
void
DIFF_h_tp_G_K(double t,
               double p,
               double & hv,
               double & dhvdt,
               double & d2hvdt2,
               double & dhvdp,
               double & d2hvdp2,
               double & d2hvdtdp)
{ // units: R, Btu/lb
  double h1, dh1dt, d2h1dt2, dhv, ddhvdt, d2dhvdt2;
  double h2, dh2dt, d2h2dt2;
  double ps, dpsdt, d2psdt2;
  double x1, dx1dt, d2x1dt2, x2, dx2dt, d2x2dt2, x4, dx4dt, d2x4dt2, abar, dabardt, d2abardt2, bd,
      dbddt, d2bddt2, b2, db2dt, d2b2dt2, b4, db4dt, d2b4dt2;
  double x1h, dx1hdt, d2x1hdt2, dx1hdp, d2x1hdp2, d2x1hdtdp, x2h, dx2hdt, d2x2hdt2, dx2hdp,
      d2x2hdp2, d2x2hdtdp, x4h, dx4hdt, d2x4hdt2, dx4hdp, d2x4hdp2, d2x4hdtdp, abarh, dabarhdt,
      d2abarhdt2, dabarhdp, d2abarhdp2, d2abarhdtdp, bdh, dbdhdt, d2bdhdt2, dbdhdp, d2bdhdp2,
      d2bdhdtdp, b2h, db2hdt, d2b2hdt2, db2hdp, d2b2hdp2, d2b2hdtdp, b4h, db4hdt, d2b4hdt2, db4hdp,
      d2b4hdp2, d2b4hdtdp;
  double hvi, dhvidt, d2hvidt2;
  //
  static const double chv1 = 310.67;
  static const double chv2 = 401.89;
  //
  DIFF_h1_t_K(t, h1, dh1dt, d2h1dt2);
  DIFF_dhv_t_K(t, dhv, ddhvdt, d2dhvdt2);
  //
  h2 = h1 + dhv;
  dh2dt = dh1dt + ddhvdt;
  d2h2dt2 = d2h1dt2 + d2dhvdt2;
  //
  DIFF_ps_t_K(t, ps, dpsdt, d2psdt2);
  DIFF_qc_t_K(t,
               ps,
               dpsdt,
               d2psdt2,
               x1,
               dx1dt,
               d2x1dt2,
               x2,
               dx2dt,
               d2x2dt2,
               x4,
               dx4dt,
               d2x4dt2,
               abar,
               dabardt,
               d2abardt2,
               bd,
               dbddt,
               d2bddt2,
               b2,
               db2dt,
               d2b2dt2,
               b4,
               db4dt,
               d2b4dt2);
  //
  DIFF_qc_tp_K(t,
                p,
                x1h,
                dx1hdt,
                d2x1hdt2,
                dx1hdp,
                d2x1hdp2,
                d2x1hdtdp,
                x2h,
                dx2hdt,
                d2x2hdt2,
                dx2hdp,
                d2x2hdp2,
                d2x2hdtdp,
                x4h,
                dx4hdt,
                d2x4hdt2,
                dx4hdp,
                d2x4hdp2,
                d2x4hdtdp,
                abarh,
                dabarhdt,
                d2abarhdt2,
                dabarhdp,
                d2abarhdp2,
                d2abarhdtdp,
                bdh,
                dbdhdt,
                d2bdhdt2,
                dbdhdp,
                d2bdhdp2,
                d2bdhdtdp,
                b2h,
                db2hdt,
                d2b2hdt2,
                db2hdp,
                d2b2hdp2,
                d2b2hdtdp,
                b4h,
                db4hdt,
                d2b4hdt2,
                db4hdp,
                d2b4hdp2,
                d2b4hdtdp);
  //
  hvi = h2 + b2 * chv1 + b4 * chv2;
  dhvidt = dh2dt + db2dt * chv1 + db4dt * chv2;
  d2hvidt2 = d2h2dt2 + d2b2dt2 * chv1 + d2b4dt2 * chv2;
  hv = hvi - b2h * chv1 - b4h * chv2;
  dhvdt = dhvidt - db2hdt * chv1 - db4hdt * chv2;
  d2hvdt2 = d2hvidt2 - d2b2hdt2 * chv1 - d2b4hdt2 * chv2;
  dhvdp = -db2hdp * chv1 - db4hdp * chv2;
  d2hvdp2 = -d2b2hdp2 * chv1 - d2b4hdp2 * chv2;
  d2hvdtdp = -d2b2hdtdp * chv1 - d2b4hdtdp * chv2;
}
//
void
DIFF_s2_t_K(double t, double & s2, double & ds2dt, double & d2s2dt2)
{ // units: R, Btu/(lb*F)
  double s1, ds1dt, d2s1dt2, dhv, ddhvdt, d2dhvdt2;
  //
  DIFF_s1_t_K(t, s1, ds1dt, d2s1dt2);
  DIFF_dhv_t_K(t, dhv, ddhvdt, d2dhvdt2);
  //
  s2 = s1 + (dhv / t);
  ds2dt = ds1dt + (ddhvdt * t - dhv) / (t * t);
  d2s2dt2 = d2s1dt2 + (d2dhvdt2 * t * t - (ddhvdt * t - dhv) * 2.0) / (t * t * t);
}
//
void
DIFF_s_tp_G_K(double t,
               double p,
               double & sv,
               double & dsvdt,
               double & d2svdt2,
               double & dsvdp,
               double & d2svdp2,
               double & d2svdtdp)
{ // units: R, Btu/(lb*F)
  double s1, ds1dt, d2s1dt2, dhv, ddhvdt, d2dhvdt2;
  double s2, ds2dt, d2s2dt2;
  double ps, dpsdt, d2psdt2;
  double x1, dx1dt, d2x1dt2, x2, dx2dt, d2x2dt2, x4, dx4dt, d2x4dt2, abar, dabardt, d2abardt2, bd,
      dbddt, d2bddt2, b2, db2dt, d2b2dt2, b4, db4dt, d2b4dt2;
  double x1h, dx1hdt, d2x1hdt2, dx1hdp, d2x1hdp2, d2x1hdtdp, x2h, dx2hdt, d2x2hdt2, dx2hdp,
      d2x2hdp2, d2x2hdtdp, x4h, dx4hdt, d2x4hdt2, dx4hdp, d2x4hdp2, d2x4hdtdp, abarh, dabarhdt,
      d2abarhdt2, dabarhdp, d2abarhdp2, d2abarhdtdp, bdh, dbdhdt, d2bdhdt2, dbdhdp, d2bdhdp2,
      d2bdhdtdp, b2h, db2hdt, d2b2hdt2, db2hdp, d2b2hdp2, d2b2hdtdp, b4h, db4hdt, d2b4hdt2, db4hdp,
      d2b4hdp2, d2b4hdtdp;
  //
  double term1, dterm1dt, d2term1dt2, dterm1dp, d2term1dp2, d2term1dtdp;
  double term2, dterm2dt, d2term2dt2, dterm2dp, d2term2dp2, d2term2dtdp;
  double term3, dterm3dt, d2term3dt2, dterm3dp, d2term3dp2, d2term3dtdp;
  double term41, dterm41dt, d2term41dt2, dterm41dp, d2term41dp2, d2term41dtdp;
  double term42, dterm42dt, d2term42dt2;
  double term4, dterm4dt, d2term4dt2, dterm4dp, d2term4dp2, d2term4dtdp;
  double term51, dterm51dt, d2term51dt2, dterm51dp, d2term51dp2, d2term51dtdp;
  double term52, dterm52dt, d2term52dt2;
  double term5, dterm5dt, d2term5dt2, dterm5dp, d2term5dp2, d2term5dtdp;
  double term61, dterm61dt, d2term61dt2, dterm61dp, d2term61dp2, d2term61dtdp;
  double term62, dterm62dt, dterm62dp, d2term62dt2, d2term62dp2, d2term62dtdp;
  double term6, dterm6dt, d2term6dt2, dterm6dp, d2term6dp2, d2term6dtdp;
  //
  double svi, dsvidt, d2svidt2;
  //
  static const double chv1 = 310.67;
  static const double chv2 = 401.89;
  static const double csv1 = 1.987180;
  static const double csv2 = 78.204;
  static const double csv3 = -17.66723;
  static const double csv4 = 24308.0;
  static const double csv5 = 156.4;
  static const double csv6 = -46.422;
  static const double csv7 = 62893.0;
  static const double csv8 = 1.987180;
  //
  DIFF_s1_t_K(t, s1, ds1dt, d2s1dt2);
  DIFF_dhv_t_K(t, dhv, ddhvdt, d2dhvdt2);
  //
  s2 = s1 + (dhv / t);
  ds2dt = ds1dt + (ddhvdt * t - dhv) / (t * t);
  d2s2dt2 = d2s1dt2 + (d2dhvdt2 * t * t - (ddhvdt * t - dhv) * 2.0) / (t * t * t);
  //
  DIFF_ps_t_K(t, ps, dpsdt, d2psdt2);
  DIFF_qc_t_K(t,
               ps,
               dpsdt,
               d2psdt2,
               x1,
               dx1dt,
               d2x1dt2,
               x2,
               dx2dt,
               d2x2dt2,
               x4,
               dx4dt,
               d2x4dt2,
               abar,
               dabardt,
               d2abardt2,
               bd,
               dbddt,
               d2bddt2,
               b2,
               db2dt,
               d2b2dt2,
               b4,
               db4dt,
               d2b4dt2);
  //
  DIFF_qc_tp_K(t,
                p,
                x1h,
                dx1hdt,
                d2x1hdt2,
                dx1hdp,
                d2x1hdp2,
                d2x1hdtdp,
                x2h,
                dx2hdt,
                d2x2hdt2,
                dx2hdp,
                d2x2hdp2,
                d2x2hdtdp,
                x4h,
                dx4hdt,
                d2x4hdt2,
                dx4hdp,
                d2x4hdp2,
                d2x4hdtdp,
                abarh,
                dabarhdt,
                d2abarhdt2,
                dabarhdp,
                d2abarhdp2,
                d2abarhdtdp,
                bdh,
                dbdhdt,
                d2bdhdt2,
                dbdhdp,
                d2bdhdp2,
                d2bdhdtdp,
                b2h,
                db2hdt,
                d2b2hdt2,
                db2hdp,
                d2b2hdp2,
                d2b2hdtdp,
                b4h,
                db4hdt,
                d2b4hdt2,
                db4hdp,
                d2b4hdp2,
                d2b4hdtdp);
  //
  double log_ps = log(ps);
  double log_p = log(p);
  double log_x1 = log(x1);
  double log_x2 = log(x2);
  double log_x4 = log(x4);
  double log_x1h = log(x1h);
  double log_x2h = log(x2h);
  double log_x4h = log(x4h);
  //
  term1 = b2 * chv1 / t;
  dterm1dt = chv1 * (db2dt * t - b2) / (t * t);
  d2term1dt2 = chv1 * (d2b2dt2 * t * t * t - (db2dt * t - b2) * 2. * t) / (t * t * t * t);
  //
  term2 = b4 * chv2 / t;
  dterm2dt = chv2 * (db4dt * t - b4) / (t * t);
  d2term2dt2 = chv2 * (d2b4dt2 * t * t * t - (db4dt * t - b4) * 2. * t) / (t * t * t * t);
  //
  term3 = (csv1 / abar) * log_ps;
  dterm3dt = (-csv1 / (abar * abar) * dabardt) * log_ps + (csv1 / abar) / ps * dpsdt;
  d2term3dt2 =
      (+csv1 * 2. / (abar * abar * abar) * dabardt * dabardt - csv1 / (abar * abar) * d2abardt2) *
          log_ps +
      (-csv1 / (abar * abar) * dabardt) / ps * dpsdt +
      ((-csv1 / (abar * abar) * dabardt) * ps - (csv1 / abar) * dpsdt) / (ps * ps) * dpsdt +
      (csv1 / abar) / ps * d2psdt2;
  //
  term41 = b2 / csv2;
  dterm41dt = db2dt / csv2;
  d2term41dt2 = d2b2dt2 / csv2;
  term42 = csv3 + csv4 / t;
  dterm42dt = -csv4 / (t * t);
  d2term42dt2 = 2. * csv4 / (t * t * t);
  term4 = term41 * term42;
  dterm4dt = dterm41dt * term42 + term41 * dterm42dt;
  d2term4dt2 =
      d2term41dt2 * term42 + dterm41dt * dterm42dt + dterm41dt * dterm42dt + term41 * d2term42dt2;
  //
  term51 = b4 / csv5;
  dterm51dt = db4dt / csv5;
  d2term51dt2 = d2b4dt2 / csv5;
  term52 = csv6 + csv7 / t;
  dterm52dt = -csv7 / (t * t);
  d2term52dt2 = 2. * csv7 / (t * t * t);
  term5 = term51 * term52;
  dterm5dt = dterm51dt * term52 + term51 * dterm52dt;
  d2term5dt2 =
      d2term51dt2 * term52 + dterm51dt * dterm52dt + dterm51dt * dterm52dt + term51 * d2term52dt2;
  //
  term61 = csv8 / abar;
  dterm61dt = -csv8 / (abar * abar) * dabardt;
  d2term61dt2 =
      +csv8 * 2. / (abar * abar * abar) * dabardt * dabardt - csv8 / (abar * abar) * d2abardt2;
  term62 = x1 * log_x1 + x2 * log_x2 + x4 * log_x4;
  dterm62dt = dx1dt * (log_x1 + 1.) + dx2dt * (log_x2 + 1.) + dx4dt * (log_x4 + 1.);
  d2term62dt2 = d2x1dt2 * (log_x1 + 1.) + dx1dt / x1 * dx1dt + d2x2dt2 * (log_x2 + 1.) +
                dx2dt / x2 * dx2dt + d2x4dt2 * (log_x4 + 1.) + dx4dt / x4 * dx4dt;
  term6 = term61 * term62;
  dterm6dt = dterm61dt * term62 + term61 * dterm62dt;
  d2term6dt2 =
      d2term61dt2 * term62 + dterm61dt * dterm62dt + dterm61dt * dterm62dt + term61 * d2term62dt2;
  //
  svi = s2 + term1 + term2 + term3 - term4 - term5 + term6;
  dsvidt = ds2dt + dterm1dt + dterm2dt + dterm3dt - dterm4dt - dterm5dt + dterm6dt;
  d2svidt2 = d2s2dt2 + d2term1dt2 + d2term2dt2 + d2term3dt2 - d2term4dt2 - d2term5dt2 + d2term6dt2;
  //
  term1 = b2h * chv1 / t;
  dterm1dt = chv1 * (db2hdt * t - b2h) / (t * t);
  d2term1dt2 = chv1 * (d2b2hdt2 * t * t * t - (db2hdt * t - b2h) * 2. * t) / (t * t * t * t);
  dterm1dp = db2hdp * chv1 / t;
  d2term1dp2 = d2b2hdp2 * chv1 / t;
  d2term1dtdp = chv1 * (d2b2hdtdp * t - db2hdp) / (t * t);
  //
  term2 = b4h * chv2 / t;
  dterm2dt = chv2 * (db4hdt * t - b4h) / (t * t);
  d2term2dt2 = chv2 * (d2b4hdt2 * t * t * t - (db4hdt * t - b4h) * 2. * t) / (t * t * t * t);
  dterm2dp = db4hdp * chv2 / t;
  d2term2dp2 = d2b4hdp2 * chv2 / t;
  d2term2dtdp = chv2 * (d2b4hdtdp * t - db4hdp) / (t * t),
  //
      term3 = (csv1 / abarh) * log_p;
  dterm3dt = (-csv1 / (abarh * abarh) * dabarhdt) * log_p;
  d2term3dt2 = (+csv1 * 2. / (abarh * abarh * abarh) * dabarhdt * dabarhdt -
                csv1 / (abarh * abarh) * d2abarhdt2) *
               log_p;
  dterm3dp = (-csv1 / (abarh * abarh) * dabarhdp) * log_p + (csv1 / abarh) / p;
  d2term3dp2 = (csv1 * 2. / (abarh * abarh * abarh) * dabarhdp * dabarhdp -
                csv1 / (abarh * abarh) * d2abarhdp2) *
                   log_p +
               (-csv1 / (abarh * abarh) * dabarhdp) / p +
               ((-csv1 / (abarh * abarh) * dabarhdp) * p - (csv1 / abarh)) / (p * p);
  d2term3dtdp = (csv1 * 2. / (abarh * abarh * abarh) * dabarhdp * dabarhdt -
                 csv1 / (abarh * abarh) * d2abarhdtdp) *
                    log_p +
                (-csv1 / (abarh * abarh) * dabarhdt) / p;
  //
  term41 = b2h / csv2;
  dterm41dt = db2hdt / csv2;
  d2term41dt2 = d2b2hdt2 / csv2;
  dterm41dp = db2hdp / csv2;
  d2term41dp2 = d2b2hdp2 / csv2;
  d2term41dtdp = d2b2hdtdp / csv2;
  term42 = csv3 + csv4 / t;
  dterm42dt = -csv4 / (t * t);
  d2term42dt2 = 2. * csv4 / (t * t * t);
  term4 = term41 * term42;
  dterm4dt = dterm41dt * term42 + term41 * dterm42dt;
  d2term4dt2 =
      d2term41dt2 * term42 + dterm41dt * dterm42dt + dterm41dt * dterm42dt + term41 * d2term42dt2;
  dterm4dp = dterm41dp * term42;
  d2term4dp2 = d2term41dp2 * term42;
  d2term4dtdp = d2term41dtdp * term42 + dterm41dp * dterm42dt;
  //
  term51 = b4h / csv5;
  dterm51dt = db4hdt / csv5;
  d2term51dt2 = d2b4hdt2 / csv5;
  dterm51dp = db4hdp / csv5;
  d2term51dp2 = d2b4hdp2 / csv5;
  d2term51dtdp = d2b4hdtdp / csv5;
  term52 = csv6 + csv7 / t;
  dterm52dt = -csv7 / (t * t);
  d2term52dt2 = 2. * csv7 / (t * t * t);
  term5 = term51 * term52;
  dterm5dt = dterm51dt * term52 + term51 * dterm52dt;
  d2term5dt2 =
      d2term51dt2 * term52 + dterm51dt * dterm52dt + dterm51dt * dterm52dt + term51 * d2term52dt2;
  dterm5dp = dterm51dp * term52;
  d2term5dp2 = d2term51dp2 * term52;
  d2term5dtdp = d2term51dtdp * term52 + dterm51dp * dterm52dt;
  //
  term61 = csv8 / abarh;
  dterm61dt = -csv8 / (abarh * abarh) * dabarhdt;
  d2term61dt2 = +csv8 * 2. / (abarh * abarh * abarh) * dabarhdt * dabarhdt -
                csv8 / (abarh * abarh) * d2abarhdt2;
  dterm61dp = -csv8 / (abarh * abarh) * dabarhdp;
  d2term61dp2 = +csv8 * 2. / (abarh * abarh * abarh) * dabarhdp * dabarhdp -
                csv8 / (abarh * abarh) * d2abarhdp2;
  d2term61dtdp = +csv8 * 2. / (abarh * abarh * abarh) * dabarhdt * dabarhdp -
                 csv8 / (abarh * abarh) * d2abarhdtdp;

  term62 = x1h * log_x1h + x2h * log_x2h + x4h * log_x4h;
  dterm62dt = dx1hdt * (log_x1h + 1.) + dx2hdt * (log_x2h + 1.) + dx4hdt * (log_x4h + 1.);
  dterm62dp = dx1hdp * (log_x1h + 1.) + dx2hdp * (log_x2h + 1.) + dx4hdp * (log_x4h + 1.);
  d2term62dt2 = d2x1hdt2 * (log_x1h + 1.) + dx1hdt / x1h * dx1hdt + d2x2hdt2 * (log_x2h + 1.) +
                dx2hdt / x2h * dx2hdt + d2x4hdt2 * (log_x4h + 1.) + dx4hdt / x4h * dx4hdt;
  d2term62dp2 = d2x1hdp2 * (log_x1h + 1.) + dx1hdp / x1h * dx1hdp + d2x2hdp2 * (log_x2h + 1.) +
                dx2hdp / x2h * dx2hdp + d2x4hdp2 * (log_x4h + 1.) + dx4hdp / x4h * dx4hdp;
  d2term62dtdp = d2x1hdtdp * (log_x1h + 1.) + dx1hdt * dx1hdp / x1h + d2x2hdtdp * (log_x2h + 1.) +
                 dx2hdt * dx2hdp / x2h + d2x4hdtdp * (log_x4h + 1.) + dx4hdt * dx4hdp / x4h;
  //
  term6 = term61 * term62;
  dterm6dt = dterm61dt * term62 + term61 * dterm62dt;
  d2term6dt2 = d2term61dt2 * term62 + 2. * dterm61dt * dterm62dt + term61 * d2term62dt2;
  dterm6dp = dterm61dp * term62 + term61 * dterm62dp;
  d2term6dp2 = d2term61dp2 * term62 + 2. * dterm61dp * dterm62dp + term61 * d2term62dp2;
  d2term6dtdp =
      d2term61dtdp * term62 + dterm61dt * dterm62dp + dterm61dp * dterm62dt + term61 * d2term62dtdp;
  //
  sv = svi - term1 - term2 - term3 + term4 + term5 - term6;
  dsvdt = dsvidt - dterm1dt - dterm2dt - dterm3dt + dterm4dt + dterm5dt - dterm6dt;
  d2svdt2 = d2svidt2 - d2term1dt2 - d2term2dt2 - d2term3dt2 + d2term4dt2 + d2term5dt2 - d2term6dt2;
  //
  dsvdp = -dterm1dp - dterm2dp - dterm3dp + dterm4dp + dterm5dp - dterm6dp;
  d2svdp2 = -d2term1dp2 - d2term2dp2 - d2term3dp2 + d2term4dp2 + d2term5dp2 - d2term6dp2;
  d2svdtdp = -d2term1dtdp - d2term2dtdp - d2term3dtdp + d2term4dtdp + d2term5dtdp - d2term6dtdp;
}
//
//-----------------------------------------------------------------------------
// derived properties
//-----------------------------------------------------------------------------
//
void
DIFF_cp_tp_L_K(double t, double p, double & cp, double & dcpdt, double & dcpdp)
{
  double h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;

  DIFF_h_tp_L_K(t, p, h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);

  cp = dhdt;
  dcpdt = d2hdt2;
  dcpdp = d2hdtdp;
}
//
void
DIFF_cp_tp_G_K(double t, double p, double & cp, double & dcpdt, double & dcpdp)
{
  double h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;

  DIFF_h_tp_G_K(t, p, h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);

  cp = dhdt;
  dcpdt = d2hdt2;
  dcpdp = d2hdtdp;
}
//
void
DIFF_cv_tp_L_K(double t, double p, double & cv, double & dcvdt, double & dcvdp)
{
  double v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  double dudt_p, d2udt2_p, dudp_t, d2udtdp, d2udp2;

  static const double atmft3_toBtu = (101325. * 0.3048 * 0.3048 * 0.3048) / 1055.05585262;

  DIFF_v_tp_L_K(t, p, v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
  DIFF_h_tp_L_K(t, p, h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);

  dudt_p = dhdt - p * dvdt * atmft3_toBtu;
  d2udt2_p = d2hdt2 - p * d2vdt2 * atmft3_toBtu;
  dudp_t = dhdp - (v + p * dvdp) * atmft3_toBtu;
  d2udtdp = d2hdtdp - (dvdt + p * d2vdtdp) * atmft3_toBtu;
  d2udp2 = d2hdp2 - (dvdp + dvdp + p * d2vdp2) * atmft3_toBtu;

  cv = dudt_p - dudp_t * dvdt / dvdp;
  dcvdt = d2udt2_p -
          (d2udtdp * dvdt / dvdp + dudp_t * (d2vdt2 * dvdp - dvdt * d2vdtdp) / (dvdp * dvdp));
  dcvdp =
      d2udtdp - (d2udp2 * dvdt / dvdp + dudp_t * (d2vdtdp * dvdp - dvdt * d2vdp2) / (dvdp * dvdp));
}
//
void
DIFF_cv_tp_G_K(double t, double p, double & cv, double & dcvdt, double & dcvdp)
{
  double v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  double dudt_p, d2udt2_p, dudp_t, d2udtdp, d2udp2;

  static const double atmft3_toBtu = (101325. * 0.3048 * 0.3048 * 0.3048) / 1055.05585262;

  DIFF_v_tp_G_K(t, p, v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
  DIFF_h_tp_G_K(t, p, h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);

  dudt_p = dhdt - p * dvdt * atmft3_toBtu;
  d2udt2_p = d2hdt2 - p * d2vdt2 * atmft3_toBtu;
  dudp_t = dhdp - (v + p * dvdp) * atmft3_toBtu;
  d2udtdp = d2hdtdp - (dvdt + p * d2vdtdp) * atmft3_toBtu;
  d2udp2 = d2hdp2 - (dvdp + dvdp + p * d2vdp2) * atmft3_toBtu;

  cv = dudt_p - dudp_t * dvdt / dvdp;
  dcvdt = d2udt2_p -
          (d2udtdp * dvdt / dvdp + dudp_t * (d2vdt2 * dvdp - dvdt * d2vdtdp) / (dvdp * dvdp));
  dcvdp =
      d2udtdp - (d2udp2 * dvdt / dvdp + dudp_t * (d2vdtdp * dvdp - dvdt * d2vdp2) / (dvdp * dvdp));
}
//
void
DIFF_w_tp_L_K(double t, double p, double & w, double & dwdt, double & dwdp)
{
  double v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double s, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp;

  static const double atmft_lb_to_ft2_s_2 = 101325. * 0.3048 / 0.45359237;

  DIFF_v_tp_L_K(t, p, v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
  DIFF_s_tp_L_K(t, p, s, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp);

  double dvdp_s = dvdp - dvdt * dsdp / dsdt;
  double dvpsdt =
      d2vdtdp - d2vdt2 * dsdp / dsdt - dvdt * (d2sdtdp * dsdt - dsdp * d2sdt2) / (dsdt * dsdt);
  double dvpsdp =
      d2vdp2 - d2vdtdp * dsdp / dsdt - dvdt * (d2sdp2 * dsdt - dsdp * d2sdtdp) / (dsdt * dsdt);
  double dpdv_s = 1. / dvdp_s;

  double sqrt_mdpdv_s = sqrt(-dpdv_s * atmft_lb_to_ft2_s_2);
  w = v * sqrt_mdpdv_s;
  dwdt = dvdt * sqrt_mdpdv_s +
         v / (2. * sqrt_mdpdv_s) * dvpsdt / (dvdp_s * dvdp_s) * atmft_lb_to_ft2_s_2;
  dwdp = dvdp * sqrt_mdpdv_s +
         v / (2. * sqrt_mdpdv_s) * dvpsdp / (dvdp_s * dvdp_s) * atmft_lb_to_ft2_s_2;
}
//
void
DIFF_w_tp_G_K(double t, double p, double & w, double & dwdt, double & dwdp)
{
  double v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double s, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp;

  static const double atmft_lb_to_ft2_s_2 = 101325. * 0.3048 / 0.45359237;

  DIFF_v_tp_G_K(t, p, v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
  DIFF_s_tp_G_K(t, p, s, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp);

  double dvdp_s = dvdp - dvdt * dsdp / dsdt;
  double dvpsdt =
      d2vdtdp - d2vdt2 * dsdp / dsdt - dvdt * (d2sdtdp * dsdt - dsdp * d2sdt2) / (dsdt * dsdt);
  double dvpsdp =
      d2vdp2 - d2vdtdp * dsdp / dsdt - dvdt * (d2sdp2 * dsdt - dsdp * d2sdtdp) / (dsdt * dsdt);
  double dpdv_s = 1. / dvdp_s;

  double sqrt_mdpdv_s = sqrt(-dpdv_s * atmft_lb_to_ft2_s_2);
  w = v * sqrt_mdpdv_s;
  dwdt = dvdt * sqrt_mdpdv_s +
         v / (2. * sqrt_mdpdv_s) * dvpsdt / (dvdp_s * dvdp_s) * atmft_lb_to_ft2_s_2;
  dwdp = dvdp * sqrt_mdpdv_s +
         v / (2. * sqrt_mdpdv_s) * dvpsdp / (dvdp_s * dvdp_s) * atmft_lb_to_ft2_s_2;
}
//
//-----------------------------------------------------------------------------
// flash routines
//-----------------------------------------------------------------------------
//
int
FLASH_prho_L_K(double p, double rho, double & t)
{
  static const double tol_v = 1.e-6;
  static const double tmin = 250. * 9. / 5.;

  double ts, dtsdp, d2tsdp2;
  double vl, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;

  DIFF_ts_p_K(p, ts, dtsdp, d2tsdp2);
  t = ts;

  double v = 1. / rho;

  int it = 0;
  double dv = 1.e12;
  while (fabs(dv) / v > tol_v)
  {
    DIFF_v_tp_L_K(t, p, vl, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
    dv = vl - v;
    t -= dv / dvdt;
    if(t < tmin) t = tmin;
    if (++it > 20)
    {
      return -1;
    }
  }
  return 0;
}
//
void
DERIV_prho_L_K(double t, double p, double & du_dp_rho, double & du_drho_p)
{
  static const double atmft3_toBtu = (101325. * 0.3048 * 0.3048 * 0.3048) / 1055.05585262;

  double v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  double dudt, dudp;
  double dp_dv_u, dp_du_v;

  DIFF_v_tp_L_K(t, p, v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
  DIFF_h_tp_L_K(t, p, h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);
  // u = h - p * v * atmft3_toBtu;
  dudt = dhdt - p * dvdt * atmft3_toBtu;
  dudp = dhdp - (v + p * dvdp) * atmft3_toBtu;

  dp_dv_u = 1. / (dvdp - dvdt * dudp / dudt);
  dp_du_v = 1. / (dudp - dudt * dvdp / dvdt);

  du_dp_rho = 1. / dp_du_v;
  du_drho_p = v * v * dp_dv_u / dp_du_v;
}
//
int
FLASH_prho_G_K(double p, double rho, double & t)
{
  static const double tol_v = 1.e-6;
  static const double tmin = 250. * 9. / 5.;

  double ts, dtsdp, d2tsdp2;
  double vl, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;

  DIFF_ts_p_K(p, ts, dtsdp, d2tsdp2);
  t = ts;

  double v = 1. / rho;

  int it = 0;
  double dv = 1.e12;
  while (fabs(dv) / v > tol_v)
  {
    DIFF_v_tp_G_K(t, p, vl, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
    dv = vl - v;
    t -= dv / dvdt;
    if(t < tmin) t = tmin;
    if (++it > 20)
    {
      return -1;
    }
  }
  return 0;
}
//
void
DERIV_prho_G_K(double t, double p, double & du_dp_rho, double & du_drho_p)
{
  static const double atmft3_toBtu = (101325. * 0.3048 * 0.3048 * 0.3048) / 1055.05585262;

  double v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  double dudt, dudp;
  double dp_dv_u, dp_du_v;

  DIFF_v_tp_G_K(t, p, v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
  DIFF_h_tp_G_K(t, p, h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);
  // u = h - p * v * atmft3_toBtu;
  dudt = dhdt - p * dvdt * atmft3_toBtu;
  dudp = dhdp - (v + p * dvdp) * atmft3_toBtu;

  dp_dv_u = 1. / (dvdp - dvdt * dudp / dudt);
  dp_du_v = 1. / (dudp - dudt * dvdp / dvdt);

  du_dp_rho = 1. / dp_du_v;
  du_drho_p = v * v * dp_dv_u / dp_du_v;
}
//
int
FLASH_ph_L_K(double p, double h, double & t)
{
  static const double tol_h = 1.e-6;
  static const double tmin = 250. * 9. / 5.;

  double ts, dtsdp, d2tsdp2;
  double hl, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;

  DIFF_ts_p_K(p, ts, dtsdp, d2tsdp2);
  t = ts;

  int it = 0;
  double dh = 1.e12;
  while (fabs(dh) > tol_h)
  {
    DIFF_h_tp_L_K(t, p, hl, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);
    dh = hl - h;
    t -= dh / dhdt;
    if(t < tmin) t = tmin;
    if (++it > 20)
    {
      return -1;
    }
  }
  return 0;
}
//
int
FLASH_ph_G_K(double p, double h, double & t)
{
  static const double tol_h = 1.e-6;
  static const double tmin = 250. * 9. / 5.;

  double ts, dtsdp, d2tsdp2;
  double hv, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;

  DIFF_ts_p_K(p, ts, dtsdp, d2tsdp2);
  t = ts;

  int it = 0;
  double dh = 1.e12;
  while (fabs(dh) > tol_h)
  {
    DIFF_h_tp_G_K(t, p, hv, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);
    dh = hv - h;
    t -= dh / dhdt;
    if(t < tmin) t = tmin;
    if (++it > 20)
    {
      return -1;
    }
  }
  return 0;
}
//
int
FLASH_ps_L_K(double p, double s, double & t)
{
  static const double tol_s = 1.e-8;
  static const double tmin = 250. * 9. / 5.;

  double ts, dtsdp, d2tsdp2;
  double sl, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp;

  DIFF_ts_p_K(p, ts, dtsdp, d2tsdp2);
  t = ts;

  int it = 0;
  double ds = 1.e12;
  while (fabs(ds) > tol_s)
  {
    DIFF_s_tp_L_K(t, p, sl, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp);
    ds = sl - s;
    t -= ds / dsdt;
    if(t < tmin) t = tmin;
    if (++it > 20)
    {
      return -1;
    }
  }
  return 0;
}
//
int
FLASH_ps_G_K(double p, double s, double & t)
{
  static const double tol_s = 1.e-8;
  static const double tmin = 250. * 9. / 5.;

  double ts, dtsdp, d2tsdp2;
  double sv, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp;

  DIFF_ts_p_K(p, ts, dtsdp, d2tsdp2);
  t = ts;

  int it = 0;
  double ds = 1.e12;
  while (fabs(ds) > tol_s)
  {
    DIFF_s_tp_G_K(t, p, sv, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp);
    ds = sv - s;
    t -= ds / dsdt;
    if(t < tmin) t = tmin;
    if (++it > 20)
    {
      return -1;
    }
  }
  return 0;
}
//
int
FLASH_vu_L_K(double v, double u, double & t, double & p)
{
  static const double tol_v = 1.e-6;
  static const double tol_u = 1.e-6;
  static const double tmin = 250. * 9. / 5.;
  static const double pmin = 1.e-5;
  static const double atmft3_toBtu = (101325. * 0.3048 * 0.3048 * 0.3048) / 1055.05585262;

  double vl, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double hl, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  double ul;

  // initial guess from normal boiling point in R
  static const double tnb = 1029.73 * 9. / 5.;
  t = tnb;
  p = 1.;

  int it = 0;
  double dv = 1.e12, du = 1.e12;
  double ddudt, ddudp, ddvdt, ddvdp, den;
  while (fabs(dv) / v > tol_v || fabs(du) > tol_u)
  {
    DIFF_v_tp_L_K(t, p, vl, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
    DIFF_h_tp_L_K(t, p, hl, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);
    ul = hl - p * vl * atmft3_toBtu;
    dv = vl - v;
    du = ul - u;
    ddvdt = dvdt;
    ddvdp = dvdp;
    ddudt = dhdt - p * dvdt * atmft3_toBtu;
    ddudp = dhdp - (vl + p * dvdp) * atmft3_toBtu;
    den = ddvdp * ddudt - ddvdt * ddudp;
    t -= (du * ddvdp - dv * ddudp) / den;
    p -= (dv * ddudt - du * ddvdt) / den;
    if(t < tmin) t = tmin;
    if (p < pmin)
      p = pmin;
    if (++it > 50)
    {
      return -1;
    }
  }
  return 0;
}
//
void
DERIV_vu_L_K(
    double t, double p, double & dp_dv_u, double & dp_du_v, double & dt_dv_u, double & dt_du_v)
{
  static const double atmft3_toBtu = (101325. * 0.3048 * 0.3048 * 0.3048) / 1055.05585262;

  double vl, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double hl, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  double dudt, dudp;

  DIFF_v_tp_L_K(t, p, vl, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
  DIFF_h_tp_L_K(t, p, hl, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);
  // ul = hl - p * vl * atmft3_toBtu;
  dudt = dhdt - p * dvdt * atmft3_toBtu;
  dudp = dhdp - (vl + p * dvdp) * atmft3_toBtu;

  dp_dv_u = 1. / (dvdp - dvdt * dudp / dudt);
  dp_du_v = 1. / (dudp - dudt * dvdp / dvdt);
  dt_dv_u = 1. / (dvdt - dvdp * dudt / dudp);
  dt_du_v = 1. / (dudt - dudp * dvdt / dvdp);
}
//
int
FLASH_vu_G_K(double v, double u, double & t, double & p)
{
  static const double tol_v = 1.e-6;
  static const double tol_u = 1.e-6;
  static const double tmin = 250. * 9. / 5.;
  static const double pmin = 1.e-5;
  static const double atmft3_toBtu = (101325. * 0.3048 * 0.3048 * 0.3048) / 1055.05585262;

  double vv, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double hv, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  double uv;

  // initial guess from normal boiling point in R
  static const double tnb = 1029.73 * 9. / 5.;
  t = tnb;
  p = 1.;

  int it = 0;
  double dv = 1.e12, du = 1.e12;
  double ddudt, ddudp, ddvdt, ddvdp, den;
  while (fabs(dv) / v > tol_v || fabs(du) > tol_u)
  {
    DIFF_v_tp_G_K(t, p, vv, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
    DIFF_h_tp_G_K(t, p, hv, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);
    uv = hv - p * vv * atmft3_toBtu;
    dv = vv - v;
    du = uv - u;
    ddvdt = dvdt;
    ddvdp = dvdp;
    ddudt = dhdt - p * dvdt * atmft3_toBtu;
    ddudp = dhdp - (vv + p * dvdp) * atmft3_toBtu;
    den = ddvdp * ddudt - ddvdt * ddudp;
    t -= (du * ddvdp - dv * ddudp) / den;
    p -= (dv * ddudt - du * ddvdt) / den;
    if(t < tmin) t = tmin;
    if (p < pmin)
      p = pmin;
    if (++it > 50)
    {
      return -1;
    }
  }
  return 0;
}
//
void
DERIV_vu_G_K(
    double t, double p, double & dp_dv_u, double & dp_du_v, double & dt_dv_u, double & dt_du_v)
{
  static const double atmft3_toBtu = (101325. * 0.3048 * 0.3048 * 0.3048) / 1055.05585262;

  double vv, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double hv, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  double dudt, dudp;

  DIFF_v_tp_G_K(t, p, vv, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
  DIFF_h_tp_G_K(t, p, hv, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);
  // uv = hv - p * vv * atmft3_toBtu;
  dudt = dhdt - p * dvdt * atmft3_toBtu;
  dudp = dhdp - (vv + p * dvdp) * atmft3_toBtu;

  dp_dv_u = 1. / (dvdp - dvdt * dudp / dudt);
  dp_du_v = 1. / (dudp - dudt * dvdp / dvdt);
  dt_dv_u = 1. / (dvdt - dvdp * dudt / dudp);
  dt_du_v = 1. / (dudt - dudp * dvdt / dvdp);
}
//
int
FLASH_vh_L_K(double v, double h, double & t, double & p)
{
  static const double tol_v = 1.e-6;
  static const double tol_h = 1.e-6;
  static const double tmin = 250. * 9. / 5.;
  static const double pmin = 1.e-5;

  double vl, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double hl, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;

  // initial guess from normal boiling point in R
  static const double tnb = 1029.73 * 9. / 5.;
  t = tnb;
  p = 1.;

  int it = 0;
  double dv = 1.e12, dh = 1.e12;
  double ddhdt, ddhdp, ddvdt, ddvdp, den;
  while (fabs(dv) / v > tol_v || fabs(dh) > tol_h)
  {
    DIFF_v_tp_L_K(t, p, vl, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
    DIFF_h_tp_L_K(t, p, hl, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);
    dv = vl - v;
    dh = hl - h;
    ddvdt = dvdt;
    ddvdp = dvdp;
    ddhdt = dhdt;
    ddhdp = dhdp;
    den = ddvdp * ddhdt - ddvdt * ddhdp;
    t -= (dh * ddvdp - dv * ddhdp) / den;
    p -= (dv * ddhdt - dh * ddvdt) / den;
    if(t < tmin) t = tmin;
    if (p < pmin)
      p = pmin;
    if (++it > 50)
    {
      return -1;
    }
  }
  return 0;
}
//
void
DERIV_vh_L_K(
    double t, double p, double & dp_dv_h, double & dp_dh_v, double & dt_dv_h, double & dt_dh_v)
{
  double vl, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double hl, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;

  DIFF_v_tp_L_K(t, p, vl, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
  DIFF_h_tp_L_K(t, p, hl, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);

  dp_dv_h = 1. / (dvdp - dvdt * dhdp / dhdt);
  dp_dh_v = 1. / (dhdp - dhdt * dvdp / dvdt);
  dt_dv_h = 1. / (dvdt - dvdp * dhdt / dhdp);
  dt_dh_v = 1. / (dhdt - dhdp * dvdt / dvdp);
}
//
int
FLASH_vh_G_K(double v, double h, double & t, double & p)
{
  static const double tol_v = 1.e-6;
  static const double tol_h = 1.e-6;
  static const double tmin = 250. * 9. / 5.;
  static const double pmin = 1.e-5;

  double vv, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double hv, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;

  // initial guess from normal boiling point in R
  static const double tnb = 1029.73 * 9. / 5.;
  t = tnb;
  p = 1.;

  int it = 0;
  double dv = 1.e12, dh = 1.e12;
  double ddhdt, ddhdp, ddvdt, ddvdp, den;
  while (fabs(dv) / v > tol_v || fabs(dh) > tol_h)
  {
    DIFF_v_tp_G_K(t, p, vv, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
    DIFF_h_tp_G_K(t, p, hv, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);
    dv = vv - v;
    dh = hv - h;
    ddvdt = dvdt;
    ddvdp = dvdp;
    ddhdt = dhdt;
    ddhdp = dhdp;
    den = ddvdp * ddhdt - ddvdt * ddhdp;
    t -= (dh * ddvdp - dv * ddhdp) / den;
    p -= (dv * ddhdt - dh * ddvdt) / den;
    if(t < tmin) t = tmin;
    if (p < pmin)
      p = pmin;
    if (++it > 50)
    {
      return -1;
    }
  }
  return 0;
}
//
void
DERIV_vh_G_K(
    double t, double p, double & dp_dv_h, double & dp_dh_v, double & dt_dv_h, double & dt_dh_v)
{
  double vv, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double hv, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;

  DIFF_v_tp_G_K(t, p, vv, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
  DIFF_h_tp_G_K(t, p, hv, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);

  dp_dv_h = 1. / (dvdp - dvdt * dhdp / dhdt);
  dp_dh_v = 1. / (dhdp - dhdt * dvdp / dvdt);
  dt_dv_h = 1. / (dvdt - dvdp * dhdt / dhdp);
  dt_dh_v = 1. / (dhdt - dhdp * dvdt / dvdp);
}
//
int
FLASH_hs_L_K(double h, double s, double & t, double & p)
{
  static const double tol_h = 1.e-6;
  static const double tol_s = 1.e-8;
  static const double tmin = 250. * 9. / 5.;
  static const double pmin = 1.e-5;

  double hl, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  double sl, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp;

  // initial guess from normal boiling point in R
  static const double tnb = 1029.73 * 9. / 5.;
  t = tnb;
  p = 1.;

  int it = 0;
  double dh = 1.e12, ds = 1.e12;
  double ddhdt, ddhdp, ddsdt, ddsdp, den;
  while (fabs(dh) > tol_h || fabs(ds) > tol_s)
  {
    DIFF_h_tp_L_K(t, p, hl, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);
    DIFF_s_tp_L_K(t, p, sl, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp);
    dh = hl - h;
    ds = sl - s;
    ddhdt = dhdt;
    ddhdp = dhdp;
    ddsdt = dsdt;
    ddsdp = dsdp;
    den = ddsdp * ddhdt - ddsdt * ddhdp;
    t -= (dh * ddsdp - ds * ddhdp) / den;
    p -= (ds * ddhdt - dh * ddsdt) / den;
    if(t < tmin) t = tmin;
    if (p < pmin)
      p = pmin;
    if (++it > 50)
    {
      return -1;
    }
  }
  return 0;
}
//
int
FLASH_hs_G_K(double h, double s, double & t, double & p)
{
  static const double tol_h = 1.e-6;
  static const double tol_s = 1.e-8;
  static const double tmin = 250. * 9. / 5.;
  static const double pmin = 1.e-5;

  double hv, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  double sv, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp;

  // initial guess from normal boiling point in R
  static const double tnb = 1029.73 * 9. / 5.;
  t = tnb;
  p = 1.;

  int it = 0;
  double dh = 1.e12, ds = 1.e12;
  double ddhdt, ddhdp, ddsdt, ddsdp, den;
  while (fabs(dh) > tol_h || fabs(ds) > tol_s)
  {
    DIFF_h_tp_G_K(t, p, hv, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);
    DIFF_s_tp_G_K(t, p, sv, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp);
    dh = hv - h;
    ds = sv - s;
    ddhdt = dhdt;
    ddhdp = dhdp;
    ddsdt = dsdt;
    ddsdp = dsdp;
    den = ddsdp * ddhdt - ddsdt * ddhdp;
    t -= (dh * ddsdp - ds * ddhdp) / den;
    p -= (ds * ddhdt - dh * ddsdt) / den;
    if(t < tmin) t = tmin;
    if (p < pmin)
      p = pmin;
    if (++it > 50)
    {
      return -1;
    }
  }
  return 0;
}
//
//-----------------------------------------------------------------------------
// transport properties
//-----------------------------------------------------------------------------
//
void
sigma_t_K(double t, double & sigma, double & dsigmadt)
{ // K, mN/m
  // dyne = 1.e-5 N; dyne/cm = mN/m
  static const double a = 115.7;
  static const double b =-0.064;
  //
  // according to NaK Handbook
  sigma = a+(t-273.15)*b;
  dsigmadt = b;
}
//
double
etal_tv_K(double t, double v)
{ // R, lb/ft3, lbm/ft-hr
    static const double tKg = 380.+273.15;
    static const double third = 1./3.;
    static const double lb_ft3_to_m3_kg=(0.3048*0.3048*0.3048) / 0.45359237;
    static const double cP_to_lbm_ft_hr=1.e-3* (3600. * 0.3048) /0.45359237;

    double eta;
    double tK = 5./9.*t;
    double rho = 1.e-3/(v*lb_ft3_to_m3_kg); //g/cm3

    if(tK<tKg) {
        eta=0.1131*pow(rho,third)*exp(680.*rho/tK); //cP=mPa*s
        eta=eta*cP_to_lbm_ft_hr;
    } else {
        eta=0.0799*pow(rho,third)*exp(978.*rho/tK);
        eta=eta*cP_to_lbm_ft_hr;
    }
    return eta;
}
//
void
etav_t_K(double t, double & etav, double & detavdt, double & d2etavdt2)
{ // R, lbm/ft-hr
  static const double a = 7.65637393e-03;
  static const double b = 1.81419228e-05;
  static const double c =-4.97899269e-10;
  //
  etav = a + t * (b + t * c);
  detavdt = b + 2. * t * c;
  d2etavdt2 = 2. * c;
}
//
void
lambdal_t_K(double t, double & lambdal, double & dlambdaldt, double & d2lambdaldt2)
{ // R, Btu/hr-ft-F
  static const double a = 0.438;
  static const double b =-2.22e-4;
  static const double c = 39.5;
  static const double d = 273.2;
  static const double W_cm_K_to_Btu_hr_ft_F = (1.e2*3600. * 0.3048 * 5. / 9.)/1055.05585262;
  static const double dtCdt=5./9.;
  //
  double tC = dtCdt*t - 273.15; //Rankine to Celsius
  //
  lambdal      = W_cm_K_to_Btu_hr_ft_F*(a + b*tC + c/(tC+d));
  dlambdaldt   = W_cm_K_to_Btu_hr_ft_F*(b - c/((tC+d)*(tC+d)))*dtCdt;
  d2lambdaldt2 = W_cm_K_to_Btu_hr_ft_F*(2.*c/((tC+d)*(tC+d)*(tC+d)))*dtCdt*dtCdt;
}
//
void
lambdav_t_K(double t, double & lambdav, double & dlambdavdt, double & d2lambdavdt2)
{ // R, Btu/hr-ft-F
  static const double a = 1.96650412e-02;
  static const double b =-5.61168099e-05;
  static const double c = 7.08532889e-08;
  static const double d =-3.83063201e-11;
  static const double e = 1.06962032e-14;
  static const double f =-1.51741453e-18;
  static const double g = 8.69047448e-23;
  //
  lambdav = a + t * (b + t * (c + t * (d + t * (e + t * (f + t * g)))));
  dlambdavdt = b + t * (2. * c + t * (3. * d + t * (4. * e + t * (5. * f + t * 6. * g))));
  d2lambdavdt2 = 2. * c + t * (6. * d + t * (12. * e + t * (20. * f + t * 30. * g)));
}
