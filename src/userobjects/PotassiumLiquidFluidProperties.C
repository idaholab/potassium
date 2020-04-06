#include "PotassiumLiquidFluidProperties.h"
#include "contrib/libPotassium/K_Golden.h"

registerMooseObject("PotassiumApp", PotassiumLiquidFluidProperties);

InputParameters
PotassiumLiquidFluidProperties::validParams()
{
  InputParameters params = SinglePhaseFluidProperties::validParams();
  params += NaNInterface::validParams();
  params.addClassDescription("Fluid properties of potassium vapor.");
  return params;
}

PotassiumLiquidFluidProperties::PotassiumLiquidFluidProperties(const InputParameters & parameters)
  : SinglePhaseFluidProperties(parameters),
    NaNInterface(this),
    _to_ft(1. / 0.3048),
    _to_m(1. / _to_ft),
    _to_ft2(_to_ft * _to_ft),
    _to_m2(1. / _to_ft2),
    _to_ft3(_to_ft * _to_ft * _to_ft),
    _to_m3(1. / _to_ft3),
    _to_lb(1. / 0.45359237),
    _to_kg(1. / _to_lb),
    _to_atm(1. / 101325.),
    _to_Pa(101325.),
    _to_R(9. / 5.),
    _to_K(5. / 9.),
    _to_Btu(1. / 1055.05585262),
    _to_J(1. / _to_Btu),
    _to_s(3600.),
    _to_ft3_lb(_to_ft3 / _to_lb),
    _to_m3_kg(1. / _to_ft3_lb),
    _to_Btu_lb(_to_Btu / _to_lb),
    _to_J_kg(1. / _to_Btu_lb),
    _to_Btu_lbR(_to_Btu / (_to_lb * _to_R)),
    _to_J_kgK(1. / _to_Btu_lbR)
{
}

Real
PotassiumLiquidFluidProperties::p_from_v_e(Real v, Real e) const
{
  v *= _to_ft3_lb;
  e *= _to_Btu_lb;

  double p, T;
  FLASH_vu_L_K(v, e, T, p);

  return p * _to_Pa;
}

void
PotassiumLiquidFluidProperties::p_from_v_e(
    Real v, Real e, Real & p, Real & dp_dv, Real & dp_de) const
{
  v *= _to_ft3_lb;
  e *= _to_Btu_lb;

  double T;
  double dp_dv_u, dp_du_v, dt_dv_u, dt_du_v;
  FLASH_vu_L_K(v, e, T, p);
  DERIV_vu_L_K(T, p, dp_dv_u, dp_du_v, dt_dv_u, dt_du_v);

  p *= _to_Pa;
  dp_dv = dp_dv_u * _to_Pa / _to_m3_kg;
  dp_de = dp_du_v * _to_Pa / _to_J_kg;
}

Real
PotassiumLiquidFluidProperties::T_from_v_e(Real v, Real e) const
{
  v *= _to_ft3_lb;
  e *= _to_Btu_lb;

  double p, T;
  FLASH_vu_L_K(v, e, T, p);

  return T * _to_K;
}

void
PotassiumLiquidFluidProperties::T_from_v_e(
    Real v, Real e, Real & T, Real & dT_dv, Real & dT_de) const
{
  v *= _to_ft3_lb;
  e *= _to_Btu_lb;

  double p;
  double dp_dv_u, dp_du_v, dt_dv_u, dt_du_v;
  FLASH_vu_L_K(v, e, T, p);
  DERIV_vu_L_K(T, p, dp_dv_u, dp_du_v, dt_dv_u, dt_du_v);

  T *= _to_K;
  dT_dv = dt_dv_u * _to_K / _to_m3_kg;
  dT_de = dt_du_v * _to_K / _to_J_kg;
}

Real
PotassiumLiquidFluidProperties::c_from_v_e(Real v, Real e) const
{
  v *= _to_ft3_lb;
  e *= _to_Btu_lb;

  double p, T;
  double w, dwdt, dwdp;
  FLASH_vu_L_K(v, e, T, p);
  DIFF_w_tp_L_K(T, p, w, dwdt, dwdp);

  return w * _to_m;
}

void
PotassiumLiquidFluidProperties::c_from_v_e(
    Real v, Real e, Real & c, Real & dc_dv, Real & dc_de) const
{
  v *= _to_ft3_lb;
  e *= _to_Btu_lb;

  double p, T;
  double w, dwdt, dwdp;
  double v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double u_, dudt, d2udt2, dudp, d2udp2, d2udtdp;
  FLASH_vu_L_K(v, e, T, p);
  DIFF_vu_tp_L_K(
      T, p, v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp, u_, dudt, d2udt2, dudp, d2udp2, d2udtdp);
  DIFF_w_tp_L_K(T, p, w, dwdt, dwdp);

  c = w * _to_m;
  dc_dv = (dwdt * dudp - dwdp * dudt) / (dvdt * dudp - dvdp * dudt) * _to_m / _to_m3_kg;
  dc_de = (dwdt * dvdp - dwdp * dvdt) / (dudt * dvdp - dudp * dvdt) * _to_m / _to_J_kg;
}

Real
PotassiumLiquidFluidProperties::cp_from_v_e(Real v, Real e) const
{
  v *= _to_ft3_lb;
  e *= _to_Btu_lb;

  double p, T;
  double cp, dcpdt, dcpdp;
  FLASH_vu_L_K(v, e, T, p);

  DIFF_cp_tp_L_K(T, p, cp, dcpdt, dcpdp);

  return cp * _to_J_kgK;
}

void
PotassiumLiquidFluidProperties::cp_from_v_e(
    Real v, Real e, Real & cp, Real & dcp_dv, Real & dcp_de) const
{
  v *= _to_ft3_lb;
  e *= _to_Btu_lb;

  double p, T;
  double dcpdt, dcpdp;
  double v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double u_, dudt, d2udt2, dudp, d2udp2, d2udtdp;
  FLASH_vu_L_K(v, e, T, p);

  DIFF_cp_tp_L_K(T, p, cp, dcpdt, dcpdp);
  DIFF_vu_tp_L_K(
      T, p, v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp, u_, dudt, d2udt2, dudp, d2udp2, d2udtdp);

  dcp_dv = (dcpdt * dudp - dcpdp * dudt) / (dvdt * dudp - dvdp * dudt);
  dcp_de = (dcpdt * dvdp - dcpdp * dvdt) / (dudt * dvdp - dudp * dvdt);

  cp *= _to_J_kgK;
  dcp_dv *= _to_J_kgK / _to_m3_kg;
  dcp_de *= _to_J_kgK / _to_J_kg;
}

Real
PotassiumLiquidFluidProperties::cv_from_v_e(Real v, Real e) const
{
  v *= _to_ft3_lb;
  e *= _to_Btu_lb;

  double p, T;
  double cv, dcvdt, dcvdp;
  FLASH_vu_L_K(v, e, T, p);

  DIFF_cv_tp_L_K(T, p, cv, dcvdt, dcvdp);

  return cv * _to_J_kgK;
}

Real
PotassiumLiquidFluidProperties::mu_from_v_e(Real v, Real e) const
{
  v *= _to_ft3_lb;
  e *= _to_Btu_lb;

  double p, T;
  double eta;
  FLASH_vu_L_K(v, e, T, p);

  eta = etal_tv_K(T, v); // eta in lbm/ft-hr

  return eta * _to_kg / (_to_s * _to_m);
}

void
PotassiumLiquidFluidProperties::mu_from_v_e(
    Real v, Real e, Real & mu, Real & dmu_dv, Real & dmu_de) const
{
  v *= _to_ft3_lb;
  e *= _to_Btu_lb;

  double p, T;
  int ierr = FLASH_vu_L_K(v, e, T, p);
  if (ierr != 0)
  {
    mu = getNaN();
    dmu_dv = getNaN();
    dmu_de = getNaN();
  }
  else
  {
    double dmudt, dmudv;
    etal_tv_K(T, v, mu, dmudt, dmudv); // eta in lbm/ft-hr

    double v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
    double u_, dudt, d2udt2, dudp, d2udp2, d2udtdp;
    DIFF_vu_tp_L_K(
        T, p, v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp, u_, dudt, d2udt2, dudp, d2udp2, d2udtdp);

    dmu_dv = dmudt * dudp / (dvdt * dudp - dvdp * dudt) + dmudv;
    dmu_de = dmudt * dvdp / (dudt * dvdp - dudp * dvdt);

    const double mu_conv = _to_kg / (_to_s * _to_m);

    mu *= mu_conv;
    dmu_dv *= mu_conv / _to_m3_kg;
    dmu_de *= mu_conv / _to_J_kg;
  }
}

Real
PotassiumLiquidFluidProperties::mu_from_p_T(Real p, Real T) const
{
  p *= _to_atm;
  T *= _to_R;

  double v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double eta;

  DIFF_v_tp_L_K(T, p, v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);

  eta = etal_tv_K(T, v); // eta in lbm/ft-hr

  return eta * _to_kg / (_to_s * _to_m);
}

Real
PotassiumLiquidFluidProperties::k_from_v_e(Real v, Real e) const
{
  v *= _to_ft3_lb;
  e *= _to_Btu_lb;

  double p, T;
  double lambda, dlambdadt, d2lambdadt2;
  FLASH_vu_L_K(v, e, T, p);

  lambdal_t_K(T, lambda, dlambdadt, d2lambdadt2); // lambda in Btu/hr-ft-F

  return lambda * _to_J / (_to_s * _to_m * _to_K);
}

void
PotassiumLiquidFluidProperties::k_from_v_e(
    Real v, Real e, Real & k, Real & dk_dv, Real & dk_de) const
{
  v *= _to_ft3_lb;
  e *= _to_Btu_lb;

  double p, T;
  FLASH_vu_L_K(v, e, T, p);

  double dkdt, d2kdt2;
  lambdal_t_K(T, k, dkdt, d2kdt2); // k in Btu/hr-ft-F

  double v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double u_, dudt, d2udt2, dudp, d2udp2, d2udtdp;
  DIFF_vu_tp_L_K(
      T, p, v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp, u_, dudt, d2udt2, dudp, d2udp2, d2udtdp);

  dk_dv = dkdt * dudp / (dvdt * dudp - dvdp * dudt);
  dk_de = dkdt * dvdp / (dudt * dvdp - dudp * dvdt);

  const double k_conv = _to_J / (_to_s * _to_m * _to_K);

  k *= k_conv;
  dk_dv *= k_conv / _to_m3_kg;
  dk_de *= k_conv / _to_J_kg;
}

Real
PotassiumLiquidFluidProperties::s_from_v_e(Real v, Real e) const
{
  v *= _to_ft3_lb;
  e *= _to_Btu_lb;

  double p, T;
  double s, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp;
  FLASH_vu_L_K(v, e, T, p);

  DIFF_s_tp_L_K(T, p, s, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp);

  return s * _to_J_kgK;
}

void
PotassiumLiquidFluidProperties::s_from_v_e(
    Real v, Real e, Real & s, Real & ds_dv, Real & ds_de) const
{
  v *= _to_ft3_lb;
  e *= _to_Btu_lb;

  double p, T;
  double dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp;
  double v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double u_, dudt, d2udt2, dudp, d2udp2, d2udtdp;
  FLASH_vu_L_K(v, e, T, p);

  DIFF_s_tp_L_K(T, p, s, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp);
  DIFF_vu_tp_L_K(
      T, p, v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp, u_, dudt, d2udt2, dudp, d2udp2, d2udtdp);

  ds_dv = (dsdt * dudp - dsdp * dudt) / (dvdt * dudp - dvdp * dudt);
  ds_de = (dsdt * dvdp - dsdp * dvdt) / (dudt * dvdp - dudp * dvdt);

  s *= _to_J_kgK;
  ds_dv *= _to_J_kgK / _to_m3_kg;
  ds_de *= _to_J_kgK / _to_J_kg;
}

Real
PotassiumLiquidFluidProperties::s_from_h_p(Real h, Real p) const
{
  p *= _to_atm;
  h *= _to_Btu_lb;

  double T;
  double s, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp;

  FLASH_ph_L_K(p, h, T);
  DIFF_s_tp_L_K(T, p, s, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp);

  return s * _to_J_kgK;
}

void
PotassiumLiquidFluidProperties::s_from_h_p(
    Real h, Real p, Real & s, Real & ds_dh, Real & ds_dp) const
{
  p *= _to_atm;
  h *= _to_Btu_lb;

  double T;
  double dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp;
  double h_, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;

  FLASH_ph_L_K(p, h, T);
  DIFF_s_tp_L_K(T, p, s, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp);
  DIFF_h_tp_L_K(T, p, h_, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);

  s *= _to_J_kgK;
  ds_dp = (dsdp - dsdt * dhdp / dhdt) * _to_J_kgK / _to_Pa;
  ds_dh = dsdt / dhdt * _to_J_kgK / _to_J_kg;
}

Real
PotassiumLiquidFluidProperties::rho_from_p_s(Real p, Real s) const
{
  p *= _to_atm;
  s *= _to_Btu_lbR;

  double T;
  double v, dvdt, d2vdt2, dvdp, d2vdp2, d2vvdtdp;

  FLASH_ps_L_K(p, s, T);
  DIFF_v_tp_L_K(T, p, v, dvdt, d2vdt2, dvdp, d2vdp2, d2vvdtdp);

  return 1. / (v * _to_m3_kg);
}

void
PotassiumLiquidFluidProperties::rho_from_p_s(
    Real p, Real s, Real & rho, Real & drho_dp, Real & drho_ds) const
{
  p *= _to_atm;
  s *= _to_Btu_lbR;

  double T;
  double v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double s_, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp;
  double dv_dp_s, dv_ds_p;

  FLASH_ps_L_K(p, s, T);
  DIFF_v_tp_L_K(T, p, v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
  DIFF_s_tp_L_K(T, p, s_, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp);

  dv_dp_s = dvdp - dvdt * dsdp / dsdt;
  dv_ds_p = dvdt / dsdt;

  v *= _to_m3_kg;
  dv_dp_s *= (_to_m3_kg / _to_Pa);
  dv_ds_p *= (_to_m3_kg / _to_J_kgK);

  rho = 1. / v;
  double drho_dv = -rho * rho;
  drho_dp = drho_dv * dv_dp_s;
  drho_ds = drho_dv * dv_ds_p;
}

Real
PotassiumLiquidFluidProperties::e_from_v_h(Real v, Real h) const
{
  v *= _to_ft3_lb;
  h *= _to_Btu_lb;

  static const double atmft3_toBtu = (101325. * 0.3048 * 0.3048 * 0.3048) / 1055.05585262;

  double p, T;
  double v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double h_, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  double u;
  FLASH_vh_L_K(v, h, T, p);

  DIFF_v_tp_L_K(T, p, v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
  DIFF_h_tp_L_K(T, p, h_, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);

  u = h - p * v * atmft3_toBtu;

  return u * _to_J_kg;
}

void
PotassiumLiquidFluidProperties::e_from_v_h(
    Real v, Real h, Real & e, Real & de_dv, Real & de_dh) const
{
  v *= _to_ft3_lb;
  h *= _to_Btu_lb;

  static const double atmft3_toBtu = (101325. * 0.3048 * 0.3048 * 0.3048) / 1055.05585262;

  double p, T;
  double v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double h_, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  double u, dudt, dudp;
  FLASH_vh_L_K(v, h, T, p);

  DIFF_v_tp_L_K(T, p, v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
  DIFF_h_tp_L_K(T, p, h_, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);

  u = h - p * v * atmft3_toBtu;
  dudt = dhdt - p * dvdt * atmft3_toBtu;
  dudp = dhdp - (v + p * dvdp) * atmft3_toBtu;

  e = u * _to_J_kg;
  de_dv = (dudt * dhdp - dudp * dhdt) / (dvdt * dhdp - dvdp * dhdt) * _to_J_kg / _to_m3_kg;
  de_dh = (dudt * dvdp - dudp * dvdt) / (dhdt * dvdp - dhdp * dvdt);
}

Real
PotassiumLiquidFluidProperties::rho_from_p_T(Real p, Real T) const
{
  p *= _to_atm;
  T *= _to_R;

  double v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;

  DIFF_v_tp_L_K(T, p, v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);

  return 1. / (v * _to_m3_kg);
}

void
PotassiumLiquidFluidProperties::rho_from_p_T(
    Real p, Real T, Real & rho, Real & drho_dp, Real & drho_dT) const
{
  p *= _to_atm;
  T *= _to_R;

  double v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;

  DIFF_v_tp_L_K(T, p, v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);

  v *= _to_m3_kg;
  dvdt *= _to_m3_kg / _to_K;
  dvdp *= _to_m3_kg / _to_Pa;

  rho = 1. / v;
  const double drho_dv = -rho * rho;
  drho_dp = drho_dv * dvdp;
  drho_dT = drho_dv * dvdt;
}

Real
PotassiumLiquidFluidProperties::e_from_p_rho(Real p, Real rho) const
{
  p *= _to_atm;
  rho /= _to_ft3_lb;
  double v = 1. / rho;

  static const double atmft3_toBtu = (101325. * 0.3048 * 0.3048 * 0.3048) / 1055.05585262;

  double T;
  double v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  double u;
  FLASH_prho_L_K(p, rho, T);
  DIFF_v_tp_L_K(T, p, v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
  DIFF_h_tp_L_K(T, p, h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);

  u = h - p * v * atmft3_toBtu;

  return u * _to_J_kg;
}

void
PotassiumLiquidFluidProperties::e_from_p_rho(
    Real p, Real rho, Real & e, Real & de_dp, Real & de_drho) const
{
  p *= _to_atm;
  rho /= _to_ft3_lb;
  double v = 1. / rho;

  static const double atmft3_toBtu = (101325. * 0.3048 * 0.3048 * 0.3048) / 1055.05585262;

  double T;
  double v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  double u, dudt, dudp;
  double de_dv;
  FLASH_prho_L_K(p, rho, T);
  DIFF_v_tp_L_K(T, p, v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
  DIFF_h_tp_L_K(T, p, h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);

  u = h - p * v * atmft3_toBtu;
  dudt = dhdt - p * dvdt * atmft3_toBtu;
  dudp = dhdp - (v + p * dvdp) * atmft3_toBtu;

  de_dp = dudp - dudt * dvdp / dvdt;
  de_dv = dudt / dvdt;
  de_drho = -de_dv * v * v;

  e = u * _to_J_kg;
  de_dp *= _to_J_kg / _to_Pa;
  de_drho *= _to_J_kg * _to_m3_kg;
}

Real
PotassiumLiquidFluidProperties::h_from_p_T(Real p, Real T) const
{
  p *= _to_atm;
  T *= _to_R;

  double h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  DIFF_h_tp_L_K(T, p, h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);

  return h * _to_J_kg;
}

void
PotassiumLiquidFluidProperties::h_from_p_T(
    Real p, Real T, Real & h, Real & dh_dp, Real & dh_dT) const
{
  p *= _to_atm;
  T *= _to_R;

  double dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  DIFF_h_tp_L_K(T, p, h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);

  h *= _to_J_kg;
  dh_dp = dhdp * _to_J_kg / _to_Pa;
  dh_dT = dhdt * _to_J_kg / _to_K;
}

Real
PotassiumLiquidFluidProperties::s_from_p_T(Real p, Real T) const
{
  p *= _to_atm;
  T *= _to_R;

  double s, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp;
  DIFF_s_tp_L_K(T, p, s, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp);

  return s * _to_J_kgK;
}

void
PotassiumLiquidFluidProperties::s_from_p_T(
    Real p, Real T, Real & s, Real & ds_dp, Real & ds_dT) const
{
  p *= _to_atm;
  T *= _to_R;

  double dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp;
  DIFF_s_tp_L_K(T, p, s, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp);

  s *= _to_J_kgK;
  ds_dp = dsdp * _to_J_kgK / _to_Pa;
  ds_dT = dsdt * _to_J_kgK / _to_K;
}

Real
PotassiumLiquidFluidProperties::p_from_h_s(Real h, Real s) const
{
  h *= _to_Btu_lb;
  s *= _to_Btu_lbR;

  double T, p;
  FLASH_hs_L_K(h, s, T, p);

  return p * _to_Pa;
}

void
PotassiumLiquidFluidProperties::p_from_h_s(
    Real h, Real s, Real & p, Real & dp_dh, Real & dp_ds) const
{
  h *= _to_Btu_lb;
  s *= _to_Btu_lbR;

  double T;
  double h_, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  double s_, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp;
  FLASH_hs_L_K(h, s, T, p);
  DIFF_h_tp_L_K(T, p, h_, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);
  DIFF_s_tp_L_K(T, p, s_, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp);

  dp_dh = 1. / (dhdp - dhdt * dsdp / dsdt);
  dp_ds = 1. / (dsdp - dsdt * dhdp / dhdt);

  p *= _to_Pa;
  dp_dh *= _to_Pa / _to_J_kg;
  dp_ds *= _to_Pa / _to_J_kgK;
}

Real
PotassiumLiquidFluidProperties::g_from_v_e(Real v, Real e) const
{
  v *= _to_ft3_lb;
  e *= _to_Btu_lb;

  static const double atmft3_toBtu = (101325. * 0.3048 * 0.3048 * 0.3048) / 1055.05585262;

  double T, p, g;
  double v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;
  double h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp;
  double s, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp;

  FLASH_vu_L_K(v, e, T, p);
  DIFF_v_tp_L_K(T, p, v_, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);
  DIFF_h_tp_L_K(T, p, h, dhdt, d2hdt2, dhdp, d2hdp2, d2hdtdp);
  DIFF_s_tp_L_K(T, p, s, dsdt, d2sdt2, dsdp, d2sdp2, d2sdtdp);

  h = e + p * v * atmft3_toBtu;
  g = h - T * s;

  return g * _to_J_kg;
}

Real
PotassiumLiquidFluidProperties::beta_from_p_T(Real p, Real T) const
{
  double rho, drho_dp, drho_dT;
  rho_from_p_T(p, T, rho, drho_dp, drho_dT);
  return -drho_dT / rho;
}

void
PotassiumLiquidFluidProperties::beta_from_p_T(
    Real p, Real T, Real & beta, Real & dbeta_dp, Real & dbeta_dT) const
{
  p *= _to_atm;
  T *= _to_R;

  double v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp;

  DIFF_v_tp_L_K(T, p, v, dvdt, d2vdt2, dvdp, d2vdp2, d2vdtdp);

  v *= _to_m3_kg;
  dvdt *= (_to_m3_kg / _to_K);
  d2vdt2 *= (_to_m3_kg / (_to_K * _to_K));
  dvdp *= (_to_m3_kg / _to_Pa);
  d2vdtdp *= (_to_m3_kg / (_to_K * _to_Pa));

  beta = dvdt / v;
  dbeta_dT = (d2vdt2 * v - dvdt * dvdt) / (v * v);
  dbeta_dp = (d2vdtdp * v - dvdt * dvdp) / (v * v);
}

// - molar mass depends on the presence of K, K2, and maybe even K4 (monomer, dimer, tetramer)
/*
Real
PotassiumLiquidFluidProperties::molarMass() const
{
  return xxx;
}
*/
Real
PotassiumLiquidFluidProperties::criticalTemperature() const
{
  return 2239.;
}

Real
PotassiumLiquidFluidProperties::criticalDensity() const
{
  return 192.;
}
