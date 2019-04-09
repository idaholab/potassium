#ifndef POTASSIUMLIQUIDFLUIDPROPERTIES_H
#define POTASSIUMLIQUIDFLUIDPROPERTIES_H

#include "SinglePhaseFluidProperties.h"
#include "NaNInterface.h"

class PotassiumLiquidFluidProperties;

template <>
InputParameters validParams<PotassiumLiquidFluidProperties>();

/**
 * Fluid properties of potassium according to Golden et al.,
 * NaK Handbook, Pasternak, Lee & Bonilla with some
 * modifications regarding the liquid phase (compressibility)
 *
 * The following are the units used in Na_Golden:
 * - pressure:           atm
 * - temperature:        R
 * - specific energy:    Btu/lb
 * - density:            lb/ft3
 * - speed of sound:     ft/s
 * - dynamic viscosity:  lb/(hr ft)
 * - conductivity:       Btu/(hr ft R)
 */
class PotassiumLiquidFluidProperties : public SinglePhaseFluidProperties,
                                       public NaNInterface
{
public:
  PotassiumLiquidFluidProperties(const InputParameters & parameters);

  virtual Real p_from_v_e(Real v, Real e) const override;
  virtual void p_from_v_e(Real v, Real e, Real & p, Real & dp_dv, Real & dp_de) const override;
  virtual Real T_from_v_e(Real v, Real e) const override;
  virtual void T_from_v_e(Real v, Real e, Real & T, Real & dT_dv, Real & dT_de) const override;
  virtual Real c_from_v_e(Real v, Real e) const override;
  virtual void c_from_v_e(Real v, Real e, Real & c, Real & dc_dv, Real & dc_de) const override;
  virtual Real cp_from_v_e(Real v, Real e) const override;
  virtual void cp_from_v_e(Real v, Real e, Real & cp, Real & dcp_dv, Real & dcp_de) const override;
  virtual Real cv_from_v_e(Real v, Real e) const override;
  virtual Real mu_from_v_e(Real v, Real e) const override;
  virtual Real mu_from_p_T(Real p, Real T) const override;
  virtual Real k_from_v_e(Real v, Real e) const override;
  virtual Real s_from_v_e(Real v, Real e) const override;
  virtual void s_from_v_e(Real v, Real e, Real & s, Real & ds_dv, Real & ds_de) const override;
  virtual Real s_from_p_T(Real p, Real T) const override;
  virtual void s_from_p_T(Real p, Real T, Real & s, Real & ds_dp, Real & ds_dT) const override;
  virtual Real s_from_h_p(Real h, Real p) const override;
  virtual void s_from_h_p(Real h, Real p, Real & s, Real & ds_dh, Real & ds_dp) const override;
  virtual Real rho_from_p_s(Real p, Real s) const override;
  virtual void
  rho_from_p_s(Real p, Real s, Real & rho, Real & drho_dp, Real & drho_ds) const override;
  virtual Real e_from_v_h(Real v, Real h) const override;
  virtual void e_from_v_h(Real v, Real h, Real & e, Real & de_dv, Real & de_dh) const override;
  virtual Real rho_from_p_T(Real p, Real T) const override;
  virtual void
  rho_from_p_T(Real p, Real T, Real & rho, Real & drho_dp, Real & drho_dT) const override;
  virtual Real e_from_p_rho(Real p, Real rho) const override;
  virtual void
  e_from_p_rho(Real p, Real rho, Real & e, Real & de_dp, Real & de_drho) const override;
  virtual Real h_from_p_T(Real p, Real T) const override;
  virtual void h_from_p_T(Real p, Real T, Real & h, Real & dh_dp, Real & dh_dT) const override;
  virtual Real p_from_h_s(Real h, Real s) const override;
  virtual void p_from_h_s(Real h, Real s, Real & p, Real & dp_dh, Real & dp_ds) const override;
  virtual Real g_from_v_e(Real v, Real e) const override;
  virtual Real beta_from_p_T(Real p, Real T) const override;
  virtual void
  beta_from_p_T(Real p, Real T, Real & beta, Real & dbeta_dp, Real & dbeta_dT) const override;
  // virtual Real molarMass() const override;
  virtual Real criticalTemperature() const override;
  virtual Real criticalDensity() const override;

protected:
  /// Conversion factor from m to ft
  const Real _to_ft;
  /// Conversion factor from ft to m
  const Real _to_m;
  /// Conversion factor from m2 to ft2
  const Real _to_ft2;
  /// Conversion factor from ft2 to m2
  const Real _to_m2;
  /// Conversion factor from m3 to ft3
  const Real _to_ft3;
  /// Conversion factor from ft3 to m3
  const Real _to_m3;
  /// Conversion factor from kg to lb
  const Real _to_lb;
  /// Conversion factor from lb to kg
  const Real _to_kg;
  /// Conversion factor from Pa to atm
  const Real _to_atm;
  /// Conversion factor from atm to Pa
  const Real _to_Pa;
  /// Conversion factor from K to R
  const Real _to_R;
  /// Conversion factor from R to K
  const Real _to_K;
  /// Conversion factor from J to Btu
  const Real _to_Btu;
  /// Conversion factor from Btu to J
  const Real _to_J;
  /// Conversion factor from hr to s
  const Real _to_s;
  /// Conversion factor from m3/kg to ft3/lb
  const Real _to_ft3_lb;
  /// Conversion factor from ft3/lb to m3/kg
  const Real _to_m3_kg;
  /// Conversion factor from J/kg to Btu/lb
  const Real _to_Btu_lb;
  /// Conversion factor from Btu/lb to J/kg
  const Real _to_J_kg;
  /// Conversion factor from J/(kg K) to Btu/(lb R)
  const Real _to_Btu_lbR;
  /// Conversion factor from Btu/(lb R) to J/(kg K)
  const Real _to_J_kgK;
};

#endif /* POTASSIUMLIQUIDFLUIDPROPERTIES_H */
