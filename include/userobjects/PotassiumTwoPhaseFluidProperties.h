#pragma once

#include "TwoPhaseFluidProperties.h"
#include "NaNInterface.h"

class PotassiumTwoPhaseFluidProperties;
class SinglePhaseFluidProperties;

template <>
InputParameters validParams<PotassiumTwoPhaseFluidProperties>();

/**
 * Two-phase potassium fluid properties
 *
 * Range of validity:
 *   1.01325 Pa (1e-5 atm) <= p <= 101.325 MPa (1000 atm)
 *   250 K <= T <= 2300 K
 *
 */
class PotassiumTwoPhaseFluidProperties : public TwoPhaseFluidProperties, public NaNInterface
{
public:
  PotassiumTwoPhaseFluidProperties(const InputParameters & parameters);

  virtual Real p_critical() const override;
  virtual Real T_sat(Real pressure) const override;
  virtual Real p_sat(Real temperature) const override;
  virtual Real dT_sat_dp(Real pressure) const override;
  virtual Real sigma_from_T(Real T) const override;
  virtual Real dsigma_dT_from_T(Real T) const override;

  virtual bool supportsPhaseChange() const override { return true; }

protected:
  // Critical pressure
  static const Real _P_critical;

protected:
  /// Conversion factor from Pa to atm
  const Real _to_atm;
  /// Conversion factor from atm to Pa
  const Real _to_Pa;
  /// Conversion factor from K to R
  const Real _to_R;
  /// Conversion factor from R to K
  const Real _to_K;
};
