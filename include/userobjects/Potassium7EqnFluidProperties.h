#ifndef POTASSIUM7EQNFLUIDPROPERTIES_H
#define POTASSIUM7EQNFLUIDPROPERTIES_H

#include "TwoPhaseFluidProperties.h"
#include "NaNInterface.h"

class Potassium7EqnFluidProperties;
class SinglePhaseFluidProperties;

template <>
InputParameters validParams<Potassium7EqnFluidProperties>();

/**
 * Potassium interface for 7-eqn model
 */
class Potassium7EqnFluidProperties : public TwoPhaseFluidProperties, public NaNInterface
{
public:
  Potassium7EqnFluidProperties(const InputParameters & parameters);

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

#endif /* POTASSIUM7EQNFLUIDPROPERTIES_H */
