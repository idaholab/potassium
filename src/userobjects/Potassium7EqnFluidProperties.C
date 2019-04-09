#include "Potassium7EqnFluidProperties.h"
#include "SinglePhaseFluidProperties.h"

void DIFF_ps_t_K(double t, double & ps, double & dpsdt, double & d2psdt2);
int DIFF_ts_p_K(double p, double & ts, double & dtsdp, double & d2tsdp2);
void sigma_t_K(double t, double & sigma, double & dsigmadt);

const Real Potassium7EqnFluidProperties::_P_critical = 15.95E+6;

registerMooseObject("PotassiumApp", Potassium7EqnFluidProperties);

template <>
InputParameters
validParams<Potassium7EqnFluidProperties>()
{
  InputParameters params = validParams<TwoPhaseFluidProperties>();
  params += validParams<NaNInterface>();
  params.addClassDescription("Fluid properties of potassium for the 7-equation model.");
  return params;
}

Potassium7EqnFluidProperties::Potassium7EqnFluidProperties(const InputParameters & parameters)
  : TwoPhaseFluidProperties(parameters),
    NaNInterface(this),
    _to_atm(1. / 101325.),
    _to_Pa(101325.),
    _to_R(9. / 5.),
    _to_K(5. / 9.)
{
  {
    std::string class_name = "PotassiumLiquidFluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    _fe_problem.addUserObject(class_name, _liquid_name, params);
  }
  _fp_liquid = &_fe_problem.getUserObject<SinglePhaseFluidProperties>(_liquid_name);

  {
    std::string class_name = "PotassiumVaporFluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    _fe_problem.addUserObject(class_name, _vapor_name, params);
  }
  _fp_vapor = &_fe_problem.getUserObject<SinglePhaseFluidProperties>(_vapor_name);
}

Real
Potassium7EqnFluidProperties::p_critical() const
{
  return _P_critical;
}

Real
Potassium7EqnFluidProperties::T_sat(Real pressure) const
{
  pressure *= _to_atm;

  static const double p0 = 1.E-8 * 0.101325; // just a lower limit, not actual triple point press.
  static const double pc = 15.95 * 0.101325;

  if (p0 <= pressure && pressure <= pc)
  {
    double ts, dtsdp, d2tsdp2;
    DIFF_ts_p_K(pressure, ts, dtsdp, d2tsdp2);
    return ts * _to_K;
  }
  else
    return getNaN();
}

Real
Potassium7EqnFluidProperties::p_sat(Real temperature) const
{
  temperature *= _to_R;

  static const double t0 = 273.15 * (9. / 5.); // just a lower limit, not actual triple point temp.
  static const double tc = 2239. * (9. / 5.);

  if (t0 < temperature && temperature < tc)
  {
    double ps, dpsdt, d2psdt2;
    DIFF_ps_t_K(temperature, ps, dpsdt, d2psdt2);
    return ps * _to_Pa;
  }
  else
    return getNaN();
}

Real
Potassium7EqnFluidProperties::dT_sat_dp(Real pressure) const
{
  pressure *= _to_atm;

  static const double p0 = 1.E-8 * 0.101325; // just a lower limit, not actual triple point press.
  static const double pc = 15.95 * 0.101325;

  if (p0 <= pressure && pressure <= pc)
  {
    double ts, dtsdp, d2tsdp2;
    DIFF_ts_p_K(pressure, ts, dtsdp, d2tsdp2);
    return dtsdp * _to_K / _to_Pa;
  }
  else
    return getNaN();
}

Real
Potassium7EqnFluidProperties::sigma_from_T(Real T) const
{
  double sigma, dsigmadt;
  sigma_t_K(T, sigma, dsigmadt);
  return sigma * 1e-3;
}

Real
Potassium7EqnFluidProperties::dsigma_dT_from_T(Real T) const
{
  double sigma, dsigmadt;
  sigma_t_K(T, sigma, dsigmadt);
  return dsigmadt * 1e-3;
}
