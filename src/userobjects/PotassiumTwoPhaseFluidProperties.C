#include "PotassiumTwoPhaseFluidProperties.h"
#include "SinglePhaseFluidProperties.h"

void DIFF_ps_t_K(double t, double & ps, double & dpsdt, double & d2psdt2);
int DIFF_ts_p_K(double p, double & ts, double & dtsdp, double & d2tsdp2);
void sigma_t_K(double t, double & sigma, double & dsigmadt);

const Real PotassiumTwoPhaseFluidProperties::_P_critical = 15.95E+6;

// Value is taken from NIST Chemistry WebBook, SRD 69.
// Original reference:
// R. E. Honing and D. A. Kramer. Vapor pressure data for the solid and liquid elements (1969)
const Real PotassiumTwoPhaseFluidProperties::_T_triple = 336.35;

// Value is taken from the following reference:
//
// O. J. Foust. Sodium-NaK Engineering Handbook, Volume 1: Sodium Chemistry and Physical Properties
// (1972). Division of Reactor Development and Technology, United States Atomic Energy Commission.
//
// Value given was 14.2 cal/g, which was converted to J/kg.
const Real PotassiumTwoPhaseFluidProperties::_L_fusion = 59412.8;

registerMooseObject("PotassiumApp", PotassiumTwoPhaseFluidProperties);

InputParameters
PotassiumTwoPhaseFluidProperties::validParams()
{
  InputParameters params = TwoPhaseFluidProperties::validParams();
  params += NaNInterface::validParams();
  params.addClassDescription("Two-phase potassium fluid properties");
  return params;
}

PotassiumTwoPhaseFluidProperties::PotassiumTwoPhaseFluidProperties(
    const InputParameters & parameters)
  : TwoPhaseFluidProperties(parameters),
    NaNInterface(this),
    _to_atm(1. / 101325.),
    _to_Pa(101325.),
    _to_R(9. / 5.),
    _to_K(5. / 9.)
{
  if (_tid == 0)
  {
    std::string class_name = "PotassiumLiquidFluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    params.set<MooseEnum>("emit_on_nan") = getParam<MooseEnum>("emit_on_nan");
    _fe_problem.addUserObject(class_name, _liquid_name, params);
  }
  _fp_liquid = &_fe_problem.getUserObject<SinglePhaseFluidProperties>(_liquid_name, _tid);

  if (_tid == 0)
  {
    std::string class_name = "PotassiumVaporFluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    params.set<MooseEnum>("emit_on_nan") = getParam<MooseEnum>("emit_on_nan");
    _fe_problem.addUserObject(class_name, _vapor_name, params);
  }
  _fp_vapor = &_fe_problem.getUserObject<SinglePhaseFluidProperties>(_vapor_name, _tid);
}

Real
PotassiumTwoPhaseFluidProperties::p_critical() const
{
  return _P_critical;
}

Real
PotassiumTwoPhaseFluidProperties::T_triple() const
{
  return _T_triple;
}

Real
PotassiumTwoPhaseFluidProperties::L_fusion() const
{
  return _L_fusion;
}

Real
PotassiumTwoPhaseFluidProperties::T_sat(Real pressure) const
{
  pressure *= _to_atm;

  static const double p0 = 1.E-8 * 0.101325; // just a lower limit, not actual triple point press.
  static const double pc = 15.95 * 0.101325;

  if (p0 <= pressure && pressure <= pc)
  {
    double ts, dtsdp, d2tsdp2;
    int ierr = DIFF_ts_p_K(pressure, ts, dtsdp, d2tsdp2);
    if (ierr != 0)
      return getNaN();
    else
      return ts * _to_K;
  }
  else
    return getNaN();
}

Real
PotassiumTwoPhaseFluidProperties::p_sat(Real temperature) const
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
PotassiumTwoPhaseFluidProperties::dT_sat_dp(Real pressure) const
{
  pressure *= _to_atm;

  static const double p0 = 1.E-8 * 0.101325; // just a lower limit, not actual triple point press.
  static const double pc = 15.95 * 0.101325;

  if (p0 <= pressure && pressure <= pc)
  {
    double ts, dtsdp, d2tsdp2;
    int ierr = DIFF_ts_p_K(pressure, ts, dtsdp, d2tsdp2);
    if (ierr != 0)
      return getNaN();
    else
      return dtsdp * _to_K / _to_Pa;
  }
  else
    return getNaN();
}

Real
PotassiumTwoPhaseFluidProperties::sigma_from_T(Real T) const
{
  double sigma, dsigmadt;
  sigma_t_K(T, sigma, dsigmadt);
  return sigma * 1e-3;
}

Real
PotassiumTwoPhaseFluidProperties::dsigma_dT_from_T(Real T) const
{
  double sigma, dsigmadt;
  sigma_t_K(T, sigma, dsigmadt);
  return dsigmadt * 1e-3;
}
