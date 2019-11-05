#include "PotassiumTwoPhaseFluidPropertiesTest.h"
#include "SinglePhaseFluidPropertiesTestUtils.h"

TEST_F(PotassiumTwoPhaseFluidPropertiesTest, test)
{
  const Real relative_perturbation = 1e-6;

  Real T = 1200.;  // K
  Real p = 101325; // Pa

  // Tsat + derivatives
  REL_TEST(_fp->T_sat(p), 1029.7303271118376, REL_TOL_SAVED_VALUE);
  {
    Real dT_dPsat = _fp->dT_sat_dp(p);

    Real dp = relative_perturbation * p;
    Real dT_dPsat_fd = (_fp->T_sat(p + dp) - _fp->T_sat(p - dp)) / (2 * dp);

    REL_TEST(dT_dPsat, dT_dPsat_fd, REL_TOL_DERIVATIVE);
  }

  // Psat
  REL_TEST(_fp->p_sat(T), 391352.20564277115, REL_TOL_SAVED_VALUE);

  {
    // sigma
    T = 1500.;
    REL_TEST(_fp->sigma_from_T(T), 0.037181600000000002, REL_TOL_SAVED_VALUE);
    DERIV_TEST_1D(_fp->sigma_from_T, _fp->dsigma_dT_from_T, T, REL_TOL_DERIVATIVE);
  }
}
