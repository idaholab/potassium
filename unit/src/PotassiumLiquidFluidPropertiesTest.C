#include "PotassiumLiquidFluidPropertiesTest.h"
#include "SinglePhaseFluidPropertiesTestUtils.h"

TEST_F(PotassiumLiquidFluidPropertiesTest, test)
{
  const Real T = 1500.;
  const Real p = 1965948.7070480157;

  const Real rho_from_p_T = _fp->rho_from_p_T(p, T);
  const Real rho = rho_from_p_T;

  const Real h_from_p_T = _fp->h_from_p_T(p, T);
  const Real h = h_from_p_T;

  const Real e_from_p_rho = _fp->e_from_p_rho(p, rho);
  const Real e = e_from_p_rho;

  const Real v = 1 / rho;

  const Real s_from_v_e = _fp->s_from_v_e(v, e);
  const Real s = s_from_v_e;

  // p
  REL_TEST(_fp->p_from_v_e(v, e), p, REL_TOL_CONSISTENCY);
  REL_TEST(_fp->p_from_h_s(h, s), p, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->p_from_v_e, v, e, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->p_from_h_s, h, s, REL_TOL_DERIVATIVE);

  // T
  REL_TEST(_fp->T_from_v_e(v, e), T, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->T_from_v_e, v, e, REL_TOL_DERIVATIVE);

  // rho
  REL_TEST(rho_from_p_T, 542.92729655194148, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->rho_from_p_s(p, s), rho_from_p_T, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->rho_from_p_T, p, T, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->rho_from_p_s, p, s, REL_TOL_DERIVATIVE);

  // e
  REL_TEST(e_from_p_rho, 1455900.029063388, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->e_from_v_h(v, h), e, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->e_from_p_rho, p, rho, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->e_from_v_h, v, h, REL_TOL_DERIVATIVE);

  // c
  const Real c = _fp->c_from_v_e(v, e);
  REL_TEST(c, 1176.5347212266745, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->c_from_v_e, v, e, REL_TOL_DERIVATIVE);

  // cp
  const Real cp = _fp->cp_from_v_e(v, e);
  REL_TEST(cp, 1059.3392531342224, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->cp_from_v_e, v, e, REL_TOL_DERIVATIVE);

  // cv
  const Real cv = _fp->cv_from_v_e(v, e);
  REL_TEST(cv, 724.81862372294802, REL_TOL_SAVED_VALUE);

  // mu
  const Real mu = _fp->mu_from_v_e(v, e);
  REL_TEST(mu, 0.000092867260153346677, REL_TOL_SAVED_VALUE);

  // k
  const Real k = _fp->k_from_v_e(v, e);
  REL_TEST(k, 19.197175558481387, REL_TOL_SAVED_VALUE);

  // s
  REL_TEST(s, 3095.3780162439684, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(_fp->s_from_h_p(h, p), s, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->s_from_v_e, v, e, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->s_from_h_p, h, p, REL_TOL_DERIVATIVE);

  // g
  REL_TEST(_fp->g_from_v_e(v, e), -3183545.978787953, REL_TOL_EXTERNAL_VALUE);

  // h
  REL_TEST(h_from_p_T, 1459521.0455779997, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->h_from_p_T, p, T, REL_TOL_DERIVATIVE);

  // beta
  const Real beta = _fp->beta_from_p_T(p, T);
  REL_TEST(beta, 0.00050075752661458332, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->beta_from_p_T, p, T, REL_TOL_DERIVATIVE);
}
