#include "PotassiumVaporFluidPropertiesTest.h"
#include "SinglePhaseFluidPropertiesTestUtils.h"

TEST_F(PotassiumVaporFluidPropertiesTest, test)
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

  // rho and v
  REL_TEST(rho_from_p_T, 7.6680283349292563, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->rho_from_p_s(p, s), rho_from_p_T, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->rho_from_p_T, p, T, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->rho_from_p_s, p, s, REL_TOL_DERIVATIVE);

  // e
  REL_TEST(e_from_p_rho, 2840324.0011762767, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->e_from_v_h(v, h), e, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->e_from_p_rho, p, rho, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->e_from_v_h, v, h, REL_TOL_DERIVATIVE);

  // c
  const Real c = _fp->c_from_v_e(v, e);
  REL_TEST(c, 608.72231247088769, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->c_from_v_e, v, e, REL_TOL_DERIVATIVE);

  // cp
  const Real cp = _fp->cp_from_v_e(v, e);
  REL_TEST(cp, 1316.6708120023293, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->cp_from_v_e, v, e, REL_TOL_DERIVATIVE);

  // cv
  const Real cv = _fp->cv_from_v_e(v, e);
  REL_TEST(cv, 802.64776711320167, REL_TOL_SAVED_VALUE);

  // mu
  const Real mu = _fp->mu_from_v_e(v, e);
  REL_TEST(mu, 0.000021913164388771186, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->mu_from_v_e, v, e, REL_TOL_DERIVATIVE);

  // k
  const Real k = _fp->k_from_v_e(v, e);
  REL_TEST(k, 0.02607040851256066, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->k_from_v_e, v, e, REL_TOL_DERIVATIVE);

  // s
  REL_TEST(s, 4186.8350242420565, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(_fp->s_from_h_p(h, p), s, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->s_from_v_e, v, e, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->s_from_h_p, h, p, REL_TOL_DERIVATIVE);

  // g
  REL_TEST(_fp->g_from_v_e(v, e), -3183545.9787879535, REL_TOL_EXTERNAL_VALUE);

  // h
  REL_TEST(h_from_p_T, 3096706.5575751313, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->h_from_p_T, p, T, REL_TOL_DERIVATIVE);

  // beta
  Real beta = _fp->beta_from_p_T(p, T);
  REL_TEST(beta, 0.0012801359946759459, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->beta_from_p_T, p, T, REL_TOL_DERIVATIVE);

  // critical parameters
  REL_TEST(_fp->criticalTemperature(), 2239., REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->criticalDensity(), 192., REL_TOL_SAVED_VALUE);
}
