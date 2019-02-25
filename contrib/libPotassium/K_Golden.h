#pragma once
//
void DIFF_ps_t_K(double t, double & ps, double & dpsdt, double & d2psdt2);
int DIFF_ts_p_K(double p, double & ts, double & dtsdp, double & d2tsdp2);
void DIFF_d1_t_K(double t, double & d1, double & dd1dt, double & d2d1dt2);
void DIFF_d_tp_L_K(double t,
                    double p,
                    double & dl,
                    double & ddldt,
                    double & d2dldt2,
                    double & ddldp,
                    double & d2dldp2,
                    double & d2dldtdp);
void DIFF_d_tp_L_K(double t,
                    double p,
                    double & dl,
                    double & ddldt,
                    double & d2dldt2,
                    double & d3dldt3,
                    double & ddldp,
                    double & d2dldp2,
                    double & d2dldtdp,
                    double & d3dldt2dp,
                    double & d3dldtdp2);
void DIFF_v1_t_K(double t, double & v1, double & dv1dt, double & d2v1dt2);
void DIFF_v_tp_L_K(double t,
                    double p,
                    double & vl,
                    double & dvldt,
                    double & d2vldt2,
                    double & dvldp,
                    double & d2vldp2,
                    double & d2vldtdp);
void DIFF_h1_t_K(double t, double & h1, double & dh1dt, double & d2h1dt2);
void DIFF_vu_tp_L_K(double t,
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
                     double & d2uldtdp);
void DIFF_h_tp_L_K(double t,
                    double p,
                    double & hl,
                    double & dhldt,
                    double & d2hldt2,
                    double & dhldp,
                    double & d2hldp2,
                    double & d2hldtdp);
void DIFF_s1_t_K(double t, double & s1, double & ds1dt, double & d2s1dt2);
void DIFF_s_tp_L_K(double t,
                    double p,
                    double & sl,
                    double & dsldt,
                    double & d2sldt2,
                    double & dsldp,
                    double & d2sldp2,
                    double & d2sldtdp);
void DIFF_qc_t_K(double t,
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
                  double & d2b4dt2);
void DIFF_qc_tp_K(double t,
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
                   double & d2b4dtdp);
void DIFF_v2_t_K(double t, double & v2, double & dv2dt, double & d2v2dt2);
void DIFF_v_tp_G_K(double t,
                    double p,
                    double & vv,
                    double & dvvdt,
                    double & d2vvdt2,
                    double & dvvdp,
                    double & d2vvdp2,
                    double & d2vvdtdp);
void DIFF_dhv_t_K(double t, double & dhv, double & ddhvdt, double & d2dhvdt2);
void DIFF_h2_t_K(double t, double & h2, double & dh2dt, double & d2h2dt2);
void DIFF_vu_tp_G_K(double t,
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
                     double & d2uvdtdp);
void DIFF_h_tp_G_K(double t,
                    double p,
                    double & hv,
                    double & dhvdt,
                    double & d2hvdt2,
                    double & dhvdp,
                    double & d2hvdp2,
                    double & d2hvdtdp);
void DIFF_s2_t_K(double t, double & s2, double & ds2dt, double & d2s2dt2);
void DIFF_s_tp_G_K(double t,
                    double p,
                    double & sv,
                    double & dsvdt,
                    double & d2svdt2,
                    double & dsvdp,
                    double & d2svdp2,
                    double & d2svdtdp);
// derived properties
void DIFF_cp_tp_L_K(double t, double p, double & cp, double & dcpdt, double & dcpdp);
void DIFF_cp_tp_G_K(double t, double p, double & cp, double & dcpdt, double & dcpdp);
void DIFF_cv_tp_L_K(double t, double p, double & cv, double & dcvdt, double & dcvdp);
void DIFF_cv_tp_G_K(double t, double p, double & cv, double & dcvdt, double & dcvdp);
void DIFF_w_tp_L_K(double t, double p, double & w, double & dwdt, double & dwdp);
void DIFF_w_tp_G_K(double t, double p, double & w, double & dwdt, double & dwdp);
// flash routines
int FLASH_prho_L_K(double p, double rho, double & t);
void DERIV_prho_L_K(double t, double p, double & du_dp_rho, double & du_drho_p);
int FLASH_prho_G_K(double p, double rho, double & t);
void DERIV_prho_G_K(double t, double p, double & du_dp_rho, double & du_drho_p);
int FLASH_ph_L_K(double p, double h, double & t);
int FLASH_ph_G_K(double p, double h, double & t);
int FLASH_ps_L_K(double p, double s, double & t);
int FLASH_ps_G_K(double p, double s, double & t);
int FLASH_vu_L_K(double v, double u, double & t, double & p);
void DERIV_vu_L_K(
    double t, double p, double & dp_dv_u, double & dp_du_v, double & dt_dv_u, double & dt_du_v);
int FLASH_vu_G_K(double v, double u, double & t, double & p);
void DERIV_vu_G_K(
    double t, double p, double & dp_dv_u, double & dp_du_v, double & dt_dv_u, double & dt_du_v);
int FLASH_vh_L_K(double v, double h, double & t, double & p);
void DERIV_vh_L_K(
    double t, double p, double & dp_dv_h, double & dp_dh_v, double & dt_dv_h, double & dt_dh_v);
int FLASH_vh_G_K(double v, double h, double & t, double & p);
void DERIV_vh_G_K(
    double t, double p, double & dp_dv_h, double & dp_dh_v, double & dt_dv_h, double & dt_dh_v);
int FLASH_hs_L_K(double h, double s, double & t, double & p);
int FLASH_hs_G_K(double h, double s, double & t, double & p);
// transport properties
void sigma_t_K(double t, double & sigma, double & dsigmadt); // K, mN/m
double etal_tv_K(double t, double v);
void etav_t_K(double t, double & etav, double & detavdt, double & d2etavdt2);
void lambdal_t_K(double t,
                  double & lambdal,
                  double & dlambdaldt,
                  double & d2lambdaldt2); // R, Btu/hr-ft-F
void lambdav_t_K(double t,
                  double & lambdav,
                  double & dlambdavdt,
                  double & d2lambdavdt2); // R, Btu/hr-ft-F

void
betat_t_K(double t, double & betat, double & dbetatdt, double & d2betatdt2, double & d3betatdt3);
