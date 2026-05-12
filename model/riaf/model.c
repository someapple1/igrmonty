#include "decs.h"
#include "coordinates.h"
#include "model_radiation.h"

#include <assert.h>

static double Rmax_record = 1.e4;

double rmax_geo = 1000.;
double rmin_geo = 0.;
double MBH_solar = 4.3e6;

double Te_unit = 1.e11;
double Ne_unit = 5.e6;

double nth0, Te0, disk_h, pow_nth, pow_T;
double keplerian_factor, infall_factor;
double r_isco;

void report_bad_input(int argc)
{
  if (argc < 2) {
    fprintf(stderr, "usage: \n");
    fprintf(stderr, "  riaf:    Ns\n");
    exit(0);
  }
}

///////////////////////////////// SUPERPHOTONS /////////////////////////////////

#define ROULETTE  1.e4
int stop_criterion(struct of_photon *ph)
{
  // stop if weight below minimum weight
  if (ph->w < WEIGHT_MIN) {
    if (monty_rand() <= 1. / ROULETTE) {
      ph->w *= ROULETTE;
    } else {
      ph->w = 0.;
      return 1;
    }
  }

  double r, h;
  bl_coord(ph->X, &r, &h);

  if (r < Rh*1.05 || r > Rmax_record) {
    return 1;
  }

  return 0;
}

int record_criterion(struct of_photon *ph)
{
  double r, h;
  bl_coord(ph->X, &r, &h);

  if (r >= Rmax_record) {
    return 1;
  }

  return 0;
}
#undef ROULETTE

#define EPS 0.01
double stepsize(double X[NDIM], double K[NDIM])
{
  double dl, dlx1, dlx2, dlx3;
  double idlx1, idlx2, idlx3;

  #define MIN(A,B) (A<B?A:B)

  dlx1 = EPS / (fabs(K[1]) + SMALL);
  dlx2 = EPS * MIN(X[2], 1. - X[2]) / (fabs(K[2]) + SMALL);
  dlx3 = EPS / (fabs(K[3]) + SMALL) ;

  #undef MIN

  idlx1 = 1. / (fabs(dlx1) + SMALL);
  idlx2 = 1. / (fabs(dlx2) + SMALL);
  idlx3 = 1. / (fabs(dlx3) + SMALL);

  dl = 1. / (idlx1 + idlx2 + idlx3);

  return dl;
}
#undef EPS

void record_super_photon(struct of_photon *ph)
{
  double lE, dx2;
  int iE, ix2, ic;

  if (isnan(ph->w) || isnan(ph->E)) {
    fprintf(stderr, "record isnan: %g %g\n", ph->w, ph->E);
    return;
  }

  // bin in X[2] BL coord while folding around the equator and check limit
  double r, th;
  bl_coord(ph->X, &r, &th);
  dx2 = M_PI/2./N_THBINS;
  if (th > M_PI/2.) {
    ix2 = (int)( (M_PI - th) / dx2 );
  } else {
    ix2 = (int)( th / dx2 );
  }
  if (ix2 < 0 || ix2 >= N_THBINS) return;

  #if CUSTOM_AVG==1
  double nu = ph->E * ME*CL*CL / HPL;
  if (nu < CA_MIN_FREQ || CA_MAX_FREQ < nu) return;
  // Get custom average bin
  dlE = (log(CA_MAX_QTY) - log(CA_MIN_QTY)) / CA_NBINS;
  lE = log(ph->QTY0);
  iE = (int) ((lE - log(CA_MIN_QTY)) / dlE + 2.5) - 2;
  if (iE < 0 || iE >= CA_NBINS) return;
  #else
  // Get energy bin (centered on iE*dlE + lE0)
  lE = log(ph->E);
  iE = (int) ((lE - lE0) / dlE + 2.5) - 2;
  if (iE < 0 || iE >= N_EBINS) return;
  #endif // CUSTOM_AVG

  // Get compton bin
  ic = ph->nscatt;
  if (ic > 3) ic = 3;

  #pragma omp atomic
  N_superph_recorded++;

  double ratio_synch = 1. - ph->ratio_brems;

  // Add superphoton to synch spectrum
  spect[ic][ix2][iE].dNdlE += ph->w * ratio_synch;
  spect[ic][ix2][iE].dEdlE += ph->w * ph->E * ratio_synch;
  spect[ic][ix2][iE].tau_abs += ph->w * ph->tau_abs * ratio_synch;
  spect[ic][ix2][iE].tau_scatt += ph->w * ph->tau_scatt * ratio_synch;
  spect[ic][ix2][iE].X1iav += ph->w * ph->X1i * ratio_synch;
  spect[ic][ix2][iE].X2isq += ph->w * (ph->X2i * ph->X2i) * ratio_synch;
  spect[ic][ix2][iE].X3fsq += ph->w * (ph->X[3] * ph->X[3]) * ratio_synch;
  spect[ic][ix2][iE].ne0 += ph->w * (ph->ne0) * ratio_synch;
  spect[ic][ix2][iE].b0 += ph->w * (ph->b0) * ratio_synch;
  spect[ic][ix2][iE].thetae0 += ph->w * (ph->thetae0) * ratio_synch;
  spect[ic][ix2][iE].nscatt += ph->w * ph->nscatt * ratio_synch;
  spect[ic][ix2][iE].nph += 1.;
  // .. to brems spectrum
  spect[ic+(N_COMPTBINS+1)][ix2][iE].dNdlE += ph->w * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].dEdlE += ph->w * ph->E * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].tau_abs += ph->w * ph->tau_abs * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].tau_scatt += ph->w * ph->tau_scatt * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].X1iav += ph->w * ph->X1i * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].X2isq += ph->w * (ph->X2i * ph->X2i) * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].X3fsq += ph->w * (ph->X[3] * ph->X[3]) * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].ne0 += ph->w * (ph->ne0) * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].b0 += ph->w * (ph->b0) * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].thetae0 += ph->w * (ph->thetae0) * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].nscatt += ph->w * ph->nscatt * ph->ratio_brems;
  spect[ic+(N_COMPTBINS+1)][ix2][iE].nph += 1.;
}

struct of_spectrum shared_spect[N_TYPEBINS][N_THBINS][N_EBINS] = { };

void omp_reduce_spect()
{
  // TODO: should shared_spect be explicitly set to zero (in all model files?)
  #pragma omp parallel
  {
    #pragma omp critical
    {
      for (int k = 0; k < N_TYPEBINS; k++) {
        for (int i = 0; i < N_THBINS; i++) {
          for (int j = 0; j < N_EBINS; j++) {
            shared_spect[k][i][j].dNdlE +=
                spect[k][i][j].dNdlE;
            shared_spect[k][i][j].dEdlE +=
                spect[k][i][j].dEdlE;
            shared_spect[k][i][j].tau_abs +=
                spect[k][i][j].tau_abs;
            shared_spect[k][i][j].tau_scatt +=
                spect[k][i][j].tau_scatt;
            shared_spect[k][i][j].X1iav +=
                spect[k][i][j].X1iav;
            shared_spect[k][i][j].X2isq +=
                spect[k][i][j].X2isq;
            shared_spect[k][i][j].X3fsq +=
                spect[k][i][j].X3fsq;
            shared_spect[k][i][j].ne0 +=
                spect[k][i][j].ne0;
            shared_spect[k][i][j].b0 +=
                spect[k][i][j].b0;
            shared_spect[k][i][j].thetae0 +=
                spect[k][i][j].thetae0;
            shared_spect[k][i][j].nscatt +=
                spect[k][i][j].nscatt;
            shared_spect[k][i][j].nph +=
                spect[k][i][j].nph;
          }
        }
      }
    } // omp critical

    #pragma omp barrier

    #pragma omp master
    {
      for (int k = 0; k < N_TYPEBINS; k++) {
        for (int i = 0; i < N_THBINS; i++) {
          for (int j = 0; j < N_EBINS; j++) {
            spect[k][i][j].dNdlE =
                shared_spect[k][i][j].dNdlE;
            spect[k][i][j].dEdlE =
                shared_spect[k][i][j].dEdlE;
            spect[k][i][j].tau_abs =
                shared_spect[k][i][j].tau_abs;
            spect[k][i][j].tau_scatt =
                shared_spect[k][i][j].tau_scatt;
            spect[k][i][j].X1iav =
                shared_spect[k][i][j].X1iav;
            spect[k][i][j].X2isq =
                shared_spect[k][i][j].X2isq;
            spect[k][i][j].X3fsq =
                shared_spect[k][i][j].X3fsq;
            spect[k][i][j].ne0 =
                shared_spect[k][i][j].ne0;
            spect[k][i][j].b0 =
                shared_spect[k][i][j].b0;
            spect[k][i][j].thetae0 =
                shared_spect[k][i][j].thetae0;
            spect[k][i][j].nscatt =
                shared_spect[k][i][j].nscatt;
            spect[k][i][j].nph =
                shared_spect[k][i][j].nph;
          }
        }
      }
    } // omp master
  } // omp parallel
}

double bias_func(double Te, double w)
{
  double bias, max;

  max = 0.5 * w / WEIGHT_MIN;

  bias = Te * Te / (5. * max_tau_scatt);
  // bias = 100. * Te * Te / (bias_norm * max_tau_scatt);

  if (bias > max)
    bias = max;

  return  bias * biasTuning;


  // TODO maybe swap this out with something in sphere_old or simplesphere ?

  // use old method with bias tuning parameter ?
  /*
  double bias, max;

  max = 0.5 * w / WEIGHT_MIN;

  if (Te > SCATTERING_THETAE_MAX) Te = SCATTERING_THETAE_MAX;
  bias = 16. * Te * Te / (5. * max_tau_scatt);

  if (bias > max) bias = max;

  return bias * biasTuning;
   */
}

void get_fluid_zone(int i, int j, int k, double *Ne, double *Thetae, double *B,
        double Ucon[NDIM], double Bcon[NDIM])
{
  // get grid zone center
  double X[4] = { 0. };
  ijktoX(i, j, k, X);

  // fill gcov
  double gcov[4][4];
  gcov_func(X, gcov);

  // create dummy vectors
  double Ucov[4] = { 0. };
  double Bcov[4] = { 0. };

  // call into analytic model
  get_fluid_params(X, gcov, Ne, Thetae, B, Ucon, Ucov, Bcon, Bcov);
}

double _get_model_Ne(double r, double th)
{
  double zc = r * cos(th);
  double rc = r * sin(th);
  return nth0 * exp(-zc * zc / 2. / rc / rc / disk_h / disk_h) *
         pow(r, pow_nth) * Ne_unit;
}

static double _get_model_Thetae(double r)
{
  return Te0 * pow(r, pow_T) * Te_unit * KBOL / (ME * CL * CL);
}

double _get_model_Bmag(double r, double th, double Ne)
{
  (void)th;
  if (Ne <= 0. || r <= 0.) {
    return 0.;
  }

  return sqrt(8. * M_PI) * 0.2 * Ne * MP * CL * CL / (6. * r);
}

double get_model_sigma(const double X[NDIM])
{
  double r, th;
  bl_coord(X, &r, &th);

  double Ne = _get_model_Ne(r, th);
  double Bmag = _get_model_Bmag(r, th, Ne);
  if (Ne <= 0. || Bmag <= 0.) {
    return 1.e100;
  }

  return Bmag * Bmag / (4. * M_PI * Ne * (MP + ME) * CL * CL);
}

double get_model_beta(const double X[NDIM])
{
  double r, th;
  bl_coord(X, &r, &th);

  double Ne = _get_model_Ne(r, th);
  double Thetae = _get_model_Thetae(r);
  double Bmag = _get_model_Bmag(r, th, Ne);
  if (Ne <= 0. || Thetae <= 0. || Bmag <= 0.) {
    return 0.;
  }

  return 8. * M_PI * Ne * Thetae * ME * CL * CL / (Bmag * Bmag);
}

static void zero_plasma(double *Ne, double *Thetae, double *B,
          double Ucon[NDIM], double Ucov[NDIM],
          double Bcon[NDIM], double Bcov[NDIM])
{
  *Ne = 0.;
  *Thetae = 0.;
  *B = 0.;
  MULOOP {
    Ucon[mu] = 0.;
    Ucov[mu] = 0.;
    Bcon[mu] = 0.;
    Bcov[mu] = 0.;
  }
}

static int hjw_bl_four_velocity(double r, double th,
          double bl_gcov[NDIM][NDIM], double bl_gcon[NDIM][NDIM],
          double bl_Ucon[NDIM])
{
  double R = r * sin(th);
  if (!isfinite(R) || R <= 0.) {
    return 0;
  }

  double l = R * sqrt(R) / (R + 1.);
  double circ_Ucov[NDIM] = { -1., 0., 0., l };
  double circ_Ucon[NDIM] = { 0. };

  MUNULOOP {
    circ_Ucon[mu] += bl_gcon[mu][nu] * circ_Ucov[nu];
  }

  double circ_norm = 0.;
  MULOOP {
    circ_norm += circ_Ucon[mu] * circ_Ucov[mu];
  }
  if (!isfinite(circ_norm) || circ_norm >= 0.) {
    return 0;
  }

  double circ_scale = sqrt(-circ_norm);
  MULOOP {
    circ_Ucon[mu] /= circ_scale;
  }

  double infall_radicand = -(1. + bl_gcon[0][0]) * bl_gcon[1][1];
  if (!isfinite(infall_radicand) || infall_radicand < 0.) {
    return 0;
  }

  double infall_Ucon[NDIM] = {
    -bl_gcon[0][0],
    sqrt(infall_radicand),
    0.,
    -bl_gcon[0][3]
  };

  if (circ_Ucon[0] == 0. || infall_Ucon[0] == 0.) {
    return 0;
  }

  double omega_orbital = circ_Ucon[3] / circ_Ucon[0];
  double omega_infall = infall_Ucon[3] / infall_Ucon[0];
  double omega = omega_orbital + keplerian_factor * (omega_infall - omega_orbital);
  double ur = circ_Ucon[1] + infall_factor * (infall_Ucon[1] - circ_Ucon[1]);

  double denom = bl_gcov[0][0] + 2. * omega * bl_gcov[0][3] +
                 omega * omega * bl_gcov[3][3];
  double ut_radicand = -(1. + bl_gcov[1][1] * ur * ur) / denom;
  if (!isfinite(ut_radicand) || ut_radicand <= 0.) {
    return 0;
  }

  bl_Ucon[0] = sqrt(ut_radicand);
  bl_Ucon[1] = ur;
  bl_Ucon[2] = 0.;
  bl_Ucon[3] = bl_Ucon[0] * omega;

  double bl_Ucov[NDIM];
  lower(bl_Ucon, bl_gcov, bl_Ucov);
  double norm = 0.;
  MULOOP {
    norm += bl_Ucon[mu] * bl_Ucov[mu];
  }

  return isfinite(norm) && fabs(norm + 1.) < 1.e-8;
}

void get_fluid_params(const double X[NDIM], double gcov[NDIM][NDIM], double *Ne,
          double *Thetae, double *B, double Ucon[NDIM],
          double Ucov[NDIM], double Bcon[NDIM],
          double Bcov[NDIM])
{
  double r, th;
  bl_coord(X, &r, &th);

  if (r < Rin || r > Rout) {
    zero_plasma(Ne, Thetae, B, Ucon, Ucov, Bcon, Bcov);
    return;
  }

  // set scalars 
  *Ne = _get_model_Ne(r, th);
  *Thetae = _get_model_Thetae(r);
  *B = _get_model_Bmag(r, th, *Ne);

  // Metrics: BL
  double bl_gcov[NDIM][NDIM], bl_gcon[NDIM][NDIM];
  gcov_bl(r, th, bl_gcov);
  gcon_func(bl_gcov, bl_gcon);

  double bl_Ucon[NDIM] = { NAN, NAN, NAN, NAN };
  int velocity_ok = hjw_bl_four_velocity(r, th, bl_gcov, bl_gcon, bl_Ucon);
  // USER PATCH: skip zones where the HJW-style four-velocity is not timelike.
  if (!velocity_ok) {
    zero_plasma(Ne, Thetae, B, Ucon, Ucov, Bcon, Bcov);
    return;
  }

  // Transform to KS coordinates,
  double ks_Ucon[NDIM];
  bl_to_ks(X, bl_Ucon, ks_Ucon);
  // then to our coordinates,
  vec_from_ks(X, ks_Ucon, Ucon);

  // and grab Ucov
  lower(Ucon, gcov, Ucov);


  // Check
#if DEBUG
  //if (r < r_isco) { fprintf(stderr, "ur = %g\n", Ucon[1]); }
  double debug_bl_Ucov[NDIM];
  double dot_U = Ucon[0]*Ucov[0] + Ucon[1]*Ucov[1] + Ucon[2]*Ucov[2] + Ucon[3]*Ucov[3];
  double sum_U = Ucon[0]+Ucon[1]+Ucon[2]+Ucon[3];
  // Following condition gets handled better above
  // (r < r_isco && fabs(Ucon[1]) < 1e-10) ||
  if (get_fluid_nu(Kcon, Ucov) == 1. ||
      fabs(fabs(dot_U) - 1.) > 1e-10 || sum_U < 0.1) {
    lower(bl_Ucon, bl_gcov, debug_bl_Ucov);
    fprintf(stderr, "RIAF model problem at r, th, phi = %g %g %g\n", r, th, X[3]);
    fprintf(stderr, "Ucon BL: %g %g %g %g\n", bl_Ucon[0], bl_Ucon[1], bl_Ucon[2], bl_Ucon[3]);
    fprintf(stderr, "Ucon KS: %g %g %g %g\n", ks_Ucon[0], ks_Ucon[1], ks_Ucon[2], ks_Ucon[3]);
    fprintf(stderr, "Ucon native: %g %g %g %g\n", Ucon[0], Ucon[1], Ucon[2], Ucon[3]);
    fprintf(stderr, "Ucov: %g %g %g %g\n", Ucov[0], Ucov[1], Ucov[2], Ucov[3]);
    fprintf(stderr, "Ubl.Ubl: %g\n", debug_bl_Ucov[0]*bl_Ucon[0]+debug_bl_Ucov[1]*bl_Ucon[1]+
                                    debug_bl_Ucov[2]*bl_Ucon[2]+debug_bl_Ucov[3]*bl_Ucon[3]);
    fprintf(stderr, "U.U: %g\n", Ucov[0]*Ucon[0]+Ucov[1]*Ucon[1]+Ucov[2]*Ucon[2]+Ucov[3]*Ucon[3]);
  }
#endif

  // HJW RIAF toroidal field.  In BL coordinates this choice is
  // orthogonal to the fluid four-velocity by construction.
  double bl_Ucov[NDIM];
  lower(bl_Ucon, bl_gcov, bl_Ucov);
  if (bl_Ucov[0] == 0.) {
    zero_plasma(Ne, Thetae, B, Ucon, Ucov, Bcon, Bcov);
    return;
  }

  double bl_Bcon[NDIM];
  bl_Bcon[0] = -bl_Ucov[3] / bl_Ucov[0];
  bl_Bcon[1] = 0.0;
  bl_Bcon[2] = 0.0;
  bl_Bcon[3] = 1.0;

  double bl_Bcov[NDIM];
  lower(bl_Bcon, bl_gcov, bl_Bcov);
  double bl_Bsq = 0.;
  MULOOP bl_Bsq += bl_Bcon[mu] * bl_Bcov[mu];
  if (!isfinite(bl_Bsq) || bl_Bsq <= 0. || *B <= 0.) {
    MULOOP {
      Bcon[mu] = 0.;
      Bcov[mu] = 0.;
    }
    *B = 0.;
    return;
  }

  double bmag_code = *B / B_unit;
  MULOOP bl_Bcon[mu] *= bmag_code / sqrt(bl_Bsq);

  // Transform to KS coordinates,
  double ks_Bcon[NDIM];
  bl_to_ks(X, bl_Bcon, ks_Bcon);
  // then to our coordinates,
  vec_from_ks(X, ks_Bcon, Bcon);
  lower(Bcon, gcov, Bcov);
}

////////////////////////////////// COORDINATES /////////////////////////////////

void gcov_func(const double X[NDIM], double gcov[NDIM][NDIM])
{
  // despite the name, get equivalent values for
  // r, th for KS coordinates
  double r, th;
  bl_coord(X, &r, &th);

  // compute ks metric
  double gcovKS[NDIM][NDIM];
  gcov_ks(r, th, gcovKS);

  // Apply coordinate transformation to code coordinates X
  double dxdX[NDIM][NDIM];
  set_dxdX(X, dxdX);

  MUNULOOP {
    gcov[mu][nu] = 0.;
    for (int lam = 0; lam < NDIM; lam++) {
      for (int kap = 0; kap < NDIM; kap++) {
        gcov[mu][nu] += gcovKS[lam][kap]*dxdX[lam][mu]*dxdX[kap][nu];
      }
    }
  }
}

double dOmega_func(int j)
{
  double dx2 = M_PI/2./N_THBINS;
  double thi = j * dx2;
  double thf = (j+1) * dx2;

  return 2.*M_PI*(-cos(thf) + cos(thi));
}

//////////////////////////////// INITIALIZATION ////////////////////////////////

#include <hdf5.h>
#include <hdf5_hl.h>

void init_data(int argc, char *argv[], Params *params)
{
  with_derefine_poles = 0; // since we've set nothing, this should default to MKS
  hslope = 1.;

  // parameter defaults
  MBH_solar = 4.3e6;
  Ne_unit = 6.e7;
  Te_unit = 1.5e11;
  //rmax_geo = ? // TODO, do these two need to be re-set if we use weird input parameters?
  //rmin_geo = ?
  a = 0.9375;
  nth0 = 1.;
  Te0 = 1.;
  disk_h = 0.5;
  pow_nth = -1.1;
  pow_T = -8.4;
  keplerian_factor = 0.5;
  infall_factor = 0.5;

  // TODO deal with this in a more clever way
  if (params->loaded && strlen(params->dump) > 0) {
  }
  biasTuning = params->biasTuning;

  model_kappa = 4.;

  // Set all the geometry for coordinates.c
  Rh = 1 + sqrt(1. - a * a);  // needed for geodesic steps

  Rin = Rh;
  Rout = rmax_geo;

  Rmax_record = 1.e4;  // this should be large enough that the source looks small

  double z1 = 1. + pow(1. - a * a, 1. / 3.) * (pow(1. + a, 1. / 3.) + pow(1. - a, 1. / 3.));
  double z2 = sqrt(3. * a * a + z1 * z1);
  r_isco = 3. + z2 - copysign(sqrt((3. - z1) * (3. + z1 + 2. * z2)), a);
  startx[0] = 0.0;
  startx[1] = log(Rin);
  startx[2] = 0.0;
  startx[3] = 0.0;

  // set units
  L_unit = GNEWT * MBH_solar * MSUN / (CL * CL);
  RHO_unit = Ne_unit * (MP + ME);
  B_unit = CL * sqrt(4.*M_PI*RHO_unit);

  // the larger this is, the thinner the surface zones -> recover low frequency behavior
  N1 = 512;
  N2 = 512;
  N3 = 1;

  dx[0] = 0.;
  dx[1] = ( log(Rout) - log(Rin) ) / N1;
  dx[2] = 1. / N2;
  dx[3] = 2. * M_PI / N3;

  stopx[0] = 1.;
  stopx[1] = startx[1]+N1*dx[1];
  stopx[2] = startx[2]+N2*dx[2];
  stopx[3] = startx[3]+N3*dx[3];

  // set other units. THESE MAY BE USED ELSEWHERE WITHOUT WARNING
  T_unit = L_unit / CL;
  M_unit = RHO_unit * pow(L_unit, 3);
  max_tau_scatt = (6. * L_unit) * RHO_unit * 0.4;

  geom = (struct of_geom**)malloc_rank2(N1, N2, sizeof(struct of_geom));
  init_geometry();

  tetrads = (struct of_tetrads***)malloc_rank3(N1, N2, N3, sizeof(struct of_tetrads));
  init_tetrads();

  n2gens = (double ***)malloc_rank3(N1, N2, N3, sizeof(double));

  fprintf(stderr, "Running HJW-style RIAF model with a=%g, nth0=%g, Te0=%g, disk_h=%g, pow_nth=%g, pow_T=%g\n",
          a, nth0, Te0, disk_h, pow_nth, pow_T);
  fprintf(stderr, "Velocity model: omega mix %g, radial mix %g\n",
          keplerian_factor, infall_factor);
}

//////////////////////////////////// OUTPUT ////////////////////////////////////

void report_spectrum(int N_superph_made, Params *params)
{
  hid_t fid = -1;

  if (params->loaded && strlen(params->spectrum) > 0) {
    fid = H5Fcreate(params->spectrum, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  } else {
    fid = H5Fcreate("spectrum.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }

  if (fid < 0) {
    fprintf(stderr, "! unable to open/create spectrum hdf5 file.\n");
    exit(-3);
  }

  h5io_add_attribute_str(fid, "/", "githash", xstr(VERSION));

  h5io_add_group(fid, "/params");

  #if CUSTOM_AVG==1
  h5io_add_data_dbl(fid, "/params/CA_MIN_FREQ", CA_MIN_FREQ);
  h5io_add_data_dbl(fid, "/params/CA_MAX_FREQ", CA_MAX_FREQ);
  h5io_add_data_dbl(fid, "/params/CA_MIN_QTY", CA_MIN_QTY);
  h5io_add_data_dbl(fid, "/params/CA_MAX_QTY", CA_MAX_QTY);
  h5io_add_data_dbl(fid, "/params/CA_NBINS", CA_NBINS);
  #endif // CUSTOM_AVG

  h5io_add_data_dbl(fid, "/params/NUCUT", NUCUT);
  h5io_add_data_dbl(fid, "/params/GAMMACUT", GAMMACUT);
  h5io_add_data_dbl(fid, "/params/NUMAX", NUMAX);
  h5io_add_data_dbl(fid, "/params/NUMIN", NUMIN);
  h5io_add_data_dbl(fid, "/params/LNUMAX", LNUMAX);
  h5io_add_data_dbl(fid, "/params/LNUMIN", LNUMIN);
  h5io_add_data_dbl(fid, "/params/DLNU", DLNU);
  h5io_add_data_dbl(fid, "/params/THETAE_MIN", THETAE_MIN);
  h5io_add_data_dbl(fid, "/params/THETAE_MAX", THETAE_MAX);
  h5io_add_data_dbl(fid, "/params/WEIGHT_MIN", WEIGHT_MIN);
  h5io_add_data_dbl(fid, "/params/KAPPA", model_kappa);
  h5io_add_data_dbl(fid, "/params/L_unit", L_unit);
  h5io_add_data_dbl(fid, "/params/T_unit", T_unit);
  h5io_add_data_dbl(fid, "/params/Thetae_unit", Thetae_unit);
  h5io_add_data_dbl(fid, "/params/Rin", Rin);
  h5io_add_data_dbl(fid, "/params/Rout", Rmax);
  h5io_add_data_dbl(fid, "/params/bias", biasTuning);

  h5io_add_data_dbl(fid, "/params/MBH_solar", MBH_solar);
  h5io_add_data_dbl(fid, "/params/a", a);
  h5io_add_data_dbl(fid, "/params/nth0", nth0);
  h5io_add_data_dbl(fid, "/params/Te0", Te0);
  h5io_add_data_dbl(fid, "/params/disk_h", disk_h);
  h5io_add_data_dbl(fid, "/params/pow_nth", pow_nth);
  h5io_add_data_dbl(fid, "/params/pow_T", pow_T);
  h5io_add_data_dbl(fid, "/params/keplerian_factor", keplerian_factor);
  h5io_add_data_dbl(fid, "/params/infall_factor", infall_factor);

  h5io_add_data_dbl(fid, "/params/Te_unit", Te_unit);
  h5io_add_data_dbl(fid, "/params/Ne_unit", Ne_unit);

  h5io_add_data_int(fid, "/params/SYNCHROTRON", SYNCHROTRON);
  h5io_add_data_int(fid, "/params/BREMSSTRAHLUNG", BREMSSTRAHLUNG);
  h5io_add_data_int(fid, "/params/COMPTON", COMPTON);
  h5io_add_data_int(fid, "/params/DIST_KAPPA", MODEL_EDF==EDF_KAPPA_FIXED?1:0);
  h5io_add_data_int(fid, "/params/N_ESAMP", N_ESAMP);
  h5io_add_data_int(fid, "/params/N_EBINS", N_EBINS);
  h5io_add_data_int(fid, "/params/N_THBINS", N_THBINS);
  h5io_add_data_int(fid, "/params/N1", N1);
  h5io_add_data_int(fid, "/params/N2", N2);
  h5io_add_data_int(fid, "/params/N3", N3);
  h5io_add_data_int(fid, "/params/Ns", Ns);

  h5io_add_data_str(fid, "/params/model", xstr(MODEL));

  // temporary data buffers
  double lnu_buf[N_EBINS];
  double dOmega_buf[N_THBINS];
  // USER PATCH: keep large output buffers on the heap so large N_THBINS does not overflow the Windows stack.
  double (*nuLnu_buf)[N_EBINS][N_THBINS] =
      malloc_rank1(N_TYPEBINS * N_EBINS * N_THBINS, sizeof(double));
  double (*tau_abs_buf)[N_EBINS][N_THBINS] =
      malloc_rank1(N_TYPEBINS * N_EBINS * N_THBINS, sizeof(double));
  double (*tau_scatt_buf)[N_EBINS][N_THBINS] =
      malloc_rank1(N_TYPEBINS * N_EBINS * N_THBINS, sizeof(double));
  double (*x1av_buf)[N_EBINS][N_THBINS] =
      malloc_rank1(N_TYPEBINS * N_EBINS * N_THBINS, sizeof(double));
  double (*x2av_buf)[N_EBINS][N_THBINS] =
      malloc_rank1(N_TYPEBINS * N_EBINS * N_THBINS, sizeof(double));
  double (*x3av_buf)[N_EBINS][N_THBINS] =
      malloc_rank1(N_TYPEBINS * N_EBINS * N_THBINS, sizeof(double));
  double (*nscatt_buf)[N_EBINS][N_THBINS] =
      malloc_rank1(N_TYPEBINS * N_EBINS * N_THBINS, sizeof(double));
  double Lcomponent_buf[N_TYPEBINS];
  double (*nph_buf)[N_EBINS][N_THBINS] =
      malloc_rank1(N_TYPEBINS * N_EBINS * N_THBINS, sizeof(double));

  // normal output routine
  double dOmega, nuLnu, tau_scatt, L, Lcomponent, dL;

  max_tau_scatt = 0.;
  L = 0.;
  dL = 0.;

  for (int j=0; j<N_THBINS; ++j) {
    // warning: this assumes geodesic X \in [0,1]
    dOmega_buf[j] = 2. * dOmega_func(j);
  }

  for (int k=0; k<N_TYPEBINS; ++k) {
    Lcomponent = 0.;
    for (int i=0; i<N_EBINS; ++i) {
      lnu_buf[i] = (i * dlE + lE0) / M_LN10;
      for (int j=0; j<N_THBINS; ++j) {

        dOmega = dOmega_buf[j];

        nuLnu = (ME * CL * CL) * (4. * M_PI / dOmega) * (1. / dlE);
        nuLnu *= spect[k][j][i].dEdlE/LSUN;

        tau_scatt = spect[k][j][i].tau_scatt/(spect[k][j][i].dNdlE + SMALL);

        nuLnu_buf[k][i][j] = nuLnu;
        tau_abs_buf[k][i][j] = spect[k][j][i].tau_abs/(spect[k][j][i].dNdlE + SMALL);
        tau_scatt_buf[k][i][j] = tau_scatt;
        x1av_buf[k][i][j] = spect[k][j][i].X1iav/(spect[k][j][i].dNdlE + SMALL);
        x2av_buf[k][i][j] = sqrt(fabs(spect[k][j][i].X2isq/(spect[k][j][i].dNdlE + SMALL)));
        x3av_buf[k][i][j] = sqrt(fabs(spect[k][j][i].X3fsq/(spect[k][j][i].dNdlE + SMALL)));
        nscatt_buf[k][i][j] = spect[k][j][i].nscatt / (spect[k][j][i].dNdlE + SMALL);

        nph_buf[k][i][j] = spect[k][j][i].nph;

        if (tau_scatt > max_tau_scatt) max_tau_scatt = tau_scatt;

        dL += ME * CL * CL * spect[k][j][i].dEdlE;
        L += nuLnu * dOmega * dlE / (4. * M_PI);
        Lcomponent += nuLnu * dOmega * dlE / (4. * M_PI);
      }
    }
    Lcomponent_buf[k] = Lcomponent;
  }

  h5io_add_group(fid, "/output");

  h5io_add_data_dbl_1d(fid, "/output/lnu", N_EBINS, lnu_buf);
  h5io_add_data_dbl_1d(fid, "/output/dOmega", N_THBINS, dOmega_buf);
  h5io_add_data_dbl_3d(fid, "/output/nuLnu", N_TYPEBINS, N_EBINS, N_THBINS, nuLnu_buf);
  h5io_add_data_dbl_3d(fid, "/output/tau_abs", N_TYPEBINS, N_EBINS, N_THBINS, tau_abs_buf);
  h5io_add_data_dbl_3d(fid, "/output/tau_scatt", N_TYPEBINS, N_EBINS, N_THBINS, tau_scatt_buf);
  h5io_add_data_dbl_3d(fid, "/output/x1av", N_TYPEBINS, N_EBINS, N_THBINS, x1av_buf);
  h5io_add_data_dbl_3d(fid, "/output/x2av", N_TYPEBINS, N_EBINS, N_THBINS, x2av_buf);
  h5io_add_data_dbl_3d(fid, "/output/x3av", N_TYPEBINS, N_EBINS, N_THBINS, x3av_buf);
  h5io_add_data_dbl_3d(fid, "/output/nscatt", N_TYPEBINS, N_EBINS, N_THBINS, nscatt_buf);
  h5io_add_data_dbl_3d(fid, "/output/nph", N_TYPEBINS, N_EBINS, N_THBINS, nph_buf);
  h5io_add_data_dbl_1d(fid, "/output/Lcomponent", N_TYPEBINS, Lcomponent_buf);

  h5io_add_data_int(fid, "/output/Nrecorded", N_superph_recorded);
  h5io_add_data_int(fid, "/output/Nmade", N_superph_made);
  h5io_add_data_int(fid, "/output/Nscattered", N_scatt);

  double Lum = L * LSUN;
  h5io_add_data_dbl(fid, "/output/L", Lum);
  h5io_add_attribute_str(fid, "/output/L", "units", "erg/s");

  free(nuLnu_buf);
  free(tau_abs_buf);
  free(tau_scatt_buf);
  free(x1av_buf);
  free(x2av_buf);
  free(x3av_buf);
  free(nscatt_buf);
  free(nph_buf);

  // diagnostic output to screen
  fprintf(stderr, "\n");

  fprintf(stderr, "MBH = %g Msun\n", MBH_solar);
  fprintf(stderr, "max_tau_scatt = %g\n", max_tau_scatt);
  fprintf(stderr, "L = %g erg/s \n", Lum);

  fprintf(stderr, "\n");

  fprintf(stderr, "N_superph_made = %d\n", N_superph_made);
  fprintf(stderr, "N_superph_scatt = %d\n", N_scatt);
  fprintf(stderr, "N_superph_recorded = %d\n", N_superph_recorded);

  H5Fclose(fid);

}


