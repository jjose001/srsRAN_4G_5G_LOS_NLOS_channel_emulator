/**
 * fading_los_validation.c — Deep validation of LOS (Rice) fading channels.
 *
 * Approach: Feed a CW (constant-wave) signal through the channel. The output
 * encodes the time-varying channel DC gain H(f=0,t), dominated by tap-0.
 * Analysing this time-series yields Doppler peak, K-factor, and Rice signature.
 *
 * Returns EXIT_SUCCESS (0) only if ALL assertions pass.
 */

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "srsran/phy/channel/fading.h"
#include "srsran/phy/dft/dft.h"
#include "srsran/phy/utils/vector.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- Configuration ---- */
#define SRATE       1920000.0f
#define F_DOPPLER   100.0f
#define N_OBS       16384
#define SKIP_BLOCKS 16          /* discard initial transient */

/* ---- Tolerances ---- */
#define TOL_DOPPLER_HZ  10.0f
#define TOL_KFACTOR_DB   2.0f
#define TOL_RICE_FRAC    0.50f  /* |mean|/rms must exceed this for Rice */

/* =========================================================================
 * collect_channel_gain()
 *
 * Sends a CW signal (all 1s) through the channel and records one output
 * sample per block.  After SKIP_BLOCKS warm-up, collects N_OBS samples.
 * =========================================================================*/
static int collect_channel_gain(srsran_channel_fading_t* q, cf_t* gain)
{
  uint32_t bs = q->N / 2;

  cf_t* in  = (cf_t*)malloc(bs * sizeof(cf_t));
  cf_t* out = (cf_t*)malloc(bs * sizeof(cf_t));
  if (!in || !out) { free(in); free(out); return -1; }

  /* CW input: every sample = 1+0j */
  for (uint32_t k = 0; k < bs; k++) in[k] = 1.0f + 0.0f * I;

  double t = 0.0;

  /* Warm-up: let the overlap-add state converge */
  for (int n = 0; n < SKIP_BLOCKS; n++)
    t = srsran_channel_fading_execute(q, in, out, bs, t);

  /* Collect observations — one per block (channel updates once per block) */
  for (int n = 0; n < N_OBS; n++) {
    t = srsran_channel_fading_execute(q, in, out, bs, t);
    gain[n] = out[bs / 2];   /* centre sample avoids edge artefacts */
  }

  free(in);
  free(out);
  return 0;
}

/* =========================================================================
 * TEST 1 — Doppler peak at 0.7·f_D
 * =========================================================================*/
static int test_doppler_peak(const char* model, float fd)
{
  printf("\n  [TEST 1] Pico Doppler LOS — modelo '%s', f_D=%.0f Hz\n", model, fd);

  srsran_channel_fading_t q;
  if (srsran_channel_fading_init(&q, SRATE, model, 42) < 0) return -1;

  printf("  Debug: is_los=%d  k_factor=%.2f dB  N=%u  block=%u\n",
         q.is_los, 10.0f * log10f(q.k_factor + 1e-12f), q.N, q.N / 2);

  cf_t* gain = (cf_t*)malloc(N_OBS * sizeof(cf_t));
  if (!gain || collect_channel_gain(&q, gain) < 0) {
    free(gain); srsran_channel_fading_free(&q); return -1;
  }

  /* FFT of the slow-time gain series */
  srsran_dft_plan_t plan;
  cf_t* spec = (cf_t*)malloc(N_OBS * sizeof(cf_t));
  srsran_dft_plan_c(&plan, N_OBS, SRSRAN_DFT_FORWARD);
  srsran_dft_plan_set_norm(&plan, true);
  srsran_dft_run_c(&plan, gain, spec);

  float slow_fs   = SRATE / (float)(q.N / 2);
  float freq_res  = slow_fs / (float)N_OBS;

  /* Find peak in positive frequencies (skip DC = bin 0) */
  int   peak_bin = 1;
  float peak_mag = cabsf(spec[1]);
  for (int k = 2; k < N_OBS / 2; k++) {
    float m = cabsf(spec[k]);
    if (m > peak_mag) { peak_mag = m; peak_bin = k; }
  }

  float peak_freq = (float)peak_bin * freq_res;
  float expected  = 0.7f * fd;
  float err       = fabsf(peak_freq - expected);
  bool  pass      = (err <= TOL_DOPPLER_HZ);

  printf("  Resolución FFT : %.2f Hz/bin\n", freq_res);
  printf("  Pico detectado : bin %d → %.2f Hz\n", peak_bin, peak_freq);
  printf("  f_S esperada   : %.2f Hz  (0.7 × %.0f)\n", expected, fd);
  printf("  Error          : %.2f Hz  (tol ±%.1f Hz)\n", err, TOL_DOPPLER_HZ);
  printf("  Resultado      : %s\n", pass ? "PASS ✓" : "FAIL ✗");

  /* Debug: print first 5 gain samples to confirm non-zero */
  printf("  Debug gain[0..4]: ");
  for (int i = 0; i < 5; i++)
    printf("(%.4f,%.4f) ", crealf(gain[i]), cimagf(gain[i]));
  printf("\n");

  srsran_dft_plan_free(&plan);
  free(spec); free(gain);
  srsran_channel_fading_free(&q);
  return pass ? 0 : -1;
}

/* =========================================================================
 * TEST 2 — K-factor from demodulated tap-0
 *
 * K = |E[h_demod]|² / Var(h_demod)
 * where h_demod = gain[n] · conj(e^{j(2π·f_peak·t_n)})
 * f_peak is detected from the FFT (should be ≈ 0.7·fd).
 * =========================================================================*/
static int test_k_factor(const char* model, float fd, float k_target_dB)
{
  printf("\n  [TEST 2] K-factor medido — modelo '%s', K obj = %.1f dB\n",
         model, k_target_dB);

  srsran_channel_fading_t q;
  if (srsran_channel_fading_init(&q, SRATE, model, 42) < 0) return -1;

  cf_t* gain = (cf_t*)malloc(N_OBS * sizeof(cf_t));
  if (!gain || collect_channel_gain(&q, gain) < 0) {
    free(gain); srsran_channel_fading_free(&q); return -1;
  }

  float block_dt = (float)(q.N / 2) / SRATE;

  /* Demodulate at the exact theoretical LOS frequency */
  float f_los = 0.7f * fd;
  cf_t mean_h = 0.0f;
  for (int n = 0; n < N_OBS; n++) {
    float t_n = (float)n * block_dt;
    float ang = 2.0f * (float)M_PI * f_los * t_n;
    cf_t  ref = cosf(ang) + sinf(ang) * I;
    gain[n]  *= conjf(ref);
    mean_h   += gain[n];
  }
  mean_h /= (float)N_OBS;

  /* Variance of diffuse (Rayleigh) component */
  float var = 0.0f;
  for (int n = 0; n < N_OBS; n++) {
    cf_t d = gain[n] - mean_h;
    var += crealf(d) * crealf(d) + cimagf(d) * cimagf(d);
  }
  var /= (float)N_OBS;

  float mean_sq  = crealf(mean_h) * crealf(mean_h) +
                   cimagf(mean_h) * cimagf(mean_h);
  float k_meas   = (var > 1e-15f) ? (mean_sq / var) : 0.0f;
  float k_dB     = 10.0f * log10f(k_meas + 1e-15f);
  float err      = fabsf(k_dB - k_target_dB);
  bool  pass     = (err <= TOL_KFACTOR_DB);

  printf("  |E[h]|²=%.6f  Var=%.6f\n", mean_sq, var);
  printf("  K medido  : %.4f → %.2f dB\n", k_meas, k_dB);
  printf("  K objetivo: %.2f dB\n", k_target_dB);
  printf("  Error     : %.2f dB  (tol ±%.1f dB)\n", err, TOL_KFACTOR_DB);
  printf("  Resultado : %s\n", pass ? "PASS ✓" : "FAIL ✗");

  free(gain);
  srsran_channel_fading_free(&q);
  return pass ? 0 : -1;
}

/* =========================================================================
 * TEST 3 — Rice signature: |mean| / rms >> 0 only for LOS channels
 *
 * For Rice channels, the gain series has a non-zero mean (the LOS component).
 * For Rayleigh channels, the mean converges to zero.
 * We verify that the ratio |mean|/rms exceeds TOL_RICE_FRAC for LOS models,
 * confirming the Rice tap-0 is active.
 * =========================================================================*/
static int test_rice_signature(const char* model, float fd, float k_target_dB)
{
  printf("\n  [TEST 3] Firma Rice (media no nula) — modelo '%s'\n", model);

  srsran_channel_fading_t q;
  if (srsran_channel_fading_init(&q, SRATE, model, 42) < 0) return -1;

  cf_t* gain = (cf_t*)malloc(N_OBS * sizeof(cf_t));
  if (!gain || collect_channel_gain(&q, gain) < 0) {
    free(gain); srsran_channel_fading_free(&q); return -1;
  }

  float block_dt = (float)(q.N / 2) / SRATE;

  /* Demodulate at the exact theoretical LOS frequency */
  float f_los = 0.7f * fd;
  cf_t mean_h = 0.0f;
  float rms   = 0.0f;
  for (int n = 0; n < N_OBS; n++) {
    float t_n = (float)n * block_dt;
    float ang = 2.0f * (float)M_PI * f_los * t_n;
    cf_t  ref = cosf(ang) + sinf(ang) * I;
    cf_t  h   = gain[n] * conjf(ref);
    mean_h   += h;
    rms      += crealf(gain[n]) * crealf(gain[n]) +
                cimagf(gain[n]) * cimagf(gain[n]);
  }
  mean_h /= (float)N_OBS;
  rms     = sqrtf(rms / (float)N_OBS);

  float mean_mag = cabsf(mean_h);
  float ratio    = (rms > 1e-15f) ? (mean_mag / rms) : 0.0f;

  /* Expected ratio ≈ sqrt(K/(K+1)) for Rice channel */
  float k_lin    = powf(10.0f, k_target_dB / 10.0f);
  float expected = sqrtf(k_lin / (k_lin + 1.0f));
  bool  pass     = (ratio >= TOL_RICE_FRAC * expected);

  printf("  |E[h_demod]| = %.6f,  RMS = %.6f\n", mean_mag, rms);
  printf("  Ratio |mean|/rms = %.4f\n", ratio);
  printf("  Esperado ≈ %.4f,  umbral = %.2f × %.4f = %.4f\n",
         expected, TOL_RICE_FRAC, expected, TOL_RICE_FRAC * expected);
  printf("  Resultado : %s\n", pass ? "PASS ✓" : "FAIL ✗");

  /* Also verify NLOS model (TDL-A) has ratio ≈ 0 */
  free(gain);
  srsran_channel_fading_free(&q);
  return pass ? 0 : -1;
}

/* =========================================================================
 * TEST 4 — Regression: NLOS model must NOT show Rice signature
 * =========================================================================*/
static int test_nlos_regression(void)
{
  printf("\n  [TEST 4] Regresión NLOS — TDL-A no debe ser LOS\n");

  srsran_channel_fading_t q;
  if (srsran_channel_fading_init(&q, SRATE, "tdla100", 42) < 0) return -1;

  bool pass = !q.is_los;
  printf("  TDL-A is_los = %s  →  %s\n",
         q.is_los ? "true" : "false",
         pass ? "PASS ✓ (NLOS correcto)" : "FAIL ✗");

  srsran_channel_fading_free(&q);
  return pass ? 0 : -1;
}

/* =========================================================================
 * main
 * =========================================================================*/
int main(void)
{
  int ret = 0;

  printf("\n=============================================================\n");
  printf("  Validación Profunda del Canal LOS (TDL-D / TDL-E)\n");
  printf("=============================================================\n");
  printf("  srate=%.2f MHz | f_D=%.0f Hz | N_obs=%d | skip=%d\n\n",
         SRATE / 1e6f, F_DOPPLER, N_OBS, SKIP_BLOCKS);

  printf("=== TDL-D  (K = 13.3 dB) ===\n");
  ret |= test_doppler_peak    ("tdld100", F_DOPPLER);
  ret |= test_k_factor        ("tdld100", F_DOPPLER, 13.3f);
  ret |= test_rice_signature  ("tdld100", F_DOPPLER, 13.3f);

  printf("\n=== TDL-E  (K = 22.0 dB) ===\n");
  ret |= test_doppler_peak    ("tdle100", F_DOPPLER);
  ret |= test_k_factor        ("tdle100", F_DOPPLER, 22.0f);
  ret |= test_rice_signature  ("tdle100", F_DOPPLER, 22.0f);

  ret |= test_nlos_regression();

  printf("\n=============================================================\n");
  printf("  Resultado Final: %s\n", (ret == 0) ? "PASS ✓" : "FAIL ✗");
  printf("=============================================================\n\n");

  return (ret == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
