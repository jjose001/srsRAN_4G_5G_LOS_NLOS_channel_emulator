/**
 * fading_nlos_validation.c — Validación de Canales Rayleigh (NLOS)
 *
 * Test de regresión para asegurar que los modelos TDL-A/B/C (NLOS)
 * mantienen su comportamiento puramente difuso (Rayleigh/Jakes) y
 * no se ven afectados por la implementación LOS añadida a TDL-D/E.
 *
 * Test 1: Espectro Doppler confinado a [-f_D, +f_D] (Jakes, sin pico LOS)
 * Test 2: Firma puramente Rayleigh (media ≈ 0, sin componente determinista)
 * Test 3: Banderas internas correctas (is_los == false, k_factor == 0)
 *
 * Metodología: Inyección CW (x=1+j0) → la salida codifica la ganancia
 * DC del canal H(f=0,t), cuyo espectro y estadísticas se analizan.
 */

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "srsran/phy/channel/fading.h"
#include "srsran/phy/dft/dft.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- Configuración ---- */
#define SRATE       1920000.0f
#define F_DOPPLER   100.0f
#define N_OBS       16384
#define SKIP_BLOCKS 16

/* ---- Tolerancias ---- */
#define TOL_RAYLEIGH_RATIO 0.15f  /* |mean|/rms debe ser < esto para Rayleigh */
#define TOL_OUT_OF_BAND    0.05f  /* máxima energía fuera de ±1.2·f_D (5%) */

/* =========================================================================
 * collect_channel_gain()
 *
 * Envía señal CW (1+j0) por el canal y recoge un sample representativo
 * de la salida por bloque.  Idéntico al método del test LOS.
 * =========================================================================*/
static int collect_channel_gain(srsran_channel_fading_t* q, cf_t* gain)
{
  uint32_t bs = q->N / 2;
  cf_t* in  = (cf_t*)malloc(bs * sizeof(cf_t));
  cf_t* out = (cf_t*)malloc(bs * sizeof(cf_t));
  if (!in || !out) { free(in); free(out); return -1; }

  for (uint32_t k = 0; k < bs; k++) in[k] = 1.0f + 0.0f * I;

  double t = 0.0;
  for (int n = 0; n < SKIP_BLOCKS; n++)
    t = srsran_channel_fading_execute(q, in, out, bs, t);

  for (int n = 0; n < N_OBS; n++) {
    t = srsran_channel_fading_execute(q, in, out, bs, t);
    gain[n] = out[bs / 2];
  }

  free(in); free(out);
  return 0;
}

/* =========================================================================
 * TEST 1 — Espectro Doppler confinado (forma U de Jakes, sin pico LOS)
 *
 * Verifica que toda la energía espectral del canal cae dentro de ±1.2·f_D.
 * Si existiera un pico LOS indebido o un error en el filtro de Jakes,
 * se detectaría como energía fuera de banda excesiva.
 * =========================================================================*/
static int test_doppler_spread(const char* model, float fd)
{
  printf("\n  [TEST 1] Densidad Espectral Doppler — modelo '%s'\n", model);

  srsran_channel_fading_t q;
  if (srsran_channel_fading_init(&q, SRATE, model, 42) < 0) return -1;

  printf("  Debug: N=%u  block=%u\n", q.N, q.N / 2);

  cf_t* gain = (cf_t*)malloc(N_OBS * sizeof(cf_t));
  if (collect_channel_gain(&q, gain) < 0) return -1;

  /* FFT de la serie temporal lenta de ganancia */
  srsran_dft_plan_t plan;
  cf_t* spec = (cf_t*)malloc(N_OBS * sizeof(cf_t));
  srsran_dft_plan_c(&plan, N_OBS, SRSRAN_DFT_FORWARD);
  srsran_dft_plan_set_norm(&plan, true);
  srsran_dft_run_c(&plan, gain, spec);

  float slow_fs  = SRATE / (float)(q.N / 2);
  float freq_res = slow_fs / (float)N_OBS;

  float total_power = 0.0f;
  float oob_power   = 0.0f;

  for (int k = 0; k < N_OBS; k++) {
    float f = (k < N_OBS / 2) ? (k * freq_res) : ((k - N_OBS) * freq_res);
    float pwr = crealf(spec[k])*crealf(spec[k]) + cimagf(spec[k])*cimagf(spec[k]);
    total_power += pwr;
    if (fabsf(f) > 1.2f * fd) oob_power += pwr;
  }

  float oob_ratio = total_power > 0 ? (oob_power / total_power) : 1.0f;
  bool pass = (oob_ratio <= TOL_OUT_OF_BAND);

  printf("  Energía confinada en ±%.0f Hz : %.2f%%\n", 1.2f * fd, (1.0f - oob_ratio) * 100.0f);
  printf("  Energía fuera de banda      : %.2f%% (tolerancia < %.0f%%)\n",
         oob_ratio * 100.0f, TOL_OUT_OF_BAND * 100.0f);
  printf("  Resultado                   : %s\n", pass ? "PASS ✓" : "FAIL ✗");

  srsran_dft_plan_free(&plan); free(spec); free(gain); srsran_channel_fading_free(&q);
  return pass ? 0 : -1;
}

/* =========================================================================
 * TEST 2 — Firma Puramente Rayleigh (Media ≈ 0)
 *
 * Un canal Rayleigh puro tiene media compleja nula (todas las fases son
 * aleatorias y se cancelan). Verificamos:
 *   a) |E[h]| / RMS ≈ 0 (media bruta)
 *   b) |E[h · e^{-j2π·0.7fD·t}]| / RMS ≈ 0 (demodulado en la frecuencia
 *      donde estaría el pico LOS si se hubiera inyectado erróneamente)
 * =========================================================================*/
static int test_rayleigh_signature(const char* model, float fd)
{
  printf("\n  [TEST 2] Firma Rayleigh (media ≈ 0) — modelo '%s'\n", model);

  srsran_channel_fading_t q;
  if (srsran_channel_fading_init(&q, SRATE, model, 42) < 0) return -1;

  cf_t* gain = (cf_t*)malloc(N_OBS * sizeof(cf_t));
  if (collect_channel_gain(&q, gain) < 0) return -1;

  /* a) Media bruta */
  cf_t raw_mean = 0.0f;
  float rms = 0.0f;
  for (int n = 0; n < N_OBS; n++) {
    raw_mean += gain[n];
    rms += crealf(gain[n])*crealf(gain[n]) + cimagf(gain[n])*cimagf(gain[n]);
  }
  raw_mean /= (float)N_OBS;
  rms = sqrtf(rms / (float)N_OBS);

  float raw_ratio = (rms > 1e-15f) ? (cabsf(raw_mean) / rms) : 0.0f;

  /* b) Demodular en 0.7·f_D para detectar componente LOS accidental */
  float block_dt = (float)(q.N / 2) / SRATE;
  float f_los = 0.7f * fd;
  cf_t dom_mean = 0.0f;
  for (int n = 0; n < N_OBS; n++) {
    float ang = 2.0f * (float)M_PI * f_los * (float)n * block_dt;
    cf_t ref = cosf(ang) + sinf(ang) * I;
    dom_mean += gain[n] * conjf(ref);
  }
  dom_mean /= (float)N_OBS;
  float dom_ratio = (rms > 1e-15f) ? (cabsf(dom_mean) / rms) : 0.0f;

  bool pass = (raw_ratio < TOL_RAYLEIGH_RATIO) && (dom_ratio < TOL_RAYLEIGH_RATIO);

  printf("  |E[h]| / RMS               = %.4f (< %.2f)\n", raw_ratio, TOL_RAYLEIGH_RATIO);
  printf("  |E[h·e^{-j0.7fD}]| / RMS   = %.4f (< %.2f)\n", dom_ratio, TOL_RAYLEIGH_RATIO);
  printf("  Resultado                   : %s\n", pass ? "PASS ✓" : "FAIL ✗");

  free(gain); srsran_channel_fading_free(&q);
  return pass ? 0 : -1;
}

/* =========================================================================
 * TEST 3 — Banderas internas (is_los == false, k_factor == 0)
 * =========================================================================*/
static int test_is_los_flag(const char* model)
{
  printf("\n  [TEST 3] Banderas internas — modelo '%s'\n", model);
  srsran_channel_fading_t q;
  if (srsran_channel_fading_init(&q, SRATE, model, 42) < 0) return -1;

  bool pass = !q.is_los && (q.k_factor == 0.0f);

  printf("  q->is_los   = %s  →  %s\n",
         q.is_los ? "true (ERROR)" : "false (OK)",
         !q.is_los ? "PASS ✓" : "FAIL ✗");
  printf("  q->k_factor = %.4f  →  %s\n",
         q.k_factor,
         (q.k_factor == 0.0f) ? "PASS ✓" : "FAIL ✗");

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
  printf("  Validación Canales NLOS Rayleigh (TDL-A / TDL-B / TDL-C)\n");
  printf("=============================================================\n");
  printf("  srate=%.2f MHz | f_D=%.0f Hz | N_obs=%d | skip=%d\n",
         SRATE / 1e6f, F_DOPPLER, N_OBS, SKIP_BLOCKS);

  const char* models[] = {"tdla100", "tdlb100", "tdlc100"};

  for (int i = 0; i < 3; i++) {
    printf("\n=== MODELO: %s ===\n", models[i]);
    ret |= test_doppler_spread(models[i], F_DOPPLER);
    ret |= test_rayleigh_signature(models[i], F_DOPPLER);
    ret |= test_is_los_flag(models[i]);
  }

  printf("\n=============================================================\n");
  printf("  Resultado Final NLOS: %s\n", (ret == 0) ? "PASS ✓" : "FAIL ✗");
  printf("=============================================================\n\n");

  return (ret == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
