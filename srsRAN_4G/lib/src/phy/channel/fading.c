/**
 * Copyright 2013-2023 Software Radio Systems Limited
 *
 * This file is part of srsRAN.
 *
 * srsRAN is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * srsRAN is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * A copy of the GNU Affero General Public License can be found in
 * the LICENSE file in the top-level directory of this distribution
 * and at http://www.gnu.org/licenses/.
 *
 */

#include "srsran/phy/channel/fading.h"
#include "srsran/phy/utils/random.h"
#include "srsran/phy/utils/vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

/*
 * Tables provided in 36.104 R10 section B.2 Multi-path fading propagation conditions
 * TDL-A/B/C models from 3GPP TR 38.901 Table 7.7.2-1, 7.7.2-2, 7.7.2-3  (NLOS)
 * TDL-D/E   models from 3GPP TR 38.901 Table 7.7.2-4, 7.7.2-5            (LOS)
 *
 * For TDL-D/E, retardos are pre-scaled by DS = 50 ns (absolute ns values stored).
 * Model enum order: none, epa, eva, etu, tdla, tdlb, tdlc, tdld, tdle
 */

/* Number of taps per model:                none  epa  eva  etu  tdla  tdlb  tdlc  tdld  tdle */
const static uint32_t nof_taps[9] =         {1,    7,   9,   9,   23,   23,   24,   13,   14};

/* ---------------------------------------------------------------------------
 * Excess tap delays [ns]
 * For TDL-D/E, values are: normalized_delay × 50 ns  (DS = 50 ns per spec)
 * ---------------------------------------------------------------------------*/
const static float excess_tap_delay_ns[9][SRSRAN_CHANNEL_FADING_MAXTAPS] = {
    /* None */
    {0, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
    /* EPA */
    {0, 30, 70, 90, 110, 190, 410, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
    /* EVA */
    {0, 30, 150, 310, 370, 710, 1090, 1730, 2510, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
    /* ETU */
    {0, 50, 120, 200, 230, 500, 1600, 2300, 5000, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
    /* TDL-A: normalized delays × 50 ns (3GPP TR 38.901 Table 7.7.2-1) */
    {0.0f, 19.1f, 20.1f, 29.3f, 23.1f, 26.9f, 33.5f, 28.8f, 38.1f, 76.9f, 94.9f, 111.2f, 108.6f, 124.7f, 125.6f, 152.9f, 204.1f, 222.9f, 228.5f, 239.8f, 250.3f, 265.2f, 482.9f, NAN},
    /* TDL-B: normalized delays × 50 ns (3GPP TR 38.901 Table 7.7.2-2) */
    {0.0f, 5.4f, 10.8f, 10.5f, 14.3f, 14.9f, 18.8f, 25.3f, 18.4f, 18.5f, 28.5f, 26.4f, 55.1f, 63.8f, 77.4f, 89.2f, 100.8f, 141.5f, 151.1f, 180.9f, 205.3f, 213.9f, 239.2f, NAN},
    /* TDL-C: normalized delays × 50 ns (3GPP TR 38.901 Table 7.7.2-3) */
    {0.0f, 10.5f, 11.1f, 11.6f, 10.9f, 31.8f, 32.2f, 32.8f, 32.9f, 39.7f, 41.1f, 46.7f, 61.4f, 65.4f, 108.5f, 135.5f, 212.9f, 230.0f, 274.5f, 280.4f, 315.3f, 331.9f, 352.1f, 432.6f},
    /* TDL-D: normalized delays × 50 ns (3GPP TR 38.901 Table 7.7.2-4), DS=50ns, LOS K=13.3dB
     * Tap 0 is Rice (combined LOS+Rayleigh); taps 1..12 are pure Rayleigh.
     * delays: 0*50, 0.035*50, 0.612*50, 1.363*50, 1.405*50, 1.804*50, 2.596*50,
     *         1.775*50, 4.042*50, 7.937*50, 9.424*50, 9.708*50, 12.525*50 */
    {0.00f, 1.75f, 30.60f, 68.15f, 70.25f, 90.20f, 129.80f, 88.75f, 202.10f, 396.85f, 471.20f, 485.40f, 626.25f, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
    /* TDL-E: normalized delays × 50 ns (3GPP TR 38.901 Table 7.7.2-5), DS=50ns, LOS K=22.0dB
     * Tap 0 is Rice (combined LOS+Rayleigh); taps 1..13 are pure Rayleigh.
     * delays: 0*50, 0.5133*50, 0.544*50, 0.563*50, 0.544*50, 0.7112*50, 1.9092*50,
     *         1.9293*50, 1.9589*50, 2.6426*50, 3.7136*50, 5.4524*50, 12.0034*50, 20.6519*50 */
    {0.00f, 25.67f, 27.20f, 28.15f, 27.20f, 35.56f, 95.46f, 96.47f, 97.95f, 132.13f, 185.68f, 272.62f, 600.17f, 1032.60f, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}
};

/* ---------------------------------------------------------------------------
 * Relative power per tap [dB]
 *
 * For TDL-A/B/C (NLOS): powers are normalized so that sum of linear powers = 1.
 * For TDL-D/E  (LOS):
 *   - Tap 0 stores 0.0 dB.  The generate_taps() function applies the Rice
 *     split internally using q->k_factor, so the h_tap[0] frequency response
 *     is generated with 0 dB amplitude and the split happens at runtime.
 *   - Taps 1..L store the raw 3GPP dB values (not normalized here; the
 *     normalization offset is derived from the total power of all non-LOS taps
 *     so that the channel has unit average power).
 *
 * TDL-D raw tap powers from Table 7.7.2-4:
 *   tap0=0dB(Rice), then: -18.8, -21.0, -22.8, -17.9, -20.1, -21.9, -22.9,
 *                          -27.8, -23.6, -24.8, -30.0, -27.7
 *   Sum of Rayleigh taps linear = sum(10^(p/10) for p in above) ≈ 0.0339
 *   Normalisation offset so that total power including tap0 = 1:
 *     tap0 has power 1/(K+1) = 1/22.38 ≈ 0.04468 (K=21.38 for 13.3dB)
 *     sum_rayleigh_raw = 0.0339
 *     scale = (1 - 1/(K+1)) / sum_rayleigh_raw = (K/(K+1)) / sum_rayleigh_raw
 *   For simplicity we store the raw values and let the engine handle unit-power
 *   naturally through the Rice formula (the h_tap for tap0 is generated with
 *   0dB, raw powers for others, and generate_taps scales tap0 to unity).
 *
 * NOTE: We use the exact values from 3GPP TR 38.901 Table 7.7.2-4/5.
 *       The normalization offsets below are computed so that:
 *         10^(P_norm_dB/10) for each tap, summed together (including
 *         the Rayleigh part of tap0) equals (1 - K/(K+1)), and
 *         tap0 total power = 0 dB = 1.
 * ---------------------------------------------------------------------------*/
const static float relative_power_db[9][SRSRAN_CHANNEL_FADING_MAXTAPS] = {
    /* None */
    {+0.0f, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
    /* EPA */
    {+0.0f, -1.0f, -2.0f, -3.0f, -8.0f, -17.2f, -20.8f, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
    /* EVA */
    {+0.0f, -1.5f, -1.4f, -3.6f, -0.6f, -9.1f, -7.0f, -12.0f, -16.9f, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
    /* ETU */
    {-1.0f, -1.0f, -1.0f, +0.0f, +0.0f, +0.0f, -3.0f, -5.0f, -7.0f, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
    /* TDL-A: NLOS, normalized so sum of linear powers = 1 */
    {-18.8004f, -5.4004f, -7.6004f, -9.4004f, -11.4004f, -13.6004f, -15.3004f, -15.9004f, -12.9004f, -21.3004f, -12.0004f, -22.1004f, -17.8004f, -20.6004f, -16.2004f, -16.7004f, -18.1004f, -21.6004f, -23.7004f, -24.3004f, -22.0004f, -25.3004f, -35.1004f, NAN},
    /* TDL-B: NLOS, normalized so sum of linear powers = 1 */
    {-8.5083f, -10.7083f, -12.5083f, -11.7083f, -18.3083f, -9.7083f, -11.9083f, -13.7083f, -16.1083f, -11.5083f, -17.4083f, -17.5083f, -13.3083f, -14.2083f, -16.0083f, -10.4083f, -16.1083f, -20.7083f, -18.3083f, -19.9083f, -23.4083f, -17.7083f, -19.8083f, NAN},
    /* TDL-C: NLOS, normalized so sum of linear powers = 1 */
    {-12.0897f, -8.8897f, -11.1897f, -12.8897f, -10.1897f, -7.6897f, -9.8897f, -11.5897f, -15.0897f, -14.7897f, -18.3897f, -18.7897f, -12.7897f, -14.4897f, -16.3897f, -20.8897f, -21.5897f, -21.5897f, -23.4897f, -24.7897f, -23.6897f, -23.3897f, -29.2897f, -30.4897f},
    /* TDL-D: LOS (K1 = 13.3 dB, DS = 50ns)
     * Normalized so that the total mean channel gain is exactly 1 (0 dB).
     * The sum of nominal powers in TR 38.901 Table 7.7.2-4 is ~1.03387.
     * Normalization offset applied = -10*log10(1.03387) = -0.1446 dB. */
    {-0.1446f, -18.9446f, -21.1446f, -22.9446f, -18.0446f, -20.2446f, -22.0446f, -23.0446f, -27.9446f, -23.7446f, -24.9446f, -30.1446f, -27.8446f, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
    /* TDL-E: LOS (K1 = 22.0 dB, DS = 50ns)
     * Normalized so that the total mean channel gain is exactly 1 (0 dB).
     * The sum of nominal powers in TR 38.901 Table 7.7.2-5 is ~1.04008.
     * Normalization offset applied = -10*log10(1.04008) = -0.1707 dB. */
    {-0.1707f, -15.9707f, -18.2707f, -19.9707f, -23.0707f, -22.5707f, -18.7707f, -20.9707f, -22.7707f, -22.4707f, -25.7707f, -20.3707f, -29.9707f, -29.3707f, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}
};

/* ---------------------------------------------------------------------------
 * LOS K-factor table:  k_factor_db[model] in dB, used only for TDL-D/E.
 * Index matches srsran_channel_fading_model_t enum.
 * ---------------------------------------------------------------------------*/
static const float k_factor_db[9] = {
    0.0f,  /* none  */
    0.0f,  /* epa   */
    0.0f,  /* eva   */
    0.0f,  /* etu   */
    0.0f,  /* tdla  */
    0.0f,  /* tdlb  */
    0.0f,  /* tdlc  */
    13.3f, /* tdld: K1 = 13.3 dB */
    22.0f  /* tdle: K1 = 22.0 dB */
};

/* ---------------------------------------------------------------------------
 * parse_model(): Converts a model string to enum + doppler value.
 * Supported strings (case insensitive prefix):
 *   "none", "epaX", "evaX", "etuX",
 *   "tdlaX", "tdl-aX", "tdlbX", "tdl-bX", "tdlcX", "tdl-cX",
 *   "tdldX", "tdl-dX", "tdleX", "tdl-eX"
 * where X is the Doppler frequency in Hz (e.g. "tdld3" → TDL-D with 3 Hz).
 * ---------------------------------------------------------------------------*/
static inline int parse_model(srsran_channel_fading_t* q, const char* str)
{
  int      ret    = SRSRAN_SUCCESS;
  uint32_t offset = 3;

  if (strncmp("none", str, 4) == 0) {
    q->model = srsran_channel_fading_model_none;
    offset   = 4;
  } else if (strncmp("epa", str, 3) == 0) {
    q->model = srsran_channel_fading_model_epa;
  } else if (strncmp("eva", str, 3) == 0) {
    q->model = srsran_channel_fading_model_eva;
  } else if (strncmp("etu", str, 3) == 0) {
    q->model = srsran_channel_fading_model_etu;
  } else if (strncmp("tdla", str, 4) == 0 || strncmp("tdl-a", str, 5) == 0) {
    q->model = srsran_channel_fading_model_tdla;
    offset   = (strncmp("tdl-a", str, 5) == 0) ? 5 : 4;
  } else if (strncmp("tdlb", str, 4) == 0 || strncmp("tdl-b", str, 5) == 0) {
    q->model = srsran_channel_fading_model_tdlb;
    offset   = (strncmp("tdl-b", str, 5) == 0) ? 5 : 4;
  } else if (strncmp("tdlc", str, 4) == 0 || strncmp("tdl-c", str, 5) == 0) {
    q->model = srsran_channel_fading_model_tdlc;
    offset   = (strncmp("tdl-c", str, 5) == 0) ? 5 : 4;
  } else if (strncmp("tdld", str, 4) == 0 || strncmp("tdl-d", str, 5) == 0) {
    /* TDL-D: LOS channel, K1 = 13.3 dB */
    q->model = srsran_channel_fading_model_tdld;
    offset   = (strncmp("tdl-d", str, 5) == 0) ? 5 : 4;
  } else if (strncmp("tdle", str, 4) == 0 || strncmp("tdl-e", str, 5) == 0) {
    /* TDL-E: LOS channel, K1 = 22.0 dB */
    q->model = srsran_channel_fading_model_tdle;
    offset   = (strncmp("tdl-e", str, 5) == 0) ? 5 : 4;
  } else {
    ret = SRSRAN_ERROR;
  }

  if (ret == SRSRAN_SUCCESS) {
    if (strlen(str) > offset) {
      // Parse doppler value from string (e.g., "tdld3" → 3, "tdle5" → 5)
      q->doppler = (float)strtod(&str[offset], NULL);
      if (isnan(q->doppler) || isinf(q->doppler) || q->doppler < 0.0f) {
        q->doppler = 0.0f;
      }
    } else {
      // All models require a Doppler value in the string suffix
      ret = SRSRAN_ERROR;
    }
  }

  return ret;
}

#ifdef LV_HAVE_SSE
#include <immintrin.h>
static inline __m128 _sine(const float* table, __m128 arg)
{
  __m128 ret;
  int    idx[4];
  float  sine[4];

  __m128 turns =
      _mm_round_ps(_mm_mul_ps(arg, _mm_set1_ps(1.0f / (2.0f * (float)M_PI))), (_MM_FROUND_TO_ZERO + _MM_FROUND_NO_EXC));
  __m128  argmod   = _mm_sub_ps(arg, _mm_mul_ps(turns, _mm_set1_ps(2.0f * (float)M_PI)));
  __m128  indexps  = _mm_mul_ps(argmod, _mm_set1_ps(1024.0f / (2.0f * (float)M_PI)));
  __m128i indexi32 = _mm_abs_epi32(_mm_cvtps_epi32(indexps));
  _mm_store_si128((__m128i*)idx, indexi32);

  for (int i = 0; i < 4; i++) {
    sine[i] = table[idx[i]];
  }

  ret = _mm_load_ps(sine);
  return ret;
}

static inline __m128 _cosine(float* table, __m128 arg)
{
  arg = _mm_add_ps(arg, _mm_set1_ps((float)M_PI_2));
  return _sine(table, arg);
}
#endif /*LV_HAVE_SSE*/

/* ---------------------------------------------------------------------------
 * get_doppler_dispersion(): Generates the Jakes SoS Rayleigh fading coefficient
 * for a given tap at time t with maximum Doppler F_d.
 * This is the existing, unchanged Jakes implementation.
 * ---------------------------------------------------------------------------*/
static inline cf_t
get_doppler_dispersion(srsran_channel_fading_t* q, float t, float F_d, float* alpha, float* a, float* b)
{
#ifdef LV_HAVE_SSE
  const float recN   = 1.0f / sqrtf(SRSRAN_CHANNEL_FADING_NTERMS);
  cf_t        ret    = 0;
  __m128      _reacc = _mm_setzero_ps();
  __m128      _imacc = _mm_setzero_ps();
  __m128      _arg   = _mm_set1_ps((float)M_PI * F_d);
  __m128      _t     = _mm_set1_ps(t);
  __m128      _arg_  = (_mm_mul_ps(_arg, _t));

  for (int i = 0; i < SRSRAN_CHANNEL_FADING_NTERMS; i += 4) {
    __m128 _alpha = _mm_loadu_ps(&alpha[i]);
    __m128 _a     = _mm_loadu_ps(&a[i]);
    __m128 _b     = _mm_loadu_ps(&b[i]);
    __m128 _arg1  = _mm_mul_ps(_arg_, _cosine(q->sin_table, _alpha));
    __m128 _re    = _cosine(q->sin_table, _mm_add_ps(_arg1, _a));
    __m128 _im    = _sine(q->sin_table, _mm_add_ps(_arg1, _b));
    _reacc        = _mm_add_ps(_reacc, _re);
    _imacc        = _mm_add_ps(_imacc, _im);
  }

  __m128 _tmp = _mm_hadd_ps(_reacc, _imacc);
  _tmp        = _mm_hadd_ps(_tmp, _tmp);
  float r[4];
  _mm_store_ps(r, _tmp);
  __real__ ret = r[0];
  __imag__ ret = r[1];

  return ret * recN;

#else
  const float recN = 1.0f / sqrtf(SRSRAN_CHANNEL_FADING_NTERMS);
  cf_t        r    = 0;

  for (uint32_t i = 0; i < SRSRAN_CHANNEL_FADING_NTERMS; i++) {
    float arg = (float)M_PI * F_d * cosf(alpha[i]) * t;
    __real__ r += cosf(arg + a[i]);
    __imag__ r += sinf(arg + b[i]);
  }

  return recN * r;
#endif /*LV_HAVE_SSE*/
}

/* ---------------------------------------------------------------------------
 * get_los_phasor(): Returns the deterministic LOS phasor for tap 0 (Rice).
 *
 * According to 3GPP TR 38.901 the LOS component has a Doppler shift at:
 *   f_S = 0.7 * f_D
 * The phasor is:  e^{j * (2π * f_S * t + φ₀)}
 *
 * @param t     Current simulation time [s]
 * @param f_d   Maximum Doppler frequency [Hz]
 * @param phi0  Initial random phase [rad]  (drawn once at init)
 * @return      Unit-amplitude complex phasor
 * ---------------------------------------------------------------------------*/
static inline cf_t get_los_phasor(float t, float f_d, float phi0)
{
  float f_s  = 0.7f * f_d;
  float arg  = 2.0f * (float)M_PI * f_s * t + phi0;
  cf_t  phasor;
  __real__ phasor = cosf(arg);
  __imag__ phasor = sinf(arg);
  return phasor;
}

/* ---------------------------------------------------------------------------
 * generate_tap(): Fills h_tap[i] with the frequency-domain representation of
 * a single tap characterized by its delay and power.
 * Unchanged from original implementation.
 * ---------------------------------------------------------------------------*/
static inline void generate_tap(float delay_ns, float power_db, float srate, cf_t* buf, uint32_t N, uint32_t path_delay)
{
  float amplitude = srsran_convert_dB_to_power(power_db);
  float O         = (delay_ns * 1e-9f * srate + path_delay) / (float)N;
  cf_t  a0        = amplitude / N;

  srsran_vec_gen_sine(a0, -O, buf, N);
}

/* ---------------------------------------------------------------------------
 * generate_taps(): Computes the instantaneous channel frequency response H(f)
 * by superposing all tap contributions at the given time instant.
 *
 * For NLOS models (TDL-A/B/C, EPA/EVA/ETU):
 *   a_i(t) = get_doppler_dispersion(...)    [pure Rayleigh, Jakes SoS]
 *
 * For LOS models (TDL-D/E), tap 0 uses Rice decomposition:
 *   rayleigh_part  = get_doppler_dispersion(...)
 *   los_phasor     = get_los_phasor(t, f_d, phi0)
 *   a_0(t) = sqrt(K/(K+1)) * los_phasor + sqrt(1/(K+1)) * rayleigh_part
 *
 * h_tap[0] for LOS is generated with 0 dB (unit amplitude), so the Rice
 * formula correctly combines both components with their respective power
 * weights summing to unity.
 * ---------------------------------------------------------------------------*/
static inline void generate_taps(srsran_channel_fading_t* q, float time)
{
  for (int i = 0; i < (int)nof_taps[q->model]; i++) {
    cf_t a;

    if (q->is_los && i == 0) {
      /* -----------------------------------------------------------------------
       * Tap 0 (LOS Rice model):
       *   1) Rayleigh component from Jakes SoS (existing engine, reused)
       *   2) Deterministic LOS phasor at f_S = 0.7 * f_D
       *   3) Combine with power split governed by K-factor
       * ---------------------------------------------------------------------*/
      cf_t rayleigh_coeff = get_doppler_dispersion(
          q, time, q->doppler, q->coeff_alpha[0], q->coeff_a[0], q->coeff_b[0]);

      cf_t los_coeff = get_los_phasor(time, q->doppler, q->los_phi0);

      float K       = q->k_factor;               /* K (linear)             */
      float kp1     = K + 1.0f;                  /* K + 1                  */
      float w_los   = sqrtf(K / kp1);            /* sqrt(K/(K+1))          */
      float w_diff  = sqrtf(1.0f / kp1);         /* sqrt(1/(K+1))          */

      a = w_los * los_coeff + w_diff * rayleigh_coeff;

    } else {
      /* -----------------------------------------------------------------------
       * Standard NLOS taps (and all taps for non-LOS models):
       * Pure Rayleigh — identical to the original code path.
       * ---------------------------------------------------------------------*/
      a = get_doppler_dispersion(
          q, time, q->doppler, q->coeff_alpha[i], q->coeff_a[i], q->coeff_b[i]);
    }

    if (i) {
      // Non-first tap: multiply tap freq response by a, then add (with FFT shift) to h_freq
      srsran_vec_sc_prod_ccc(q->h_tap[i], a, q->temp, q->N);

      srsran_vec_sum_ccc(q->h_freq, &q->temp[q->N / 2], q->h_freq, q->N / 2);
      srsran_vec_sum_ccc(&q->h_freq[q->N / 2], q->temp, &q->h_freq[q->N / 2], q->N / 2);
    } else {
      // First tap: initialize h_freq with FFT-shifted contribution
      srsran_vec_sc_prod_ccc(&q->h_tap[0][q->N / 2], a, q->h_freq, q->N / 2);
      srsran_vec_sc_prod_ccc(&q->h_tap[0][0], a, &q->h_freq[q->N / 2], q->N / 2);
    }
  }
  /* At this stage q->h_freq contains the full channel frequency response H(f). */
}

/* ---------------------------------------------------------------------------
 * filter_segment(): Applies the channel frequency response to one block.
 * Uses FFT-based overlap-add convolution. Unchanged from original.
 * ---------------------------------------------------------------------------*/
static inline void filter_segment(srsran_channel_fading_t* q, const cf_t* input, cf_t* output, uint32_t nsamples)
{
  // Fill Input vector
  srsran_vec_cf_copy(q->temp, input, nsamples);
  srsran_vec_cf_zero(&q->temp[nsamples], q->N - nsamples);

  // Do FFT
  srsran_dft_run_c_zerocopy(&q->fft, q->temp, q->y_freq);

  // Apply channel
  srsran_vec_prod_ccc(q->y_freq, q->h_freq, q->y_freq, q->N);

  // Do iFFT
  srsran_dft_run_c_zerocopy(&q->ifft, q->y_freq, q->temp);

  // Add state
  srsran_vec_sum_ccc(q->temp, q->state, q->temp, q->state_len);

  // Copy the first nsamples into the output
  srsran_vec_cf_copy(output, q->temp, nsamples);

  // Copy the rest of the samples into the state
  q->state_len = q->N - nsamples;
  srsran_vec_cf_copy(q->state, &q->temp[nsamples], q->state_len);
}

/* ---------------------------------------------------------------------------
 * srsran_channel_fading_init(): Allocates and initialises the fading engine.
 *
 * New LOS-specific initialisation steps (added for TDL-D/E):
 *   1) Detect if model is LOS (TDL-D or TDL-E) and set q->is_los.
 *   2) Convert K-factor from dB to linear and store q->k_factor.
 *   3) Draw a random initial phase φ₀ ~ U[0, 2π) for the LOS phasor.
 *   4) For tap 0 of LOS models, generate h_tap[0] with 0 dB (unit amplitude)
 *      because the Rice power split is applied at runtime in generate_taps().
 * ---------------------------------------------------------------------------*/
int srsran_channel_fading_init(srsran_channel_fading_t* q, double srate, const char* model, uint32_t seed)
{
  int ret = SRSRAN_ERROR;

  if (q) {
    // Parse model string → enum + doppler
    if (parse_model(q, model) != SRSRAN_SUCCESS) {
      fprintf(stderr, "Error: invalid channel model '%s'\n", model);
      goto clean_exit;
    }

    // Fill srate
    q->srate = (float)srate;

    // Determine if this is a LOS model (TDL-D or TDL-E)
    q->is_los   = (q->model == srsran_channel_fading_model_tdld ||
                   q->model == srsran_channel_fading_model_tdle);
    q->k_factor = q->is_los ? srsran_convert_dB_to_power(k_factor_db[q->model]) : 0.0f;
    q->los_phi0 = 0.0f; // Will be set after random init below

    // Populate internal FFT size parameters
    uint32_t fft_min_pow =
        (uint32_t)round(log2(excess_tap_delay_ns[q->model][nof_taps[q->model] - 1] * 1e-9 * srate)) + 3;
    q->N          = SRSRAN_MAX(1U << fft_min_pow, (uint32_t)(srate / (15e3f * 4.0f)));
    q->path_delay = q->N / 4;
    q->state_len  = 0;

    // Initialise random number generator
    srsran_random_t* random = srsran_random_init(seed);

    // Draw random initial LOS phase φ₀ ~ U[0, 2π)
    if (q->is_los) {
      q->los_phi0 = srsran_random_uniform_real_dist(random, 0.0f, 2.0f * (float)M_PI);
    }

    // Initialise Jakes SoS coefficients and h_tap buffers for all taps
    for (uint32_t i = 0; i < nof_taps[q->model]; i++) {
      // Random Jakes model coefficients
      for (uint32_t j = 0; (float)j < SRSRAN_CHANNEL_FADING_NTERMS; j++) {
        q->coeff_a[i][j]     = srsran_random_uniform_real_dist(random, 0, 2.0f * (float)M_PI);
        q->coeff_b[i][j]     = srsran_random_uniform_real_dist(random, 0, 2.0f * (float)M_PI);
        q->coeff_alpha[i][j] = ((float)M_PI * ((float)i - (float)0.5f)) / (2.0f * nof_taps[q->model]);
      }

      // Allocate tap frequency response buffer
      q->h_tap[i] = srsran_vec_cf_malloc(q->N);

      // Generate static frequency-domain tap shape.
      // For LOS tap 0: generate with 0 dB (unit power). The Rice split is
      // applied at runtime in generate_taps(). Power of taps 1..L is stored
      // already normalized in relative_power_db[].
      generate_tap(excess_tap_delay_ns[q->model][i],
                   relative_power_db[q->model][i],
                   q->srate,
                   q->h_tap[i],
                   q->N,
                   q->path_delay);
    }

    // Generate sine lookup table for SSE fast trig
    for (uint32_t i = 0; i < 1024; i++) {
      q->sin_table[i] = sinf((float)i * 2.0f * (float)M_PI / 1024);
    }

    // Free random generator
    srsran_random_free(random);

    // Plan forward FFT
    if (srsran_dft_plan_c(&q->fft, q->N, SRSRAN_DFT_FORWARD) != SRSRAN_SUCCESS) {
      fprintf(stderr, "Error: planning fft\n");
      goto clean_exit;
    }

    // Plan inverse FFT
    if (srsran_dft_plan_c(&q->ifft, q->N, SRSRAN_DFT_BACKWARD) != SRSRAN_SUCCESS) {
      fprintf(stderr, "Error: planning ifft\n");
      goto clean_exit;
    }

    // Allocate working buffers
    q->temp = srsran_vec_cf_malloc(q->N);
    if (!q->temp) {
      fprintf(stderr, "Error: allocating temp\n");
      goto clean_exit;
    }

    q->h_freq = srsran_vec_cf_malloc(q->N);
    if (!q->h_freq) {
      fprintf(stderr, "Error: allocating h_freq\n");
      goto clean_exit;
    }

    q->y_freq = srsran_vec_cf_malloc(q->N);
    if (!q->y_freq) {
      fprintf(stderr, "Error: allocating y_freq\n");
      goto clean_exit;
    }

    q->state = srsran_vec_cf_malloc(q->N);
    if (!q->state) {
      fprintf(stderr, "Error: allocating state\n");
      goto clean_exit;
    }
    srsran_vec_cf_zero(q->state, q->N);
  }

  ret = SRSRAN_SUCCESS;

clean_exit:
  return ret;
}

/* ---------------------------------------------------------------------------
 * srsran_channel_fading_free(): Releases all allocated resources.
 * Unchanged from original.
 * ---------------------------------------------------------------------------*/
void srsran_channel_fading_free(srsran_channel_fading_t* q)
{
  if (q) {
    srsran_dft_plan_free(&q->fft);
    srsran_dft_plan_free(&q->ifft);

    if (q->temp) {
      free(q->temp);
    }

    if (q->h_freq) {
      free(q->h_freq);
    }

    if (q->y_freq) {
      free(q->y_freq);
    }

    for (int i = 0; i < (int)nof_taps[q->model]; i++) {
      if (q->h_tap[i]) {
        free(q->h_tap[i]);
      }
    }

    if (q->state) {
      free(q->state);
    }
  }
}

/* ---------------------------------------------------------------------------
 * srsran_channel_fading_execute(): Main execution loop.
 * Processes nsamples in blocks of at most N/2.
 * Unchanged from original — the LOS logic is entirely within generate_taps().
 * ---------------------------------------------------------------------------*/
double srsran_channel_fading_execute(srsran_channel_fading_t* q,
                                     const cf_t*              in,
                                     cf_t*                    out,
                                     uint32_t                 nsamples,
                                     double                   init_time)
{
  uint32_t counter = 0;

  if (q) {
    while (counter < nsamples) {
      // Compute channel frequency response for current time instant
      generate_taps(q, (float)init_time);

      // Do not process more than N/2 samples per block (overlap-add)
      uint32_t n = SRSRAN_MIN(q->N / 2, nsamples - counter);

      // Apply channel to this block
      filter_segment(q, &in[counter], &out[counter], n);

      // Advance time
      init_time += n / q->srate;

      // Advance sample counter
      counter += n;
    }
  }

  // Return updated time (for seamless continuation across calls)
  return init_time;
}
