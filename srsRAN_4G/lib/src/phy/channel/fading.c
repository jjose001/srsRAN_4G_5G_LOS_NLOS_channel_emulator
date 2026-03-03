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
 * TDL models from 3GPP TR 38.901 Table 7.7.2-1, 7.7.2-2, 7.7.2-3
 */
const static uint32_t nof_taps[7] = {1, 7, 9, 9, 23, 23, 24};

const static float excess_tap_delay_ns[7][SRSRAN_CHANNEL_FADING_MAXTAPS] = {
    /* None */ {0, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
    /* EPA  */ {0, 30, 70, 90, 110, 190, 410, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
    /* EVA  */ {0, 30, 150, 310, 370, 710, 1090, 1730, 2510, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
    /* ETU  */ {0, 50, 120, 200, 230, 500, 1600, 2300, 5000, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
    /* TDL-A: normalized delays * 50 ns (3GPP TR 38.901 Table 7.7.2-1) */
    {0.0f, 19.1f, 20.1f, 29.3f, 23.1f, 26.9f, 33.5f, 28.8f, 38.1f, 76.9f, 94.9f, 111.2f, 108.6f, 124.7f, 125.6f, 152.9f, 204.1f, 222.9f, 228.5f, 239.8f, 250.3f, 265.2f, 482.9f, NAN},
    /* TDL-B: normalized delays * 50 ns (3GPP TR 38.901 Table 7.7.2-2) */
    {0.0f, 5.4f, 10.8f, 10.5f, 14.3f, 14.9f, 18.8f, 25.3f, 18.4f, 18.5f, 28.5f, 26.4f, 55.1f, 63.8f, 77.4f, 89.2f, 100.8f, 141.5f, 151.1f, 180.9f, 205.3f, 213.9f, 239.2f, NAN},
    /* TDL-C: normalized delays * 50 ns (3GPP TR 38.901 Table 7.7.2-3) */
    {0.0f, 10.5f, 11.1f, 11.6f, 10.9f, 31.8f, 32.2f, 32.8f, 32.9f, 39.7f, 41.1f, 46.7f, 61.4f, 65.4f, 108.5f, 135.5f, 212.9f, 230.0f, 274.5f, 280.4f, 315.3f, 331.9f, 352.1f, 432.6f}};

const static float relative_power_db[7][SRSRAN_CHANNEL_FADING_MAXTAPS] = {
    /* None */ {+0.0f, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
    /* EPA  */ {+0.0f, -1.0f, -2.0f, -3.0f, -8.0f, -17.2f, -20.8f, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
    /* EVA  */ {+0.0f, -1.5f, -1.4f, -3.6f, -0.6f, -9.1f, -7.0f, -12.0f, -16.9f, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
    /* ETU  */ {-1.0f, -1.0f, -1.0f, +0.0f, +0.0f, +0.0f, -3.0f, -5.0f, -7.0f, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
    /* TDL-A: 3GPP TR 38.901 Table 7.7.2-1 (all taps use Rayleigh fading) - NORMALIZED: sum of linear powers = 1 (0 dB total gain) */
    {-18.8004f, -5.4004f, -7.6004f, -9.4004f, -11.4004f, -13.6004f, -15.3004f, -15.9004f, -12.9004f, -21.3004f, -12.0004f, -22.1004f, -17.8004f, -20.6004f, -16.2004f, -16.7004f, -18.1004f, -21.6004f, -23.7004f, -24.3004f, -22.0004f, -25.3004f, -35.1004f, NAN},
    /* TDL-B: 3GPP TR 38.901 Table 7.7.2-2 (all taps use Rayleigh fading) - NORMALIZED: sum of linear powers = 1 (0 dB total gain) */
    {-8.5083f, -10.7083f, -12.5083f, -11.7083f, -18.3083f, -9.7083f, -11.9083f, -13.7083f, -16.1083f, -11.5083f, -17.4083f, -17.5083f, -13.3083f, -14.2083f, -16.0083f, -10.4083f, -16.1083f, -20.7083f, -18.3083f, -19.9083f, -23.4083f, -17.7083f, -19.8083f, NAN},
    /* TDL-C: 3GPP TR 38.901 Table 7.7.2-3 (all taps use Rayleigh fading) - NORMALIZED: sum of linear powers = 1 (0 dB total gain) */
    {-12.0897f, -8.8897f, -11.1897f, -12.8897f, -10.1897f, -7.6897f, -9.8897f, -11.5897f, -15.0897f, -14.7897f, -18.3897f, -18.7897f, -12.7897f, -14.4897f, -16.3897f, -20.8897f, -21.5897f, -21.5897f, -23.4897f, -24.7897f, -23.6897f, -23.3897f, -29.2897f, -30.4897f}};

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
  } else {
    ret = SRSRAN_ERROR;
  }

  if (ret == SRSRAN_SUCCESS) {
    if (strlen(str) > offset) {
      // Parse doppler value from string (e.g., "epa5" -> 5, "tdla70" -> 70)
      q->doppler = (float)strtod(&str[offset], NULL);
      if (isnan(q->doppler) || isinf(q->doppler) || q->doppler < 0.0f) {
        q->doppler = 0.0f;
      }
    } else {
      // EPA, EVA, ETU, TDL-A, TDL-B, TDL-C all require doppler value
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

static inline void generate_tap(float delay_ns, float power_db, float srate, cf_t* buf, uint32_t N, uint32_t path_delay)
{
  float amplitude = srsran_convert_dB_to_power(power_db);
  float O         = (delay_ns * 1e-9f * srate + path_delay) / (float)N;
  cf_t  a0        = amplitude / N;

  srsran_vec_gen_sine(a0, -O, buf, N);
}

static inline void generate_taps(srsran_channel_fading_t* q, float time)
{
  // Generate taps
  for (int i = 0; i < nof_taps[q->model]; i++) {
    // Compute phase for the doppler dispersion
    cf_t a = get_doppler_dispersion(q, time, q->doppler, q->coeff_alpha[i], q->coeff_a[i], q->coeff_b[i]);

    if (i) {
      // Copy tap frequency response
      srsran_vec_sc_prod_ccc(q->h_tap[i], a, q->temp, q->N);

      // Add to frequency response, shifts FFT at same time
      srsran_vec_sum_ccc(q->h_freq, &q->temp[q->N / 2], q->h_freq, q->N / 2);
      srsran_vec_sum_ccc(&q->h_freq[q->N / 2], q->temp, &q->h_freq[q->N / 2], q->N / 2);
    } else {
      // Copy tap frequency response
      srsran_vec_sc_prod_ccc(&q->h_tap[i][q->N / 2], a, q->h_freq, q->N / 2);
      srsran_vec_sc_prod_ccc(&q->h_tap[i][0], a, &q->h_freq[q->N / 2], q->N / 2);
    }
  }
  // at this stage, q->h_freq should contain the frequency response
}

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

int srsran_channel_fading_init(srsran_channel_fading_t* q, double srate, const char* model, uint32_t seed)
{
  int ret = SRSRAN_ERROR;

  if (q) {
    // Parse model
    if (parse_model(q, model) != SRSRAN_SUCCESS) {
      fprintf(stderr, "Error: invalid channel model '%s'\n", model);
      goto clean_exit;
    }

    // Fill srate
    q->srate = (float)srate;

    // Populate internal parameters
    uint32_t fft_min_pow =
        (uint32_t)round(log2(excess_tap_delay_ns[q->model][nof_taps[q->model] - 1] * 1e-9 * srate)) + 3;
    q->N          = SRSRAN_MAX(1U << fft_min_pow, (uint32_t)(srate / (15e3f * 4.0f)));
    q->path_delay = q->N / 4;
    q->state_len  = 0;

    // Initialise random number
    srsran_random_t* random = srsran_random_init(seed);

    // Initialise values for each tap
    for (uint32_t i = 0; i < nof_taps[q->model]; i++) {
      // Random Jakes model Coeffients
      for (uint32_t j = 0; (float)j < SRSRAN_CHANNEL_FADING_NTERMS; j++) {
        q->coeff_a[i][j]     = srsran_random_uniform_real_dist(random, 0, 2.0f * (float)M_PI);
        q->coeff_b[i][j]     = srsran_random_uniform_real_dist(random, 0, 2.0f * (float)M_PI);
        q->coeff_alpha[i][j] = ((float)M_PI * ((float)i - (float)0.5f)) / (2.0f * nof_taps[q->model]);
      }

      // Allocate tap frequency response
      q->h_tap[i] = srsran_vec_cf_malloc(q->N);

      // Generate tap frequency response
      generate_tap(
          excess_tap_delay_ns[q->model][i], relative_power_db[q->model][i], q->srate, q->h_tap[i], q->N, q->path_delay);
    }

    // Generate sine Table
    for (uint32_t i = 0; i < 1024; i++) {
      q->sin_table[i] = sinf((float)i * 2.0f * (float)M_PI / 1024);
    }

    // Free random
    srsran_random_free(random);

    // Plan FFT
    if (srsran_dft_plan_c(&q->fft, q->N, SRSRAN_DFT_FORWARD) != SRSRAN_SUCCESS) {
      fprintf(stderr, "Error: planning fft\n");
      goto clean_exit;
    }

    // Plan iFFT
    if (srsran_dft_plan_c(&q->ifft, q->N, SRSRAN_DFT_BACKWARD) != SRSRAN_SUCCESS) {
      fprintf(stderr, "Error: planning ifft\n");
      goto clean_exit;
    }

    // Allocate memory
    q->temp = srsran_vec_cf_malloc(q->N);
    if (!q->temp) {
      fprintf(stderr, "Error: allocating h_freq\n");
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
      fprintf(stderr, "Error: allocating y_freq\n");
      goto clean_exit;
    }
    srsran_vec_cf_zero(q->state, q->N);
  }

  ret = SRSRAN_SUCCESS;

clean_exit:
  return ret;
}

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

    for (int i = 0; i < nof_taps[q->model]; i++) {
      if (q->h_tap[i]) {
        free(q->h_tap[i]);
      }
    }

    if (q->state) {
      free(q->state);
    }
  }
}

double srsran_channel_fading_execute(srsran_channel_fading_t* q,
                                     const cf_t*              in,
                                     cf_t*                    out,
                                     uint32_t                 nsamples,
                                     double                   init_time)
{
  uint32_t counter = 0;

  if (q) {
    while (counter < nsamples) {
      // Generate taps
      generate_taps(q, (float)init_time);

      // Do not process more than N/2 samples
      uint32_t n = SRSRAN_MIN(q->N / 2, nsamples - counter);

      // Execute
      filter_segment(q, &in[counter], &out[counter], n);

      // Increment time
      init_time += n / q->srate;

      // Increment counter
      counter += n;
    }
  }

  // Return time
  return init_time;
}
