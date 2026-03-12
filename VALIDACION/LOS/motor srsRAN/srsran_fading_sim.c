#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "srsran/phy/channel/fading.h"
#include "srsran/phy/utils/vector.h"

int main() {
    srsran_channel_fading_t q;
    float srate = 23040000;
    const int num_iters = 100000;
    
    // Init TDL-D LOS with Doppler 100Hz
    if(srsran_channel_fading_init(&q, srate, "tdld100", 1234) != SRSRAN_SUCCESS) {
        printf("Error init.\n");
        return -1;
    }
    
    int N = q.N;
    printf("Initialized TDL-D channel. N=%d\n", N);
    
    cf_t *in = srsran_vec_cf_malloc(N);
    cf_t *out = srsran_vec_cf_malloc(N);
    
    float *pdp_linear = (float*)calloc(N, sizeof(float));
    float *y_history_abs = (float*)malloc((size_t)num_iters * N * sizeof(float));
    
    if (!y_history_abs || !pdp_linear || !in || !out) {
        printf("Error: Memory allocation failed.\n");
        return -1;
    }
    
    double current_time = 0;
    double time_step = 20e-3; // Increase time step Slightly for better decorrelation
    
    for (int m = 0; m < num_iters; m++) {
        srsran_vec_cf_zero(in, N);
        in[0] = 1.0f; // Input impulse
        
        // Advance time for this realization
        current_time += time_step;
        
        // Reset state carefully to avoid tail leakage but keep internal consistency
        srsran_vec_cf_zero(q.state, N);
        q.state_len = 0;
        
        srsran_channel_fading_execute(&q, in, out, N, current_time);
        
        for (int n = 0; n < N; n++) {
            float i = __real__ out[n];
            float q_val = __imag__ out[n];
            float mag2 = i*i + q_val*q_val;
            
            // Protection against numerical stability issues
            if (!isfinite(mag2)) mag2 = 0.0f;
            
            y_history_abs[(size_t)m * N + n] = sqrtf(mag2);
            pdp_linear[n] += mag2;
        }
        
        if (m % 20000 == 0) printf("Progress: %d%%\n", (m*100)/num_iters);
    }
    
    // Average power
    for (int n = 0; n < N; n++) {
        pdp_linear[n] /= num_iters;
    }
    
    float total_pwr = 0;
    for (int n = 0; n < N; n++) total_pwr += pdp_linear[n];
    printf("Energia Total C (Promediada): %.6f\n", total_pwr);
    
    // Write results
    FILE *fp = fopen("c_pdp_linear.bin", "wb");
    fwrite(pdp_linear, sizeof(float), N, fp);
    fclose(fp);
    
    fp = fopen("c_y_history.bin", "wb");
    fwrite(y_history_abs, sizeof(float), (size_t)num_iters * N, fp);
    fclose(fp);
    
    printf("Simulacion completada. Datos guardados.\n");
    
    free(pdp_linear);
    free(y_history_abs);
    free(in);
    free(out);
    srsran_channel_fading_free(&q);
    
    return 0;
}
