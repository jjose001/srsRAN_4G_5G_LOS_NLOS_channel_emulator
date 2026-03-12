#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "srsran/phy/channel/fading.h"
#include "srsran/phy/utils/vector.h"

int main() {
    srsran_channel_fading_t q;
    float srate = 23040000;
    int num_iters = 10000;
    
    // Init TDL-A with a high Doppler (say 300Hz) to ensure fast fading and uncorrelated samples
    if(srsran_channel_fading_init(&q, srate, "tdla300", 1234) != SRSRAN_SUCCESS) {
        printf("Error init.\n");
        return -1;
    }
    
    int N = q.N;
    printf("Initialized TDL-A channel. N=%d\n", N);
    
    cf_t *in = srsran_vec_cf_malloc(N);
    cf_t *out = srsran_vec_cf_malloc(N);
    
    float *pdp_linear = (float*)calloc(N, sizeof(float));
    float *y_history_abs = (float*)malloc(num_iters * N * sizeof(float)); // We will store |y[n]|
    
    double init_time = 0;
    double time_step = 14e-3; // Jump in time by 14 ms to ensure COMPLETE uncorrelatedness (since coherence time ~ 1/fD ~ 1/300 = 3ms)
    
    for (int m = 0; m < num_iters; m++) {
        srsran_vec_cf_zero(in, N);
        
        // Emulate an impulse in frequency domain (X[k]=1 for all k).
        // The IFFT of this is just an impulse in time domain at index 0.
        // Oh wait, the `in` variable of `srsran_channel_fading_execute` is TIME DOMAIN. 
        // So I just send an impulse time-domain signal.
        // Wait, fading execute filters `in` (time domain) to `out` (time domain).
        // The python script did np.fft.ifft(X_k * H_k) = np.fft.ifft(1 * H_k) = h_n.
        // So taking the output of the channel for an impulse input is perfectly equivalent to h_n!
        // But the input impulse `in` has energy 1. Wait, fading_execute maintains a state.
        // Is passing `in[0] = 1.0` equivalent to all ones in X_k? Yes, because DFT(delta) = ones.
        
        // Wait, fading_execute uses `filter_segment` which adds the FFT of `in` with `h_freq`.
        // If we do `in[0] = 1.0`, its FFT is all 1/N or 1. Let's see how `srsran_dft_run_c_zerocopy` handles scaling.
        // Typically DFT is not scaled, so DFT(delta) = 1. Let's multiply `in[0]` by N, or just `in[0] = 1.0` and we'll normalize later.
        in[0] = 1.0f; 
        
        // We artificially advance time to make samples perfectly independent like the Python test
        init_time += time_step; 
        
        // In order to not carry over 'state' from previous simulation jump, we must clear state?
        // Yes, the fading channel keeps a state for overlap-add. Since we jumped in time, we should clear it.
        srsran_vec_cf_zero(q.state, N);
        q.state_len = 0;
        
        // Execute channel for N samples
        srsran_channel_fading_execute(&q, in, out, N, init_time);
        
        for (int n = 0; n < N; n++) {
            float mag2 = SRSRAN_CSQABS(out[n]); // raw power output
            if(!isnormal(mag2) && mag2 != 0.0f) {
                printf("Garbage at m=%d n=%d mag2=%f\n", m, n, mag2);
            }
            y_history_abs[m * N + n] = sqrtf(mag2);
            pdp_linear[n] += mag2;
        }
    }
    
    for (int n = 0; n < N; n++) {
        pdp_linear[n] /= num_iters;
    }
    
    // Calculate total energy
    float energia_estimada = 0;
    for(int n=0; n<N; n++) energia_estimada += pdp_linear[n];
    printf("Energia Total C: %.6f\n", energia_estimada);
    
    // Save to binary files for python to plot
    FILE *fp = fopen("c_pdp_linear.bin", "wb");
    fwrite(pdp_linear, sizeof(float), N, fp);
    fclose(fp);
    
    fp = fopen("c_y_history.bin", "wb");
    fwrite(y_history_abs, sizeof(float), num_iters * N, fp);
    fclose(fp);
    
    free(pdp_linear);
    free(y_history_abs);
    free(in);
    free(out);
    srsran_channel_fading_free(&q);
    return 0;
}
