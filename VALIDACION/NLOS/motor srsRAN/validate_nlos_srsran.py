import numpy as np
import matplotlib.pyplot as plt
import os

def main():
    # --- Especificaciones del sistema ---
    srate = 23.04e6  # sample rate 23.04 MHz
    num_iters = 10000
    
    # Calcular N (igual que fading.c)
    max_delay_ns = 482.9
    fft_min_pow = int(np.round(np.log2(max_delay_ns * 1e-9 * srate))) + 3
    N = int(max(1 << fft_min_pow, srate / 60000.0))
    path_delay_samples = N // 4

    # --- Taps TDL-A (Teóricos) ---
    delays_ns = np.array([
        0.0, 19.1, 20.1, 29.3, 23.1, 26.9, 33.5, 28.8, 38.1, 76.9, 
        94.9, 111.2, 108.6, 124.7, 125.6, 152.9, 204.1, 222.9, 228.5, 
        239.8, 250.3, 265.2, 482.9
    ])
    powers_db = np.array([
        -18.8004, -5.4004, -7.6004, -9.4004, -11.4004, -13.6004, 
        -15.3004, -15.9004, -12.9004, -21.3004, -12.0004, -22.1004, 
        -17.8004, -20.6004, -16.2004, -16.7004, -18.1004, -21.6004, 
        -23.7004, -24.3004, -22.0004, -25.3004, -35.1004
    ])
    P_linear = 10 ** (powers_db / 10.0)
    delays_samples = delays_ns * 1e-9 * srate

    # --- Lectura de los datos Nativos de srsRAN (generados por C) ---
    if not os.path.exists("c_pdp_linear.bin") or not os.path.exists("c_y_history.bin"):
        print("Error: No se encontraron los archivos binarios de C.")
        return
        
    pdp_linear = np.fromfile("c_pdp_linear.bin", dtype=np.float32)
    y_history_abs = np.fromfile("c_y_history.bin", dtype=np.float32).reshape((num_iters, N))
    
    energia_nativa = np.sum(pdp_linear)
    energia_teorica = np.sum(P_linear)
    
    pdp_est_db = 10.0 * np.log10(pdp_linear + 1e-12)
    energia_estimada = energia_nativa
    
    # --- Análisis de Fluctuación (Histograma Rayleigh) ---
    n_peak = np.argmax(pdp_linear)
    amplitudes_peak = y_history_abs[:, n_peak]
    
    # Parámetro sigma de Rayleigh
    sigma = np.sqrt(pdp_linear[n_peak] / 2.0)
    x_rayleigh = np.linspace(0, np.max(amplitudes_peak), 200)
    pdf_rayleigh = (x_rayleigh / sigma**2) * np.exp(-(x_rayleigh**2) / (2 * sigma**2))

    # --- Graficación Conjunta ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Subplot 1: PDP
    ax1.plot(np.arange(N), pdp_est_db, marker='o', linestyle='-', color='tab:green', 
             linewidth=1.5, markersize=5, label='PDP srsRAN Nativo C')
             
    n_teoricos = delays_samples + path_delay_samples
    ax1.stem(n_teoricos, powers_db, linefmt='r--', markerfmt='rX', 
             basefmt=" ", label='Taps Teóricos (3GPP TDL-A)')

    textos = (
        "Validación de Energía:\n"
        f"∑ Taps Teóricos = {energia_teorica:.4f}\n"
        f"∑ PDP srsRAN = {energia_estimada:.4f}"
    )
    ax1.text(0.02, 0.05, textos, transform=ax1.transAxes, fontsize=11,
             bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'))

    min_xlim = path_delay_samples - 5
    max_xlim = np.max(n_teoricos) + 5
    ax1.set_xlim(min_xlim, max_xlim)
    ax1.set_ylim(-45, 5)
    
    ax1.set_title('Validación PDP - TDL-A NLOS (srsRAN Native f_c)', fontsize=14)
    ax1.set_xlabel('Índice de Muestra $n$', fontsize=12)
    ax1.set_ylabel('Potencia (dB)', fontsize=12)
    ax1.legend(loc='upper right', fontsize=11)
    ax1.grid(True, which='major', linestyle='-', alpha=0.5)
    ax1.minorticks_on()
    ax1.grid(True, which='minor', linestyle=':', alpha=0.3)
    
    # Subplot 2: Histograma Rayleigh
    ax2.hist(amplitudes_peak, bins=50, density=True, alpha=0.6, color='tab:green', edgecolor='black', 
             label=f'Histog. srsRAN $|y[n={n_peak}, m]|$')
    ax2.plot(x_rayleigh, pdf_rayleigh, 'r-', linewidth=2, label='PDF Teórica Rayleigh')
    
    ax2.set_title(f'Validación Fluctuación Rayleigh nativa (Muestra $n={n_peak}$)', fontsize=14)
    ax2.set_xlabel('Amplitud $|y|$', fontsize=12)
    ax2.set_ylabel('Densidad de Probabilidad', fontsize=12)
    ax2.legend(loc='upper right', fontsize=11)
    ax2.grid(True, linestyle='--', alpha=0.6)

    fig.tight_layout()
    output_filename = 'Validacion_NLOS_TDLA_srsRAN.png'
    fig.savefig(output_filename, dpi=300, format='png', bbox_inches='tight')
    print(f"Grafica combinada guardada en {output_filename}")

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    main()
