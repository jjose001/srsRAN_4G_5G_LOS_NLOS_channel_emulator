import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import rice

def main():
    # --- Especificaciones del sistema ---
    srate = 23.04e6  # sample rate 23.04 MHz
    num_iters = 100000
    
    # Calcular N (igual que fading.c)
    max_delay_ns = 626.25
    fft_min_pow = int(np.round(np.log2(max_delay_ns * 1e-9 * srate))) + 3
    N = int(max(1 << fft_min_pow, srate / 60000.0))
    path_delay_samples = N // 4

    # --- Taps TDL-D (LOS, DS=50ns) ---
    delays_ns = np.array([
        0.00, 1.75, 30.60, 68.15, 70.25, 90.20, 129.80, 88.75, 202.10, 396.85, 471.20, 485.40, 626.25
    ])
    # Utilizamos las tablas exactas precalculadas nativas de srsRAN_4G/lib/src/phy/channel/fading.c
    powers_db = np.array([
        -0.1446, -18.9446, -21.1446, -22.9446, -18.0446, -20.2446, -22.0446, -23.0446, -27.9446, -23.7446, -24.9446, -30.1446, -27.8446
    ])
    P_linear = 10 ** (powers_db / 10.0)
    delays_samples = delays_ns * 1e-9 * srate

    # --- Lectura de los datos Nativos de srsRAN (generados por C) ---
    if not os.path.exists("c_pdp_linear.bin") or not os.path.exists("c_y_history.bin"):
        print("Error: No se encontraron los archivos binarios de C.")
        return
        
    pdp_linear = np.fromfile("c_pdp_linear.bin", dtype=np.float32)
    # Normalización unitaria (1.0) sugerida por el usuario
    pdp_linear = pdp_linear / np.sum(pdp_linear)
    y_history_abs = np.fromfile("c_y_history.bin", dtype=np.float32).reshape((num_iters, N))
    # Normalizamos también las muestras individuales por el mismo factor para mantener consistencia
    total_energy_raw = np.fromfile("c_pdp_linear.bin", dtype=np.float32).sum()
    y_history_abs = y_history_abs / np.sqrt(total_energy_raw)
    
    energia_nativa = np.sum(pdp_linear)
    # Taps teóricos normalizados a 1.0 también para la línea de referencia
    P_linear_norm = P_linear / np.sum(P_linear)
    powers_db_norm = 10 * np.log10(P_linear_norm + 1e-12)
    
    pdp_est_db = 10.0 * np.log10(pdp_linear + 1e-12)
    energia_estimada = energia_nativa
    
    # --- Análisis de Fluctuación (Histograma Rice) ---
    n_peak = np.argmax(pdp_linear)
    amplitudes_peak = y_history_abs[:, n_peak]
    
    # Para el primer tap en TDL-D (Ricean tap), la potencia total de la muestra es P_peak.
    # El K-factor es K = 13.3 dB
    K_db = 13.3
    K_linear = 10**(K_db / 10.0)
    Potencia_total_tap = pdp_linear[n_peak]
    
    # La potencia se divide en LOS (ν^2) y NLOS dispersión (2σ^2)
    # K = ν^2 / (2σ^2)  y Potencia_total_tap = ν^2 + 2σ^2
    # ν^2 = Potencia_total_tap * K / (K + 1)
    # 2σ^2 = Potencia_total_tap / (K + 1)
    
    nu = np.sqrt(Potencia_total_tap * K_linear / (K_linear + 1.0))
    sigma = np.sqrt((Potencia_total_tap / (K_linear + 1.0)) / 2.0)
    
    # scipy.stats.rice toma parametro b = ν/σ y la escala de la distribución es σ.
    b_param = nu / sigma
    x_rice = np.linspace(0, np.max(amplitudes_peak)*1.2, 200)
    pdf_rice = rice.pdf(x_rice, b_param, scale=sigma)

    # --- Graficación Conjunta ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Subplot 1: PDP
    ax1.plot(np.arange(N), pdp_est_db, marker='o', linestyle='-', color='tab:orange', 
             linewidth=1.5, markersize=5, label='PDP srsRAN Nativo C')
             
    n_teoricos = delays_samples + path_delay_samples
    ax1.stem(n_teoricos, powers_db_norm, linefmt='b--', markerfmt='bX', 
             basefmt=" ", label='Taps Teóricos (3GPP Normalizados)')

    textos = (
        "Energía:\n"
        f"∑ Taps Teóricos = {np.sum(P_linear_norm):.4f}\n"
        f"∑ PDP srsRAN = {energia_estimada:.4f}"
    )
    ax1.text(0.02, 0.05, textos, transform=ax1.transAxes, fontsize=11,
             bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'))

    min_xlim = path_delay_samples - 5
    max_xlim = np.max(n_teoricos) + 5
    ax1.set_xlim(min_xlim, max_xlim)
    ax1.set_ylim(-45, 5)
    
    ax1.set_title('Validación PDP - TDL-D LOS (srsRAN Native f_c)', fontsize=14)
    ax1.set_xlabel('Índice de Muestra $n$', fontsize=12)
    ax1.set_ylabel('Potencia (dB)', fontsize=12)
    ax1.legend(loc='upper right', fontsize=11)
    ax1.grid(True, which='major', linestyle='-', alpha=0.5)
    ax1.minorticks_on()
    ax1.grid(True, which='minor', linestyle=':', alpha=0.3)
    
    # Subplot 2: Histograma Rice
    ax2.hist(amplitudes_peak, bins=50, density=True, alpha=0.6, color='tab:orange', edgecolor='black', 
             label=f'Histog. srsRAN $|y[n={n_peak}, m]|$')
    ax2.plot(x_rice, pdf_rice, 'b-', linewidth=2, label='PDF Teórica Rice ($K=13.3$dB)')
    
    ax2.set_title(f'Validación Distribución Rice nativa (Muestra $n={n_peak}$)', fontsize=14)
    ax2.set_xlabel('Amplitud $|y|$', fontsize=12)
    ax2.set_ylabel('Densidad de Probabilidad', fontsize=12)
    ax2.legend(loc='upper right', fontsize=11)
    ax2.grid(True, linestyle='--', alpha=0.6)

    fig.tight_layout()
    output_filename = 'Validacion_LOS_TDLD_srsRAN.png'
    fig.savefig(output_filename, dpi=300, format='png', bbox_inches='tight')
    print(f"Grafica combinada guardada en {output_filename}")

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    main()
