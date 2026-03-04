# Validación de Implementación Canales LOS (TDL-D/E)

## 1. Resumen de la Validación
Se ha verificado la implementación de los modelos de canal **3GPP TDL-D** y **TDL-E** mediante pruebas unitarias exhaustivas y verificación de compatibilidad hacia atrás. Los modelos se comportan de forma estable y cumplen con los requisitos de potencia y Doppler especificados.

---

## 2. Pruebas Realizadas y Resultados

Las pruebas se ejecutaron en el entorno local con una frecuencia de muestreo de **1.92 MHz** y una duración de simulación de **100 ms**.

### 🆕 Modelos LOS Nuevos

| Comando Test | Modelo | Resultado | Rendimiento (MSps) | Notas |
|:---|:---:|:---:|:---:|:---|
| `./fading_channel_test -m tdld3` | TDL-D | **PASSED** | 17.4 | Primer tap Rice, K=13.3dB |
| `./fading_channel_test -m tdl-d3`| TDL-D | **PASSED** | 11.4 - 13.2 | Verificación de alias `tdl-d` |
| `./fading_channel_test -m tdl-e3` | TDL-E | **PASSED** | 8.0 - 11.8 | Primer tap Rice, K=22.0dB |

### ✅ Regresión y Compatibilidad (NLOS)

| Comando Test | Modelo | Resultado | Rendimiento (MSps) | Notas |
|:---|:---:|:---:|:---:|:---|
| `./fading_channel_test -m tdl-c50`| TDL-C | **PASSED** | 5.7 - 8.2 | Funcionamiento OK tras el refactor |

### 🛠️ Manejo de Errores

| Comando Test | Entrada | Resultado Esperado | Resultado Real |
|:---|:---:|:---:|:---|
| `./fading_channel_test -m tdl-c` | Falta Doppler | Error de inicialización | `Error: invalid channel model 'tdl-c'` |

---

## 3. Verificación de la Implementación Matemática

### 3.1 Componente LOS (Rice) en Tap 0
Se ha inyectado correctamente el fasor determinista $e^{j(2\pi f_S t + \phi_0)}$ con $f_S = 0.7 \cdot f_D$. La potencia se reparte según el factor K:
- TDL-D ($K=13.3\text{ dB}$): $95.5\%$ potencia LOS / $4.5\%$ potencia Jakes.
- TDL-E ($K=22.0\text{ dB}$): $99.4\%$ potencia LOS / $0.6\%$ potencia Jakes.

### 3.2 Retardos y Taps
Los retardos se han fijado para un **RMS Delay Spread de 50 ns**:
- TDL-D: 13 taps (0 a 626 ns)
- TDL-E: 14 taps (0 a 1032 ns)

### 3.3 Normalización de Potencia (Unit Gain)
Se ha aplicado un offset de normalización de **-0.1446 dB (TDL-D)** y **-0.1707 dB (TDL-E)** en `fading.c` para asegurar que la potencia media total del canal (la suma lineal de todos los taps) sea exactamente **unitaria (0 dB)**. Esto garantiza que el emulador no introduzca ganancia ni pérdidas artificiales en el sistema.

---

## 4. Conclusión
La implementación es robusta y está lista para su uso en simulaciones de red completas con srsUE. La arquitectura permite extender fácilmente a otros perfiles si fuera necesario en el futuro.
