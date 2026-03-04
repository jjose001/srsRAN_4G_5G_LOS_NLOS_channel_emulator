# Implementación de Canales LOS TDL-D / TDL-E en srsRAN_4G

**Fecha:** 2026-03-04  
**Repositorio:** `srsRAN_TDL_LOS_JOAO/srsRAN_4G`  
**Versión:** 1.0

---

## 1. Objetivo

Añadir soporte completo para los modelos de canal **TDL-D** y **TDL-E** (canales LOS según 3GPP TR 38.901) al emulador de canal interno de `srsRAN_4G`, manteniendo **compatibilidad total hacia atrás** con los modelos NLOS ya existentes (EPA/EVA/ETU/TDL-A/B/C).

---

## 2. Contexto y Motivación

Los modelos **TDL-D** y **TDL-E** representan escenarios de propagación con **Línea de Vista (LOS)** entre transmisor y receptor. A diferencia de los modelos NLOS (TDL-A/B/C), en LOS el **primer tap** sigue una distribución **Riciana** (mezcla de componente determinista + componente difusa Rayleigh), mientras que el resto de taps mantienen el comportamiento Rayleigh habitual.

### 2.1 Ecuación del Tap 0 (Rice)

$$a_0(t) = \sqrt{P_0} \left( \sqrt{\frac{K}{K+1}} \cdot e^{j(2\pi f_S t + \phi_0)} + \sqrt{\frac{1}{K+1}} \cdot g_0(t) \right)$$

- **K**: K-factor lineal del perfil (TDL-D: 13.3 dB = 21.38; TDL-E: 22.0 dB = 158.49)
- **f_S = 0.7 · f_D**: Frecuencia Doppler del pico LOS (según TR 38.901)
- **φ_0**: Fase inicial uniforme en [0, 2π), generada aleatoriamente en el arranque
- **g_0(t)**: Proceso Rayleigh (Jakes SoS, ya existente en el motor)
- **P_0 = 0 dB**: Potencia total normalizada del tap 0 (según nota de la tabla 3GPP)

### 2.2 Taps 1..L (Rayleigh, igual que modelos NLOS)

$$a_n(t) = \sqrt{P_n} \cdot g_n(t), \quad n \geq 1$$

---

## 3. Tablas de Parámetros (DS = 50 ns)

### TDL-D — 13 Taps Rayleigh + Tap 0 Rice (total 13 índices, tap[0] = Rice)

| Idx | Delay norm | Delay (ns) | Power (dB) | Tipo  |
|-----|-----------|------------|------------|-------|
| 0   | 0.000     | 0.00       | 0.0 total  | Rice (K=13.3 dB) |
| 1   | 0.035     | 1.75       | -18.8      | Rayleigh |
| 2   | 0.612     | 30.60      | -21.0      | Rayleigh |
| 3   | 1.363     | 68.15      | -22.8      | Rayleigh |
| 4   | 1.405     | 70.25      | -17.9      | Rayleigh |
| 5   | 1.804     | 90.20      | -20.1      | Rayleigh |
| 6   | 2.596     | 129.80     | -21.9      | Rayleigh |
| 7   | 1.775     | 88.75      | -22.9      | Rayleigh |
| 8   | 4.042     | 202.10     | -27.8      | Rayleigh |
| 9   | 7.937     | 396.85     | -23.6      | Rayleigh |
| 10  | 9.424     | 471.20     | -24.8      | Rayleigh |
| 11  | 9.708     | 485.40     | -30.0      | Rayleigh |
| 12  | 12.525    | 626.25     | -27.7      | Rayleigh |

**Nota:** El tap 0 tiene K₁ = 13.3 dB y potencia media de 0 dB.

### TDL-E — 14 Taps Rayleigh + Tap 0 Rice (total 14 índices, tap[0] = Rice)

| Idx | Delay norm | Delay (ns) | Power (dB) | Tipo  |
|-----|-----------|------------|------------|-------|
| 0   | 0.0000    | 0.00       | 0.0 total  | Rice (K=22.0 dB) |
| 1   | 0.5133    | 25.67      | -15.8      | Rayleigh |
| 2   | 0.5440    | 27.20      | -18.1      | Rayleigh |
| 3   | 0.5630    | 28.15      | -19.8      | Rayleigh |
| 4   | 0.5440    | 27.20      | -22.9      | Rayleigh |
| 5   | 0.7112    | 35.56      | -22.4      | Rayleigh |
| 6   | 1.9092    | 95.46      | -18.6      | Rayleigh |
| 7   | 1.9293    | 96.47      | -20.8      | Rayleigh |
| 8   | 1.9589    | 97.95      | -22.6      | Rayleigh |
| 9   | 2.6426    | 132.13     | -22.3      | Rayleigh |
| 10  | 3.7136    | 185.68     | -25.6      | Rayleigh |
| 11  | 5.4524    | 272.62     | -20.2      | Rayleigh |
| 12  | 12.0034   | 600.17     | -29.8      | Rayleigh |
| 13  | 20.6519   | 1032.60    | -29.2      | Rayleigh |

**Nota:** El tap 0 tiene K₁ = 22.0 dB y potencia media de 0 dB.

---

## 4. Ficheros Modificados

Solo se modifican **2 ficheros del motor de fading**. El resto del stack no requiere cambios gracias a la arquitectura existente.

### 4.1 `lib/include/srsran/phy/channel/fading.h`

**Cambios:**
1. Añadir `srsran_channel_fading_model_tdld` y `srsran_channel_fading_model_tdle` al enum.
2. Añadir 3 campos a la struct `srsran_channel_fading_t` para el estado LOS:
   - `bool is_los` — activa la rama Rice en tap 0
   - `float k_factor` — K lineal del tap 0
   - `float los_phi0` — Fase inicial aleatoria φ₀

### 4.2 `lib/src/phy/channel/fading.c`

**Cambios:**
1. Ampliar `nof_taps[9]` con 2 entradas (TDL-D: 13, TDL-E: 14).
2. Ampliar `excess_tap_delay_ns[9][...]` con las tablas TDL-D/E (DS=50 ns).
3. Ampliar `relative_power_db[9][...]` con las potencias de los taps Rayleigh. Para asegurar una **ganancia media total de 1.0 (0 dB)**, se ha aplicado un offset de normalización a cada tabla:
   - **TDL-D:** Offset de -0.1446 dB aplicado a todos los taps.
   - **TDL-E:** Offset de -0.1707 dB aplicado a todos los taps.
4. Añadir `parse_model()` entries para `tdld`/`tdl-d` y `tdle`/`tdl-e`.
5. Añadir helper `get_los_phasor(t, f_d, phi0) → cf_t`.
6. Modificar `generate_taps()`: para `q->is_los` en tap `i==0`, combinar componente LOS + Rayleigh escaladas.
7. Modificar `srsran_channel_fading_init()`: detectar modelo LOS, inicializar `is_los`, `k_factor`, `los_phi0`.

---

## 5. Arquitectura de la Solución (Reutilización Máxima)

```
fading_execute()
    └─► generate_taps(q, time)
            │
            ├─ tap i=0, is_los=true:
            │     rayleigh_coeff = get_doppler_dispersion(...)  ← REUTILIZADO
            │     los_coeff      = get_los_phasor(t, 0.7*fd, phi0)  ← NUEVO
            │     a = sqrt(P0) * (sqrt(K/(K+1))*los_coeff + sqrt(1/(K+1))*rayleigh_coeff)
            │
            └─ tap i>0 (o is_los=false para NLOS):
                  a = get_doppler_dispersion(...)  ← SIN CAMBIO
```

Todo el pipeline FFT/Overlap-Add de `filter_segment()` se reutiliza sin cambios.

---

## 6. Compatibilidad

| Modelo | Estado |
|--------|--------|
| `none` | ✅ Sin cambios |
| `epaX`, `evaX`, `etuX` | ✅ Sin cambios |
| `tdlaX`, `tdlbX`, `tdlcX` | ✅ Sin cambios |
| `tdldX` | 🆕 TDL-D LOS (Rice tap 0, K=13.3 dB) |
| `tdleX` | 🆕 TDL-E LOS (Rice tap 0, K=22.0 dB) |

Uso en `ue.conf`:
```ini
[channel.ul.fading]
enable = true
model  = tdld3       ; TDL-D con Doppler 3 Hz
; model = tdle5      ; TDL-E con Doppler 5 Hz
```

---

## 7. Estrategia de Verificación

1. **Compilación limpia** con `make srsue -j$(nproc)` desde el directorio `build/`.
2. **Test unitario** `fading_channel_test` con modelos `tdld` y `tdle`.
3. **Regresión NLOS** con `tdla`, `epa` para confirmar backward compatibility.
4. **Integración** en sesión srsUE con fading LOS activo.

Ver `validation.md` para resultados completos.
