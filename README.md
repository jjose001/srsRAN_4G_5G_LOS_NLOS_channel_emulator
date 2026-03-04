# Implementation of LOS TDL-D / TDL-E Channels in srsRAN_4G

**Date:** 2026-03-04  
**Repository:** `srsRAN_TDL_LOS_JOAO/srsRAN_4G`  
**Version:** 1.0

---

## 1. Objective

Add full support for the **TDL-D** and **TDL-E** channel models (**LOS channels** as per 3GPP TR 38.901) to the internal channel emulator in `srsRAN_4G`, while maintaining **full backwards compatibility** with the existing NLOS models (EPA/EVA/ETU/TDL-A/B/C).

---

## 2. Context and Motivation

The **TDL-D** and **TDL-E** models represent propagation scenarios with **Line of Sight (LOS)** between transmitter and receiver. Unlike the NLOS models (TDL-A/B/C), under LOS the **first tap** follows a **Rician** distribution (mixture of a deterministic component + a diffuse Rayleigh component), while the remaining taps retain the usual Rayleigh behaviour.

### 2.1 Tap 0 Equation (Rice)

$$a_0(t) = \sqrt{P_0} \left( \sqrt{\frac{K}{K+1}} \cdot e^{j(2\pi f_S t + \phi_0)} + \sqrt{\frac{1}{K+1}} \cdot g_0(t) \right)$$

- **K**: linear K-factor of the profile (TDL-D: 13.3 dB = 21.38; TDL-E: 22.0 dB = 158.49)
- **f_S = 0.7 · f_D**: Doppler frequency of the LOS peak (as per TR 38.901)
- **φ_0**: initial phase uniformly distributed in [0, 2π), randomly generated at start-up
- **g_0(t)**: Rayleigh process (Jakes SoS, already present in the engine)
- **P_0 = 0 dB**: normalised total power of tap 0 (as per the 3GPP table note)

### 2.2 Taps 1..L (Rayleigh, same as NLOS models)

$$a_n(t) = \sqrt{P_n} \cdot g_n(t), \quad n \geq 1$$

---

## 3. Parameter Tables (DS = 50 ns)

### TDL-D — Tap 0 Rice + remaining Rayleigh taps (total 13 indices, tap[0] = Rice)

| Idx | Delay norm | Delay (ns) | Power (dB) | Type  |
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

**Note:** Tap 0 has K₁ = 13.3 dB and an average power of 0 dB.

### TDL-E — Tap 0 Rice + remaining Rayleigh taps (total 14 indices, tap[0] = Rice)

| Idx | Delay norm | Delay (ns) | Power (dB) | Type  |
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

**Note:** Tap 0 has K₁ = 22.0 dB and an average power of 0 dB.

---

## 4. Modified Files

Only **2 files** in the fading engine are modified. The rest of the stack requires no changes thanks to the existing architecture.

### 4.1 `lib/include/srsran/phy/channel/fading.h`

**Changes:**
1. Add `srsran_channel_fading_model_tdld` and `srsran_channel_fading_model_tdle` to the enum.
2. Add 3 fields to the `srsran_channel_fading_t` struct for the LOS state:
   - `bool is_los` — enables the Rice branch on tap 0
   - `float k_factor` — linear K for tap 0
   - `float los_phi0` — random initial phase φ₀

### 4.2 `lib/src/phy/channel/fading.c`

**Changes:**
1. Extend `nof_taps[9]` with 2 entries (TDL-D: 13, TDL-E: 14).
2. Extend `excess_tap_delay_ns[9][...]` with the TDL-D/E tables (DS=50 ns).
3. Extend `relative_power_db[9][...]` with the Rayleigh tap powers. To ensure a **total mean gain of 1.0 (0 dB)**, a normalisation offset has been applied to each table:
   - **TDL-D:** -0.1446 dB offset applied to all taps.
   - **TDL-E:** -0.1707 dB offset applied to all taps.
4. Add `parse_model()` entries for `tdld`/`tdl-d` and `tdle`/`tdl-e`.
5. Add helper `get_los_phasor(t, f_d, phi0) → cf_t`.
6. Modify `generate_taps()`: for `q->is_los` at tap `i==0`, combine scaled LOS + Rayleigh components.
7. Modify `srsran_channel_fading_init()`: detect LOS model, initialise `is_los`, `k_factor`, `los_phi0`.

---

## 5. Solution Architecture (Maximum Reuse)

fading_execute()
└─► generate_taps(q, time)
│
├─ tap i=0, is_los=true:
│ rayleigh_coeff = get_doppler_dispersion(...) ← REUSED
│ los_coeff = get_los_phasor(t, 0.7*fd, phi0) ← NEW
│ a = sqrt(P0) * (sqrt(K/(K+1))*los_coeff + sqrt(1/(K+1))*rayleigh_coeff)
│
└─ tap i>0 (or is_los=false for NLOS):
a = get_doppler_dispersion(...) ← UNCHANGED


The whole FFT/Overlap-Add pipeline in `filter_segment()` is reused without changes.

---

## 6. Compatibility

| Model | Status |
|--------|--------|
| `none` | ✅ Unchanged |
| `epaX`, `evaX`, `etuX` | ✅ Unchanged |
| `tdlaX`, `tdlbX`, `tdlcX` | ✅ Unchanged |
| `tdldX` | 🆕 TDL-D LOS (Rice tap 0, K=13.3 dB) |
| `tdleX` | 🆕 TDL-E LOS (Rice tap 0, K=22.0 dB) |

Usage in `ue.conf`:
```ini
[channel.ul.fading]
enable = true
model  = tdld3       ; TDL-D with 3 Hz Doppler
; model = tdle5      ; TDL-E with 5 Hz Doppler
```

---

## 7. Verification Strategy

- Clean build with `make srsue -j$(nproc)` from the `build/` directory.
- Unit test `fading_channel_test` with `tdld` and `tdle` models.
- NLOS regression with `tdla`, `epa` to confirm backwards compatibility.
- Integration in an srsUE session with LOS fading enabled.

See `validation.md` for the full results.
