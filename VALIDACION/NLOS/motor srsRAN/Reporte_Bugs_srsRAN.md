# Reporte de Optimización y Corrección de Bugs en srsRAN (Motor de Fading TDL)

Este documento detalla el análisis, descubrimiento y resolución de dos *bugs* críticos encontrados matemáticamente en el motor nativo de canales con desvanecimiento (Fading) del emulador de `srsRAN_4G/5G`. 

El problema se identificó durante la validación estricta del perfil de retardo de potencia (Power Delay Profile - PDP) de las tramas TDL-A de la 3GPP TR 38.901. Al inyectar la señal perfecta a través del simulador físico `fading.c` para extraer los históricos en crudo, la energía conservada y la densidad espectral experimentaban severas oscilaciones (un factor de escala anómalo `~x17.7`).

## 1. Bug de Módulo vs Potencia Lineal
**Ubicación:** `srsRAN_4G/lib/src/phy/channel/fading.c` (Línea ~323)

### Análisis del Problema
Cuando el motor generaba los *Taps* del canal frecuencial (`generate_tap()`), los decibelios tabulados por la norma 3GPP (almacenados en `power_db`) se convertían a la base de cálculo de esta forma original:
```c
// [CÓDIGO ORIGINAL - ERROR]
float amplitude = srsran_convert_dB_to_power(power_db);
cf_t a0 = amplitude / N;
```

La función `srsran_convert_dB_to_power()` entrega la **potencia lineal** ($P_{lineal} = 10^{\text{dB}/10}$). Sin embargo, al modelar la onda matemáticamente en el dominio IQ en variables de tipo factor complejo `cf_t`, estamos trabajando en el dominio del **voltaje (o amplitud)** y no de la potencia. Al inyectar la potencia como si fuera la amplitud, la potencia real efectiva se convertía en el *cuadrado* de la potencia deseada ($(P_{lineal})^2$), destrozando completamente el balance y normalización unitaria (1.0) del PDP teorico de la 3GPP.

### Solución
Se ha aplicado la matemática física correcta al generador de la envolvente, extrayendo explícitamente la **raíz cuadrada** de la potencia lineal para obtener la amplitud analítica correcta:
```c
// [CÓDIGO CORREGIDO]
// Módulo H_tap[k]: La amplitud debe ser la RAÍZ CUADRADA de la potencia lineal
float amplitude = sqrtf(srsran_convert_dB_to_power(power_db));
float O         = (delay_ns * 1e-9f * srate + path_delay) / (float)N;
cf_t  a0        = amplitude / N;
```

---

## 2. Bug Estocástico del Método de Jakes (Desvanecimiento Doppler)
**Ubicación:** `srsRAN_4G/lib/src/phy/channel/fading.c` (Línea ~475)

### Análisis del Problema
Para emular un canal Rayleigh independiente sin línea de visión directa (NLOS), `srsRAN` implementa el principio "*Jakes' superposition of sinusoids*". Para que la envolvente decaiga siguiendo estrictamente una curva estadística Chi-Cuadrado de 2 grados de libertad, se deben sumar múltiples osciladores ($N_{terms}=16$) generando un *fading* constructivo-destructivo empleando retardos equiespaciados ($\alpha_j$).  

Al auditar la inicialización matemática, descubrimos el siguiente error de indexado maestro:
```c
// [CÓDIGO ORIGINAL - ERROR]
for (uint32_t i = 0; i < nof_taps[q->model]; i++) {
    for (uint32_t j = 0; (float)j < SRSRAN_CHANNEL_FADING_NTERMS; j++) {
        // ... inicializando variables `a` y `b` ...
        q->coeff_alpha[i][j] = ((float)M_PI * ((float)i - (float)0.5f)) / (2.0f * nof_taps[q->model]);
    }
}
```

¡Efectivamente! El desvío de fase `coeff_alpha` estaba empleando `i` (índice maestro del tap de retardo) en lugar de `j` (índice maestro del oscilador iterativo de Jakes). 
En este código anómalo, para un Tap concreto, la superposición de todos los $16$ senos electromagnéticos de simulación estaban **completamente sincronizados** en la misma fase temporal. El oscilador estocástico estaba roto y, sin promediarse geométricamente, forzaban estancamientos en niveles fijos lejanos al `1.0` empírico esperado para un desvanecimiento uncorrelated.

### Solución
La ecuación matemática de Jakes exige distribuir progresivamente los osciladores internos basándose en el total de osciladores, no en base a cuántos "taps" compongan el modelo general TDL-A. 
Se procedió a corregir la fórmula empleando el divisor estático y el índice iterativo adecuado de sub-estado:

```c
// [CÓDIGO CORREGIDO]
for (uint32_t j = 0; (float)j < SRSRAN_CHANNEL_FADING_NTERMS; j++) {
    q->coeff_a[i][j]     = srsran_random_uniform_real_dist(random, 0, 2.0f * (float)M_PI);
    q->coeff_b[i][j]     = srsran_random_uniform_real_dist(random, 0, 2.0f * (float)M_PI);
    q->coeff_alpha[i][j] = ((float)M_PI * ((float)j - (float)0.5f)) / (2.0f * SRSRAN_CHANNEL_FADING_NTERMS);
}
```

## Conclusión

Tras parchear ambos módulos base en `C`, **recompilamos** la librería central `Libsrsran_phy`. Seguidamente, se ejecutó una validación cruda C-to-Python mediante `srsran_fading_sim.c` forzando $10,000$ iteraciones físicas. 

Al observar la conservación de la energía y comparar la curva de densidad espectral:
**Sin ningún tipo de factor corrector, la dispersión Rayleigh converge matemáticamente a 1.0 y el Histograma de la muestra $n_{peak}$ traza de forma impecable la Función de Densidad de la Probabilidad Estándar.**
