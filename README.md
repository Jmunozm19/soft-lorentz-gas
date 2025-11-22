# Soft Lorentz Gas – Diffusion Simulations  
### Proyecto de Mecánica Estadística Avanzada (PUC)

Este repositorio contiene el código, figuras y notebooks utilizados para estudiar el comportamiento difusivo en un **Soft Lorentz Gas** bidimensional con potencial suavizado tipo Fermi.  
El núcleo del proyecto es el archivo:

- `src/simulate_soft_lorentz_diffusion_fixed.py`

donde se implementa la dinámica de partículas, el cálculo del **Mean Squared Displacement (MSD)** y la estimación del **coeficiente de difusión** \(D(E)\).

---

# 📦 Método Numérico y Estructura del Código

Esta sección describe la estructura del código utilizado para simular la dinámica del **Soft Lorentz Gas** y calcular el **coeficiente de difusión** \( D(E) \).  
Todo el código del proyecto está disponible en este repositorio.

---

## 🔧 1. Estructura general del código

El proyecto está organizado en tres componentes principales:

### **1. Módulo de dinámica**  
📄 `src/simulate_soft_lorentz_diffusion_fixed.py`  
Contiene toda la implementación del modelo físico y del algoritmo numérico:

- definición del potencial suave tipo Fermi y su derivada,  
- cálculo del potencial total y las fuerzas,  
- implementación de condiciones periódicas reales,  
- integrador Velocity Verlet,  
- cálculo del MSD,  
- estimación del coeficiente de difusión \( D(E) \).

---

### **2. Notebook de análisis**  
📄 `notebooks/Resultados_finales.ipynb`  

Funciona como *interfaz de usuario* del proyecto:

- realiza barridos en energía,  
- ejecuta la simulación usando el módulo de dinámica,  
- genera las gráficas (MSD vs t, mesetas, D(E), ajustes log–log),  
- guarda resultados y figuras.

---

## 🧠 2. Módulo de dinámica: funciones esenciales

### 🔸 `V_fermi(r)` y `dVdr_fermi(V)`
Implementan el potencial suave tipo Fermi y su derivada radial.

### 🔸 `potential_at(r)`
Suma la contribución del potencial Fermi para todos los centros de la red.

### 🔸 `force_at(r)`
Calcula la fuerza total:

\[
\mathbf{F}(\mathbf r) = -\nabla V_{\text{tot}}(\mathbf r)
\]

### 🔸 `wrap_to_cell(r)`
Implementa **condiciones periódicas exactas**, devolviendo:
- la posición envuelta en la celda unitaria,  
- el desplazamiento de red acumulado.

### 🔸 `sample_positions(N, E)`
Genera posiciones iniciales uniformes garantizando energía cinética positiva.

### 🔸 `sample_velocities(N, E, V_at_r)`
Asigna velocidades iniciales compatibles con energía total fija \(E\).

### 🔸 `simulate_msd(E, ...)`
Ejecuta:
1. muestreo de posiciones y velocidades,  
2. integración temporal con Verlet,  
3. condiciones periódicas,  
4. reconstrucción absoluta,  
5. cálculo del **MSD**.

Retorna: tiempos y valores de MSD.

### 🔸 `estimate_D_from_msd(times, msd)`
Ajusta:

\[
\mathrm{MSD}(t) \approx 4Dt
\]

y obtiene:
- \(D\),  
- error de ajuste,  
- pendiente,  
- intercepto.

---

## 📓3. Notebooks y figuras

- `notebooks/Resultados_finales.ipynb`:
  - define energías,  
  - ejecuta simulaciones,  
  - extrae \( D(E) \),  
  - genera curvas y figuras finales.
## 📊 Resultados principales

### 1. Meseta difusiva (ambos regímenes)
<p align="center">
  <img src="figures/MSD_over_4t_both.png" width="600">
  <br>
  <em>MSD/(4t) mostrando regiones cercanas al umbral y altas energías.</em>
</p>

### 2. Región cercana al umbral \(E \to 1^{+}\)
<p align="center">
  <img src="figures/DE_near1_fit.png" width="600">
  <br>
  <em>Ajuste en escala log–log: \(D \propto (E - 1)^{b}\).</em>
</p>

### 3. Alta energía – dependencia \(D(E)\)
<p align="center">
  <img src="figures/DE_loglog_fit_highE.png" width="600">
  <br>
  <em>Ajuste lineal en log–log para el régimen de energías grandes.</em>
</p>

---

## 📈 MSD detallado en cada régimen

### MSD vs t (región \(E \to 1^{+}\))
<p align="center">
  <img src="figures/MSD_vs_t_near1.png" width="600">
  <br>
  <em>Comportamiento lineal del MSD en la cola temporal para energías cercanas al umbral.</em>
</p>

### MSD vs t (alta energía)
<p align="center">
  <img src="figures/MSD_vs_t_alta_energia.png" width="600">
  <br>
  <em>Comportamiento lineal bien definido en la región difusiva para energías grandes.</em>
</p>

---

## 🧭 Geometría del sistema

### Campo de potencial en la celda unitaria
<p align="center">
  <img src="figures/soft_lorentz_unitcell.png" width="600">
  <br>
  <em>Mapa del potencial suavizado tipo Fermi en la red triangular.</em>
</p>

### Perfil 1D del potencial en la base de la celda
<p align="center">
  <img src="figures/soft_lorentz_profile.png" width="600">
  <br>
  <em>Visualización del mínimo del canal y del ancho efectivo accesible para distintas energías.</em>
</p>

---

## 📝 Resumen

Este repositorio reúne las herramientas necesarias para simular y analizar el comportamiento difusivo en el **Soft Lorentz Gas** mediante dinámica clásica en un potencial periódico suavizado.

Para utilizar el código y recrear los resultados del proyecto, se recomienda descargar y ejecutar los siguientes archivos de Python del repositorio:

- `src/simulate_soft_lorentz_diffusion_fixed.py`  
- `notebooks/Resultados_finales.ipynb` 

El archivo `simulate_soft_lorentz_diffusion_fixed.py` contiene las funciones esenciales del modelo:

- Definición de la red triangular y sus vectores de base  
- Implementación del potencial suave tipo Fermi y su fuerza asociada  
- Manejo de condiciones periódicas reales mediante “wrapping” y acumulación de desplazamientos de red  
- Integración temporal con el esquema de Velocity Verlet  
- Cálculo del **Mean Squared Displacement (MSD)** y estimación del coeficiente de difusión efectivo \(D(E)\)

El notebook `Resultados_finales.ipynb` muestra cómo:

- Ejecutar las simulaciones para distintos valores de energía  
- Graficar MSD vs. tiempo  
- Ajustar el régimen difusivo para extraer \(D\)  
- Construir la curva \(D(E)\) y comparar tendencias con la literatura

---

### 🔬 Extensibilidad del código

El diseño modular del programa permite usar este repositorio como **base para estudiar otros sistemas**.  
En particular, es posible:

- Cambiar el potencial (por ejemplo, otros perfiles suaves o duros)  
- Modificar la geometría de la red (cuadrada, hexagonal, etc.) alterando los vectores `a1` y `a2`  
- Incorporar nuevos términos en la dinámica (campos externos, masas distintas, etc.)

De este modo, quien descargue los archivos de Python puede reutilizar la estructura general del código para explorar **modelos de difusión y transporte en otros medios periódicos**, sin tener que reescribir desde cero el framework numérico.

---

## 👥 Autores

- **Carlos Alberto Meza Morales** – Pontificia Universidad Católica de Chile  
- **Jesús David Muñoz Muñoz** – Pontificia Universidad Católica de Chile  

Proyecto desarrollado en el contexto del curso **Mecánica Estadística Avanzada (FIM8451-1)**.

---

## 📜 Licencia

El código se distribuye bajo licencia MIT (o la que se defina en el archivo `LICENSE`).  
Puedes usarlo, modificarlo y distribuirlo citando adecuadamente este repositorio.
