# Soft Lorentz Gas – Diffusion Simulations  
### Proyecto de Mecánica Estadística Avanzada (PUC)

Este repositorio contiene el código, figuras y notebooks utilizados para estudiar el comportamiento difusivo en un **Soft Lorentz Gas** bidimensional con potencial suavizado tipo Fermi.  
El núcleo del proyecto es el archivo:

- `src/simulate_soft_lorentz_diffusion_fixed.py`

donde se implementa la dinámica de partículas, el cálculo del **Mean Squared Displacement (MSD)** y la estimación del **coeficiente de difusión** \(D(E)\).

---

## 🧠 Idea general del código

El código simula un conjunto de partículas puntuales que se mueven en un potencial periódico generado por una **red triangular** de “discos suaves”.  
Para una energía total \(E\), las partículas obedecen la dinámica clásica:

- \\( \dot{\mathbf r} = \mathbf v \\)
- \\( \dot{\mathbf v} = \mathbf F(\mathbf r) = -\nabla V(\mathbf r) \\)

donde el potencial total es una suma de potenciales de tipo Fermi centrados en cada punto de la red.

A partir de las trayectorias se calcula:

\[
\text{MSD}(t) = \big\langle \lvert \mathbf r(t) - \mathbf r(0) \rvert^2 \big\rangle,
\]

y en régimen difusivo se ajusta la ley:

\[
\text{MSD}(t) \approx 4 D\, t,
\]

para obtener el coeficiente de difusión efectivo \(D\).

---

## ⚙️ Parámetros globales y geometría de la red

Al inicio del archivo se definen los parámetros geométricos y de potencial:

```python
r0 = 1.0
sigma = 0.01
w = 0.05
L = 2.0*r0 + w

a1 = np.array([L, 0.0])
a2 = np.array([L/2.0, np.sqrt(3.0)*L/2.0])
A = np.column_stack((a1, a2))
Ainv = np.linalg.inv(A)
```

- `r0`: radio efectivo del scatterer suave.  
- `sigma`: suavidad del borde del potencial Fermi (controla qué tan “blando” es el disco).  
- `w`: *gap* (distancia mínima entre discos).  
- `L = 2*r0 + w`: periodo básico de la red.  
- `a1`, `a2`: vectores base de una **red triangular**.  
- `A`: matriz cuyas columnas son `a1` y `a2`.  
- `Ainv`: inversa de `A`, usada para pasar de coordenadas cartesianas a coordenadas de red \\( (u,v) \\).

También se construye un conjunto finito de puntos de la red cercana:

```python
def lattice_points(nrange=2):
    pts = []
    for i in range(-nrange, nrange+1):
        for j in range(-nrange, nrange+1):
            pts.append(i*a1 + j*a2)
    return np.array(pts, dtype=float)

LP = lattice_points(2)
```

- `lattice_points(nrange)`: genera todos los puntos de red \\( i\,\mathbf a_1 + j\,\mathbf a_2 \\) con \\( i,j \in [-nrange, nrange] \\).
- `LP`: arreglo de centros de potencial que se usan para sumar las contribuciones al potencial total y a la fuerza.

---

## 📉 Potencial de tipo Fermi y su derivada

El **potencial suave** asociado a un disco se define como:

```python
def V_fermi(r, r0=r0, sigma=sigma):
    return 1.0 / (1.0 + np.exp((r - r0)/sigma))
```

Matemáticamente:

\[
V(r) = \frac{1}{1 + e^{(r - r_0)/\sigma}}.
\]

- Para \\( r \ll r_0 \\): el potencial es cercano a 1.  
- Para \\( r \gg r_0 \\): decae suavemente hacia 0.  
- `sigma` controla la pendiente en la región de borde.

La derivada radial efectiva (usada para la fuerza) se implementa como:

```python
def dVdr_fermi(V, sigma=sigma):
    return -(1.0/sigma) * V * (1.0 - V)
```

Usando la identidad:

\[
\frac{dV}{dr} = -\frac{1}{\sigma} V(1 - V).
\]

---

## 🧲 Potencial total y fuerza

### Potencial total en un punto

```python
def potential_at(r):
    r = np.atleast_2d(r)
    Vtot = np.zeros(r.shape[0])
    for R in LP:
        d = r - R
        dist = np.linalg.norm(d, axis=1)
        Vtot += V_fermi(dist)
    return Vtot
```

- Entrada: `r` puede ser un solo punto 2D o un arreglo `N x 2` de posiciones.
- Para cada centro de la red `R` en `LP`:
  - Se calcula el vector \\( \mathbf d = \mathbf r - \mathbf R \\) y su norma.
  - Se suma la contribución Fermi \\( V_\text{fermi}(|\mathbf d|) \\).
- Salida: un vector de potenciales totales \\( V_\text{tot}(\mathbf r_i) \\).

### Fuerza total

```python
def force_at(r):
    r = np.atleast_2d(r)
    F = np.zeros_like(r)
    for R in LP:
        d = r - R
        dist = np.linalg.norm(d, axis=1) + 1e-14
        V = V_fermi(dist)
        dVdr = dVdr_fermi(V)
        coeff = -dVdr / dist
        F += (coeff[:,None] * d)
    return F
```

- Calcula la fuerza \\( \mathbf F = -\nabla V \\) usando la derivada radial:  
  \\( \mathbf F = -\frac{dV}{dr} \frac{\mathbf d}{r} \\).  
- `1e-14` evita divisiones por cero cuando `dist` es muy pequeño.  
- La fuerza total es la suma sobre todos los centros en `LP`.

---

## 🔁 Condiciones periódicas: `wrap_to_cell`

Para imponer periodicidad y trabajar en una celda fundamental, se define:

```python
def wrap_to_cell(r):
    uv = (Ainv @ r.T).T
    n = np.floor(uv)
    uv_wrapped = uv - n
    r_wrapped = (A @ uv_wrapped.T).T
    shifts = n[...,0][:,None]*a1 + n[...,1][:,None]*a2
    return r_wrapped, shifts
```

Pasos:

1. Se mapean las posiciones cartesianas `r` a coordenadas de red `uv` usando `Ainv`.  
2. Se separan las partes enteras `n = floor(uv)` (cuántas celdas se ha salido la partícula).  
3. `uv_wrapped = uv - n` restringe las coordenadas al intervalo `[0,1)` (celda fundamental).  
4. Se reconstruyen las posiciones envueltas: `r_wrapped = A @ uv_wrapped`.  
5. Se calculan los **desplazamientos de red** efectivos (`shifts`), que luego se acumulan para reconstruir la trayectoria absoluta.

Este mecanismo permite que la dinámica se ejecute siempre dentro de la celda fundamental, pero al mismo tiempo conservar la información de cuántas celdas se han cruzado para medir desplazamientos reales.

---

## 🎯 Muestreo de condiciones iniciales

### Posiciones iniciales: `sample_positions`

```python
def sample_positions(N, E, max_tries=10000, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    u = rng.random((N,)); v = rng.random((N,))
    r = u[:,None]*a1 + v[:,None]*a2
    V = potential_at(r)
    tries = 0
    mask = (E - V) <= 0.0
    while np.any(mask):
        M = np.sum(mask)
        u = rng.random((M,)); v = rng.random((M,))
        r_new = u[:,None]*a1 + v[:,None]*a2
        V_new = potential_at(r_new)
        r[mask] = r_new; V[mask] = V_new
        mask = (E - V) <= 0.0
        tries += 1
        if tries > max_tries:
            raise RuntimeError("Could not place enough initial positions for E={}".format(E))
    return r, V
```

- Se eligen posiciones `r` al azar en la celda fundamental combinando `a1` y `a2` con coeficientes uniformes `u, v ∈ [0,1)`.  
- Se evalúa el potencial `V(r)`.  
- Se **rechazan** las posiciones que violan \\( E - V(r) > 0 \\) (es decir, donde la energía cinética sería negativa).  
- Se repite hasta que todas las partículas tengan \\( K = E - V(r) > 0 \\) o se alcance `max_tries`.

Resultado:

- `r`: posiciones iniciales aceptadas.  
- `V`: potencial en esas posiciones.

### Velocidades iniciales: `sample_velocities`

```python
def sample_velocities(N, E, V_at_r, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    K = E - V_at_r
    vmag = np.sqrt(2.0 * K)
    phi = 2.0*np.pi*rng.random((N,))
    return np.column_stack((vmag*np.cos(phi), vmag*np.sin(phi)))
```

- Se calcula la energía cinética local: \\( K_i = E - V(\mathbf r_i) \\).  
- Magnitud de la velocidad: \\( v_i = \sqrt{2 K_i} \\) (masa tomada igual a 1).  
- La dirección del vector velocidad se elige uniforme en \\( [0,2\pi) \\).  
- Se devuelve un arreglo `N x 2` de velocidades iniciales.

---

## ⏱️ Dinámica y cálculo del MSD: `simulate_msd`

```python
def simulate_msd(E, N=200, dt=1e-3, steps=20000, sample_stride=20, rng=None):
    if rng is None:
        rng = np.random.default_rng(2025)
    r, V = sample_positions(N, E, rng=rng)
    v = sample_velocities(N, E, V, rng=rng)
    F = force_at(r)

    # FIX: keep cumulative lattice shift and absolute positions
    L_cum = np.zeros_like(r)
    r_abs = r.copy()
    r_abs0 = r_abs.copy()

    times = [0.0]; msd = [0.0]
    for k in range(1, steps+1):
        v_half = v + 0.5*dt*F
        r_new = r + dt*v_half
        r_wrapped, shifts = wrap_to_cell(r_new)
        L_cum += shifts
        r = r_wrapped
        F_new = force_at(r)
        v = v_half + 0.5*dt*F_new
        F = F_new

        if (k % sample_stride)==0:
            r_abs = r + L_cum
            t = k*dt
            disp = r_abs - r_abs0
            msd_t = np.mean(np.sum(disp*disp, axis=1))
            times.append(t); msd.append(msd_t)
    return np.array(times), np.array(msd)
```

### Entradas principales

- `E`: energía total fijada para todas las partículas.  
- `N`: número de partículas que se simulan en paralelo.  
- `dt`: paso de tiempo del integrador.  
- `steps`: número total de pasos de integración.  
- `sample_stride`: cada cuántos pasos se registra el MSD.  
- `rng`: generador de números aleatorios (permite reproducibilidad).

### Flujo del algoritmo

1. **Inicialización**  
   - `sample_positions`: se generan posiciones iniciales con \\( K>0 \\).  
   - `sample_velocities`: se asignan velocidades compatibles con la energía total.  
   - Se evalúa la fuerza inicial `F = force_at(r)`.

2. **Acumulador de desplazamientos de red**  
   - `L_cum`: vector que acumula los desplazamientos de red (cuántas celdas se han cruzado en cada dirección).  
   - `r_abs0`: posiciones absolutas iniciales (sirven como referencia para el MSD).

3. **Integración temporal – esquema Velocity Verlet**

   Para cada paso `k`:

   - Actualización de velocidad a mitad de paso:
     \\( \mathbf v_{1/2} = \mathbf v + \frac{dt}{2} \mathbf F \\).
   - Predicción de nueva posición:
     \\( \mathbf r_\text{new} = \mathbf r + dt\,\mathbf v_{1/2} \\).
   - Se aplica `wrap_to_cell`:
     - `r_wrapped`: posición dentro de la celda fundamental.  
     - `shifts`: desplazamiento de red asociado a ese cruce de frontera.  
   - Se acumula `L_cum += shifts`.  
   - Se recalcula la fuerza con `F_new = force_at(r_wrapped)`.  
   - Se corrige la velocidad al final del paso:
     \\( \mathbf v_\text{new} = \mathbf v_{1/2} + \frac{dt}{2} \mathbf F_\text{new} \\).

4. **Cálculo del MSD**

   Cada `sample_stride` pasos:

   - Se reconstruyen posiciones absolutas:  
     \\( \mathbf r_\text{abs} = \mathbf r_\text{wrapped} + L_\text{cum} \\).  
   - Se calcula el desplazamiento respecto a las posiciones iniciales:  
     \\( \Delta \mathbf r_i = \mathbf r_{\text{abs},i}(t) - \mathbf r_{\text{abs},i}(0) \\).  
   - Se calcula:
     \\( \text{MSD}(t) = \frac{1}{N}\sum_i |\Delta \mathbf r_i|^2 \\).  
   - Se almacenan los valores en `times` y `msd`.

5. **Salida**

   - `times`: arreglo con los tiempos muestreados.  
   - `msd`: arreglo con los valores correspondientes de MSD.

---

## 📏 Ajuste lineal y cálculo de D: `estimate_D_from_msd`

```python
def estimate_D_from_msd(times, msd, frac_window=0.5):
    n = len(times)
    start = int(n*(1.0 - frac_window))
    if start < 1: start = 1
    t_fit = times[start:]; m_fit = msd[start:]
    A = np.vstack([t_fit, np.ones_like(t_fit)]).T
    slope, intercept = np.linalg.lstsq(A, m_fit, rcond=None)[0]
    D = slope/4.0
    resid = m_fit - (slope*t_fit + intercept)
    s_err = np.sqrt(np.mean(resid**2))
    denom = np.sum((t_fit - np.mean(t_fit))**2) + 1e-12
    slope_err = s_err / np.sqrt(denom)
    D_err = slope_err / 4.0
    return D, D_err, slope, intercept
```

### ¿Qué hace esta función?

1. **Ventana de ajuste**  
   - Usa sólo una fracción final de los datos (`frac_window`, por defecto 50%) para ajustar la recta en el régimen supuestamente difusivo.

2. **Ajuste por mínimos cuadrados**  
   - Ajusta \\( \text{MSD}(t) \approx \text{slope} \cdot t + \text{intercept} \\).  
   - De la pendiente se extrae:
     \\( D = \text{slope} / 4 \\).

3. **Estimación de error**  
   - Calcula los residuos del ajuste y a partir de ellos una estimación del error cuadrático medio.  
   - Propaga esa incertidumbre a la pendiente y luego a `D_err`.

### Salida

- `D`: estimación del coeficiente de difusión.  
- `D_err`: estimación de la incertidumbre en `D`.  
- `slope`, `intercept`: parámetros de la recta ajustada.

---

## 💻 Ejemplos de uso

### 1. Simular el MSD para una energía fija

```python
from simulate_soft_lorentz_diffusion_fixed import simulate_msd, estimate_D_from_msd

E = 3.0
times, msd = simulate_msd(E, N=200, dt=1e-3, steps=20000, sample_stride=20)

D, D_err, slope, intercept = estimate_D_from_msd(times, msd)

print(f"E = {E}")
print(f"D = {D:.4f} ± {D_err:.4f}")
```

### 2. Barrido en energía y construcción de D(E)

```python
import numpy as np
from simulate_soft_lorentz_diffusion_fixed import simulate_msd, estimate_D_from_msd

energies = np.linspace(1.0, 7.0, 13)
D_vals = []
D_errs = []

for E in energies:
    times, msd = simulate_msd(E, N=200, dt=1e-3, steps=20000, sample_stride=20)
    D, D_err, _, _ = estimate_D_from_msd(times, msd)
    D_vals.append(D)
    D_errs.append(D_err)
    print(f"E = {E:.2f} -> D = {D:.4f} ± {D_err:.4f}")
```

Con estos datos se pueden generar figuras de **D vs E** y comparar con los resultados reportados en la literatura.

---

## 📓 Notebooks y figuras

- `notebooks/Resultados_finales.ipynb`:
  - Contiene ejemplos de corridas, gráficas de MSD(t), ajustes lineales y la curva D(E).  
  - Sirve como cuaderno de trabajo donde se documentan los parámetros utilizados y se guardan las figuras finales.

- `figures/`:
  - Se recomienda guardar aquí:
    - MSD vs. t para diferentes energías.
    - Ajustes lineales de la parte difusiva.
    - Gráfica D(E) (con barras de error).
    - Mapas de potencial 2D si se generan.

---

## 👥 Autores

- **Carlos Alberto Meza Morales** – Pontificia Universidad Católica de Chile  
- **Jesús David Muñoz Muñoz** – Pontificia Universidad Católica de Chile  

Proyecto desarrollado en el contexto del curso **Mecánica Estadística Avanzada (PHYS-4035)**.

---

## 📜 Licencia

El código se distribuye bajo licencia MIT (o la que se defina en el archivo `LICENSE`).  
Puedes usarlo, modificarlo y distribuirlo citando adecuadamente este repositorio.
