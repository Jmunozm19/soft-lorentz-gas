# Soft Lorentz Gas – Diffusion Simulations  
### Proyecto de Mecánica Estadística Avanzada (PUC)

Este repositorio contiene el código, figuras y notebooks utilizados para reproducir y estudiar el comportamiento difusivo en el **Soft Lorentz Gas** con potencial suavizado tipo Fermi.

## Contenido del repositorio

/src  
    simulate_soft_lorentz_diffusion_fixed.py  
/notebooks  
    Resultados_finales.ipynb  
/figures  
    ...  

## Modelo físico

Potencial Fermi:
V(r) = 1 / (1 + exp((r - r0)/sigma))

## Ejecución

from simulate_soft_lorentz_diffusion_fixed import simulate_msd
times, msd = simulate_msd(3.0)

## Autores

Carlos Alberto Meza Morales  
Jesús David Muñoz Muñoz  
