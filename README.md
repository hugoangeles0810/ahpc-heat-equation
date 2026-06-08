# Heat Equation — Proyecto AHPC

Análisis de performance de un código **MPI** que resuelve la **ecuación de calor**
(una ecuación diferencial parcial, EDP) en paralelo. Proyecto del curso *Applied High
Performance Computing*, Posgrado UTEC 2026-I.

## Contexto

La ecuación de calor modela la transferencia/difusión de temperatura en un dominio a lo
largo del tiempo. Se resuelve numéricamente discretizando el dominio en una malla y
actualizando cada punto en función de sus vecinos. El trabajo consiste en evaluar cómo
escala esta solución paralela al variar el número de procesos y el tamaño del problema.

> 📖 Para una explicación didáctica de la física, la discretización y
> cómo se mapea cada concepto a la implementación en C y MPI ver [`docs/heat-equation.md`](docs/heat-equation.md).

## Objetivos

- **a) PRAM:** elaborar el modelo PRAM del algoritmo y determinar su complejidad teórica.
- **b) Mediciones:** medir tiempos de ejecución y compararlos con la complejidad teórica,
  para distinto número de procesos (`p`) y tamaño del problema (`n`).
- **c) Análisis de performance:** desarrollar software que calcule y grafique las métricas
  (tiempo de ejecución, **speedup**, **eficiencia**) y concluir sobre la escalabilidad.

## Equipo

- _Hugo Angeles_
- _Christian Cajusol_
- _Jhomar Yurivilca_
- _Francisco Meza_
