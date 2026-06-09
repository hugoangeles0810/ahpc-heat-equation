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

## Estructura del repositorio

```
├── src/                 # solvers en C
│   ├── heat-big.c       # versión serial de referencia
│   └── heat-mpi-big.c   # versión paralela con MPI
├── scripts/             # cómo ejecutar / medir el solver
│   ├── run_local.sh     # smoke test rápido con mpirun (laptop)
│   ├── benchmark.slurm   # barrido p × n en el clúster (SLURM)
│   └── collect_results.sh # parsea la salida del solver → CSV
├── analysis/            # objetivo (c): métricas y gráficas
│   └── plot_performance.py # CSV → tiempo / speedup / eficiencia
├── results/             # CSVs y figuras finales (versionados)
├── docs/                # guía didáctica + figuras
└── Makefile             # make → bin/heat-big, bin/heat-mpi-big
```

## Uso

```bash
make                      # compila ambos binarios en bin/
scripts/run_local.sh 4    # corre el solver MPI con 4 procesos (local)
sbatch scripts/benchmark.slurm           # barrido de escalabilidad en el clúster
python analysis/plot_performance.py      # genera las gráficas desde results/benchmark.csv
```

## Referencias

Material de apoyo a considerar para el objetivo (c) — definición de métricas y
metodología de medición de performance:

- **[Serial vs. Parallel Programming with MPI — NMSU HPC](https://hpc.nmsu.edu/discovery/mpi/serial-vs-parallel/)**
  Caso práctico serial vs. MPI: muestra cómo, para problemas pequeños, el overhead
  de inicialización/comunicación de MPI hace que la versión serial sea más rápida, y
  cómo la versión paralela gana terreno al crecer `n`. Útil para interpretar por qué
  el speedup y la eficiencia dependen del tamaño del problema en nuestro barrido `p × n`.
- **[Measuring performance — LAMMPS docs](https://docs.lammps.org/Speed_measure.html)**
  Distingue performance serial (velocidad de un proceso) de eficiencia paralela, y
  define la **eficiencia** como `rendimiento(N procesos) / (rendimiento(1 proceso) × N)`.
  Discute *strong* vs *weak scaling* y la degradación de la eficiencia al aumentar
  procesos cuando el trabajo por proceso es pequeño frente a la comunicación — la base
  metodológica de las métricas que calcula `analysis/plot_performance.py`.

## Equipo

- _Hugo Angeles_
- _Christian Cajusol_
- _Jhomar Yurivilca_
- _Francisco Meza_
