#!/usr/bin/env python3
"""
Genera las figuras didacticas para docs/heat-equation.md.

Reproduce EXACTAMENTE el solver de heat-big.c / heat-mpi-big.c:
  - Dominio: cuadrado unitario [0,1] x [0,1], malla 80 x 80.
  - Esquema: Euler explicito (Jacobi) con stencil de 5 puntos.
  - dt = min(dx^2, dy^2) / 4   (limite de estabilidad).
  - Condiciones de frontera (Dirichlet) identicas al codigo C:
        interior            phi = 0
        columna derecha     phi = 1                (borde caliente)
        fila superior/inf.  phi = k*dx             (rampa lineal 0 -> 1)
        columna izquierda   phi = 0
  - Convergencia: max|dphi| < eps  (estado estacionario => ecuacion de Laplace).

Salida: docs/img/*.png y docs/img/heat-convergence.gif

Uso:  python3 docs/generate_figures.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation

HERE = os.path.dirname(os.path.abspath(__file__))
IMG = os.path.join(HERE, "img")
os.makedirs(IMG, exist_ok=True)

# ----- parametros, identicos a heat-big.c -----
imax = kmax = 80           # imax = direccion y (filas), kmax = direccion x (columnas)
eps = 1e-8
itmax = 20000
dx = 1.0 / kmax
dy = 1.0 / imax
dx2i = 1.0 / dx**2
dy2i = 1.0 / dy**2
dt = min(dx**2, dy**2) / 4.0


def init_field():
    """phi[i, k] : i = y (fila), k = x (columna). Mismas BC que el codigo C."""
    phi = np.zeros((imax + 1, kmax + 1))
    k = np.arange(kmax + 1)
    phi[:, kmax] = 1.0          # columna derecha (x = 1): caliente
    phi[0, :] = k * dx          # fila superior: rampa 0 -> 1
    phi[imax, :] = k * dx       # fila inferior: rampa 0 -> 1
    phi[:, kmax] = 1.0          # re-fija esquinas derechas a 1 (como el codigo)
    return phi


def step(phi):
    """Un paso Jacobi explicito. Devuelve (phi_nuevo, max|dphi|)."""
    lap = (
        (phi[2:, 1:-1] + phi[:-2, 1:-1] - 2 * phi[1:-1, 1:-1]) * dy2i
        + (phi[1:-1, 2:] + phi[1:-1, :-2] - 2 * phi[1:-1, 1:-1]) * dx2i
    )
    dphi = lap * dt
    phi[1:-1, 1:-1] += dphi
    return phi, np.max(dphi)


def solve(capture_at=None):
    """Itera hasta convergencia. Captura instantaneas en las iteraciones dadas."""
    phi = init_field()
    snaps, hist = {}, []
    capture = set(capture_at or [])
    if 0 in capture:
        snaps[0] = phi.copy()
    it = 0
    for it in range(1, itmax + 1):
        phi, dmax = step(phi)
        hist.append(dmax)
        if it in capture:
            snaps[it] = phi.copy()
        if dmax < eps:
            break
    snaps[it] = phi.copy()   # estado final
    return phi, snaps, hist, it


def save_heatmap(phi, title, fname):
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    # extent: x (k) de 0..1, y (i) de 0..1; origin lower para que y crezca hacia arriba
    im = ax.imshow(phi, origin="lower", extent=[0, 1, 0, 1],
                   cmap="inferno", aspect="equal", vmin=0, vmax=1)
    ax.set_xlabel("x  (indice k)")
    ax.set_ylabel("y  (indice i)")
    ax.set_title(title)
    cb = fig.colorbar(im, ax=ax, shrink=0.85)
    cb.set_label("temperatura  phi")
    fig.tight_layout()
    fig.savefig(os.path.join(IMG, fname), dpi=110)
    plt.close(fig)
    print("  ->", fname)


def main():
    print("Resolviendo y generando figuras...")
    frames_at = [0, 5, 20, 60, 150, 400, 1000, 3000, 8000]
    phi, snaps, hist, it_conv = solve(capture_at=frames_at)
    print(f"  convergio en {it_conv} iteraciones (eps={eps})")

    # 1) condicion inicial
    save_heatmap(snaps[0], "Condicion inicial (t = 0)", "heat-initial.png")

    # 2) estado estacionario
    save_heatmap(phi, f"Estado estacionario (convergio en {it_conv} iter.)",
                 "heat-steady.png")

    # 3) curva de convergencia (max|dphi| vs iteracion, escala log)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.semilogy(range(1, len(hist) + 1), hist, color="#c0392b")
    ax.axhline(eps, ls="--", color="gray", label=f"eps = {eps:g}")
    ax.set_xlabel("iteracion")
    ax.set_ylabel("max |dphi|  (cambio maximo por paso)")
    ax.set_title("Convergencia hacia el estado estacionario")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(IMG, "heat-convergence-curve.png"), dpi=110)
    plt.close(fig)
    print("  -> heat-convergence-curve.png")

    # 4) mosaico de instantaneas
    order = [k for k in frames_at if k in snaps] + (
        [it_conv] if it_conv not in frames_at else [])
    fig, axes = plt.subplots(2, 5, figsize=(15, 6.2))
    for ax, key in zip(axes.flat, order):
        ax.imshow(snaps[key], origin="lower", extent=[0, 1, 0, 1],
                  cmap="inferno", aspect="equal", vmin=0, vmax=1)
        lbl = "inicial" if key == 0 else f"iter {key}"
        if key == it_conv:
            lbl += " (final)"
        ax.set_title(lbl, fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
    for ax in axes.flat[len(order):]:
        ax.axis("off")
    fig.suptitle("Difusion del calor: de la condicion inicial al equilibrio",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(IMG, "heat-stages.png"), dpi=100)
    plt.close(fig)
    print("  -> heat-stages.png")

    # 5) animacion GIF
    print("  generando animacion (puede tardar)...")
    phi_a = init_field()
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    im = ax.imshow(phi_a, origin="lower", extent=[0, 1, 0, 1],
                   cmap="inferno", aspect="equal", vmin=0, vmax=1)
    ax.set_xlabel("x (k)"); ax.set_ylabel("y (i)")
    fig.colorbar(im, ax=ax, shrink=0.85, label="phi")
    title = ax.set_title("iteracion 0")

    state = {"phi": phi_a}
    steps_per_frame = 25
    n_frames = 90

    def update(frame):
        for _ in range(steps_per_frame):
            state["phi"], _ = step(state["phi"])
        im.set_array(state["phi"])
        title.set_text(f"iteracion {frame * steps_per_frame}")
        return im, title

    anim = animation.FuncAnimation(fig, update, frames=n_frames,
                                   interval=80, blit=False)
    gif = os.path.join(IMG, "heat-convergence.gif")
    anim.save(gif, writer=animation.PillowWriter(fps=12))
    plt.close(fig)
    print("  -> heat-convergence.gif")
    print("Listo.")


if __name__ == "__main__":
    main()
