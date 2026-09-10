from fdtd_core import FDTD2D, C0
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

sim = FDTD2D(300, 300, 0.00001)

# Viele Punktquellen mit gleicher Phase entlang einer Linie ueberlagern sich
# zu einer ebenen Welle, die sich in i-Richtung ausbreitet.
N_LAMBDA = 20                       # Gitterpunkte pro Wellenlaenge
wellenlaenge = N_LAMBDA * sim.d
f0 = C0 / wellenlaenge
i_quelle = sim.npml + 2
for j in range(1, sim.ny - 1):
    sim.add_soft_source(i_quelle, j, lambda t: np.sin(2 * np.pi * f0 * t))

X, Y = np.meshgrid(np.arange(sim.nx), np.arange(sim.ny), indexing="ij")

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.set_zlim(-1.0, 1.0)


def update(_frame):
    sim.step()
    ax.clear()
    ax.set_zlim(-1.0, 1.0)
    ax.set_xlabel("i")
    ax.set_ylabel("j")
    ax.set_title(f"t = {sim.t:.3e} s")
    print(sim.bx_i)
    # rstride/cstride: nur jeden 4. Gitterpunkt zeichnen, sonst ist ein
    # Redraw von 300x300 Punkten pro Frame zu langsam fuer eine Animation.
    return ax.plot_surface(X, Y, sim.Ez, cmap="RdBu", vmin=-1.0, vmax=1.0,
                            rstride=4, cstride=4)


ani = animation.FuncAnimation(fig, update, frames=600, interval=20)
plt.show()
