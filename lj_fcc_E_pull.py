# ============================================================
# Cristal FCC 3-D · Lennard-Jones · baño NVT
# Tracción uniaxial con fix addforce para medir E(S,X)
# ============================================================
from lammps import lammps
import numpy as np

# ---------- Parámetros ----------
T_target = 100.0
rho      = 0.8442
n_cells  = 6           # 6×6×6 celdas FCC → 864 átomos
steps_eq = 20000       # equilibración
steps_pr = 50000       # producción
dt       = 0.002
dumpstr  = 1000
f_atom   = 0.1         # fuerza externa por átomo (en ε/σ)

lmp = lammps()
c   = lmp.command

# ---------- Construcción y potencial ----------
c("units lj")
c("atom_style atomic")
c("boundary p p p")
c(f"lattice fcc {rho}")
c(f"region box block 0 {n_cells} 0 {n_cells} 0 {n_cells}")
c("create_box 1 box")
c("create_atoms 1 box")
c("pair_style lj/cut 2.5")
c("pair_coeff 1 1 1.0 1.0 2.5")
c("mass 1 1.0")

# ---------- Velocidades y minimización ----------
c("velocity all create 1.0 12345")
c(f"timestep {dt}")
c("minimize 1e-4 1e-6 100 1000")

# ---------- Baño NVT ----------
c(f"fix nvt all nvt temp 1.0 {T_target} 2.0")

# ---------- Tracción uniaxial ----------
c("group face region box")               # todos los átomos (demo)
c(f"fix pull face addforce {f_atom} 0.0 0.0")

# ---------- Variables para E ----------
c("variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)")   # |f_i|
c("compute  Fsum all reduce sum v_fmag")          # Σ|f_i|
c("variable f_sum equal c_Fsum")
c("variable f_obs equal f_pull")                  # fuerza externa total

# ---------- Dumps y termodinámica ----------
c(f"dump d1 all custom {dumpstr} dump.lj id type x y z fx fy fz")
c("thermo 500")

# ---------- Equilibración ----------
c(f"run {steps_eq}")

# ---------- Producción ----------
samples = steps_pr // dumpstr
E_vals  = np.zeros(samples)

for i in range(samples):
    c(f"run {dumpstr}")
    f_sum = lmp.extract_variable("f_sum", None, 0)
    f_obs = lmp.extract_variable("f_obs", None, 0)
    E_vals[i] = f_obs / f_sum

np.savetxt("E_vs_t.txt", E_vals)
print(f"E medio  = {E_vals.mean():.4f}")
print(f"Desviac. = {E_vals.std():.4f}")
