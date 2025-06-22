# ============================================================
# Cristal FCC 3-D · Lennard-Jones · baño NVT
# Tracción uniaxial suave → cálculo de E(S,X)
# ============================================================
from lammps import lammps
import numpy as np

# ---------- parámetros ----------
T_target = 100.0          # K (unidades LJ)
rho      = 0.8442         # densidad reducida
n_cells  = 6              # 6×6×6 celdas → 864 átomos
steps_eq = 20000          # equilibración
steps_pr = 200000         # producción (más largo para mejor estadística)
dt       = 0.002
dumpstr  = 2000           # guardar cada 2 000 pasos
f_atom   = 0.05           # fuerza por átomo en +x (ε/σ)

lmp = lammps()
c   = lmp.command

# ---------- construcción ----------
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

# ---------- arranque + minimización ----------
c("velocity all create 1.0 12345")
c(f"timestep {dt}")
c("minimize 1e-4 1e-6 100 1000")

# ---------- baño NVT ----------
c(f"fix nvt all nvt temp 1.0 {T_target} 2.0")

# ---------- tracción uniaxial ----------
c("group face region box")                # todos los átomos
c(f"fix pull face addforce {f_atom} 0.0 0.0")

# ---------- variables para E ----------
c("variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)")   # |f_i|
c("compute  Fsum all reduce sum v_fmag")          # Σ|f_i|
c("variable f_sum equal c_Fsum")
c("variable f_obs equal abs(f_pull)")             # |F_ext|

# ---------- dumps + termo ----------
c(f"dump d1 all custom {dumpstr} dump.lj id type x y z fx fy fz")
c("thermo 1000")

# ---------- equilibración ----------
c(f"run {steps_eq}")

# ---------- producción ----------
samples = steps_pr // dumpstr
E_vals  = np.zeros(samples)

for _ in range(samples):
    c(f"run {dumpstr}")
    f_sum = lmp.extract_variable('f_sum', None, 0)
    f_obs = lmp.extract_variable('f_obs', None, 0)
    E_vals[_] = f_obs / f_sum

np.savetxt("E_vs_t.txt", E_vals)
print(f"E medio  = {E_vals.mean():.4f}")
print(f"Desviac. = {E_vals.std():.4f}")
