# ============================================================
# Cristal FCC · Lennard-Jones · baño NVT
# Cálculo del operador E(S,X)  – versión clase lammps (sin PyLammps)
# ============================================================
from lammps import lammps
import numpy as np

# ----- parámetros -----
T_target = 100.0
rho      = 0.8442
n_cells  = 6
steps_eq = 20000
steps_pr = 50000
dt       = 0.002
dumpstr  = 1000

lmp = lammps()
c   = lmp.command           # alias

# ----- construcción -----
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
c("velocity all create 1.0 12345")
c(f"timestep {dt}")
c("minimize 1e-4 1e-6 100 1000")
c(f"fix nvt all nvt temp 1.0 {T_target} 2.0")

# ----- variables para E -----
c('variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)')      # magnitud |f_i|
c('compute Fsum all reduce sum v_fmag')              # Σ|f_i|
c('compute Fx   all reduce sum fx')                  # Σ f_ix
c('variable f_sum equal c_Fsum')
c('variable f_obs equal abs(c_Fx)')

c(f'dump d1 all custom {dumpstr} dump.lj id type x y z fx fy fz')
c('thermo 500')

# ----- equilibración -----
c(f'run {steps_eq}')

# ----- producción -----
samples = steps_pr // dumpstr
E_vals  = np.zeros(samples)

for i in range(samples):
    c(f'run {dumpstr}')
    f_sum = lmp.extract_variable('f_sum', None, 0)
    f_obs = lmp.extract_variable('f_obs', None, 0)
    E_vals[i] = f_obs / f_sum

np.savetxt('E_vs_t.txt', E_vals)
print(f'E medio  = {E_vals.mean():.4f}')
print(f'Desviac. = {E_vals.std():.4f}')
