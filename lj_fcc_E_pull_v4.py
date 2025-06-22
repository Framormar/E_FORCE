from lammps import lammps
import numpy as np

# -------- parámetros --------
T_target = 100.0
rho, n_cells = 0.8442, 6
dt, steps_eq, steps_pr, dumpstr = 0.002, 20000, 400000, 4000
f_atom = 0.05                      # fuerza externa/átomo (+x)

lmp = lammps()
c   = lmp.command                  # alias correcto

# --- sistema LJ FCC ---
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

# --- tracción uniaxial ---
c(f"variable f_atom equal {f_atom}")
c("group face region box")
c("fix pull face addforce ${f_atom} 0.0 0.0")

# --- variables internas para E ---
c('variable fx_int atom "fx - v_f_atom"')
c('variable fy_int atom "fy"')
c('variable fz_int atom "fz"')
c('variable fmag_i atom "sqrt(v_fx_int*v_fx_int + v_fy_int*v_fy_int + v_fz_int*v_fz_int)"')

c("compute  Fsum all reduce sum v_fmag_i")
c("variable f_sum equal c_Fsum")
c("variable f_obs equal abs(f_pull)")

# --- salida ---
c(f"dump d1 all custom {dumpstr} dump.lj id type x y z fx fy fz")
c("thermo 1000")

# --- equilibración ---
c(f"run {steps_eq}")

# --- producción ---
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
