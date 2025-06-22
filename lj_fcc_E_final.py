from lammps import lammps
import numpy as np

# ---------------- PARÁMETROS ----------------
T_target = 100.0            # K  (unidades LJ)
rho, n_cells = 0.8442, 6    # densidad y celdas FCC
dt = 0.002
steps_eq, steps_pr, dumpstr = 20_000, 400_000, 4_000
f_atom = 0.05               # fuerza externa por átomo (+x) en ε/σ

# ---------------- INICIALIZACIÓN ------------
lmp = lammps()
c   = lmp.command            # alias para enviar comandos

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

# ---------------- TRACCIÓN ------------------
c(f"variable f_atom equal {f_atom}")       # fuerza por átomo
c("group face region box")                 # (usa todo el cristal por simplicidad)
c("fix pull face addforce ${f_atom} 0.0 0.0")

# nº de átomos bajo tracción
c("variable Nface equal count(face)")

# ---------------- FUERZAS INTERNAS ----------
c('variable fx_int atom "fx - v_f_atom"')
c('variable fmag_i atom "abs(v_fx_int)"')           # |Fx interna|
c("compute  Fsum all reduce sum v_fmag_i")
c("variable f_sum equal c_Fsum")

# ---------------- FUERZA OBSERVABLE ---------
c('variable f_obs equal "v_f_atom * v_Nface"')       # F_ext total

# ---------------- SALIDA --------------------
c(f"dump d1 all custom {dumpstr} dump.lj id type x y z fx")
c("thermo 1000")

# ---------------- EQUILIBRACIÓN -------------
c(f"run {steps_eq}")

# ---------------- PRODUCCIÓN ----------------
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
