"""
Barre cinco valores de tracción externa (0.02 … 0.10 ε/σ por átomo)
y registra E(f_obs) para un cristal FCC Lennard-Jones 3-D (NVT).
Salida: archivo E_vs_pull.txt con columnas f_atom, E_mean, E_std
"""
from lammps import lammps
import numpy as np

# --- parámetros fijos ---
T       = 100.0           # K (unidades LJ)
rho     = 0.8442          # densidad reducida
n_cells = 6               # 6×6×6 celdas → 864 átomos
dt      = 0.002
eq      = 10000           # pasos de equilibración
pr      = 80000           # pasos de producción
dump    = 4000            # stride de medición
f_scan  = np.linspace(0.02, 0.10, 5)   # fuerzas/átomo a probar

results = []

for f_atom in f_scan:
    lmp = lammps()
    c   = lmp.command

    # --- construcción del cristal ---
    c("units lj")
    c("atom_style atomic")
    c("boundary p p p")
    c(f"lattice fcc {rho}")
    c(f"region box block 0 {n_cells} 0 {n_cells} 0 {n_cells}")
    c("create_box 1 box")
    c("create_atoms 1 box")
    c("pair_style lj/cut 2.5")
    c("pair_coeff 1 1 1 1 2.5")
    c("mass 1 1")
    c("velocity all create 1.0 12345")
    c(f"timestep {dt}")
    c("minimize 1e-4 1e-6 100 1000")
    c(f"fix nvt all nvt temp 1.0 {T} 2.0")

    # --- tracción externa ---
    c(f"variable f_atom equal {f_atom}")
    c("group face region box")               # (todo el cristal para demo)
    c("fix pull face addforce ${f_atom} 0 0")
    c("variable Nface equal count(face)")

    # --- fuerzas internas (sólo módulo y sólo en 'face') ---
    c('variable fx_int atom "fx - v_f_atom*(gmask(face)>0)"')
    c('variable fmag_i atom "abs(v_fx_int)"')
    c('compute  Fsum face reduce sum v_fmag_i')
    c('variable f_sum equal c_Fsum')

    # --- fuerza observable ---
    c('variable f_obs equal "v_f_atom * v_Nface"')

    # --- termodinámica y corridas ---
    c("thermo 1000")
    c(f"run {eq}")

    samples = pr // dump
    E_vals  = np.zeros(samples)

    for i in range(samples):
        c(f"run {dump}")
        f_sum = lmp.extract_variable("f_sum", None, 0)
        f_obs = lmp.extract_variable("f_obs", None, 0)
        E_vals[i] = f_obs / f_sum

    results.append((f_atom, E_vals.mean(), E_vals.std()))
    lmp.close()

np.savetxt("E_vs_pull.txt", results, header="f_atom  E_mean  E_std")
print("✔️ Escaneo completado -> archivo E_vs_pull.txt")
