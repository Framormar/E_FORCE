"""
Fase D: correlaciona el operador emergente E con el módulo de Young E_Young
en un cristal FCC Lennard-Jones 3-D.

Salida: E_vs_Young.txt con columnas
rho  E_emergent  E_Young
"""
import numpy as np
from lammps import lammps

# Parámetros globales
n_cells  = 6
T        = 100.0      # temperatura (unidades LJ)
dt       = 0.002
eq_steps = 15000      # pasos de equilibración
f_atom   = 0.05       # fuerza externa por átomo (+x)
strain   = 0.01       # deformación uniaxial para medir Young

rho_list = [0.84, 0.86, 0.88, 0.90]

results = []

for rho in rho_list:
    # ── 1) Medir E emergente ────────────────────────────────────────
    L = lammps(); c = L.command
    c("units lj"); c("atom_style atomic"); c("boundary p p p")
    c(f"lattice fcc {rho}")
    c(f"region box block 0 {n_cells} 0 {n_cells} 0 {n_cells}")
    c("create_box 1 box")
    c("create_atoms 1 box")
    c("pair_style lj/cut 2.5"); c("pair_coeff 1 1 1 1 2.5")
    c("mass 1 1"); c("velocity all create 1.0 12345")
    c(f"timestep {dt}")
    c("minimize 1e-4 1e-6 100 1000")
    c(f"fix nvt all nvt temp 1.0 {T} 2.0")
    # configurar medición de E emergente
    c(f"variable f_atom equal {f_atom}")
    c("group face region box")
    c("fix pull face addforce ${f_atom} 0 0")
    c("variable Nface equal count(face)")
    c('variable fx_int atom "fx - v_f_atom"')
    c('variable fmag atom "abs(v_fx_int)"')
    c("compute Fsum all reduce sum v_fmag")
    c("variable f_sum equal c_Fsum")
    c('variable f_obs equal "v_f_atom * v_Nface"')
    c("thermo 0")
    c(f"run {eq_steps}")
    f_sum = L.extract_variable("f_sum", None, 0)
    f_obs = L.extract_variable("f_obs", None, 0)
    E_em = f_obs / f_sum
    L.close()

    # ── 2) Medir módulo de Young con cambio de caja y get_thermo ───
    L = lammps(); c = L.command
    c("units lj"); c("atom_style atomic"); c("boundary p p p")
    c(f"lattice fcc {rho}")
    c(f"region box block 0 {n_cells} 0 {n_cells} 0 {n_cells}")
    c("create_box 1 box")
    c("create_atoms 1 box")
    c("pair_style lj/cut 2.5"); c("pair_coeff 1 1 1 1 2.5")
    c("mass 1 1"); c("velocity all create 1.0 12345")
    c(f"timestep {dt}")
    c("minimize 1e-4 1e-6 100 1000")
    # aplicar deformación uniaxial en x
    c(f"change_box all x scale {1+strain} y scale 1.0 z scale 1.0 remap units box")
    # configurar para sacar la componente Pxx
    c("thermo_style custom pxx")
    c("thermo 1")
    c("run 0")
    pxx = L.get_thermo("pxx")
    # Módulo de Young: E = -pxx / strain
    E_young = -pxx / strain
    L.close()

    results.append((rho, E_em, E_young))

# Guardar correlación
np.savetxt("E_vs_Young.txt", results,
           header="rho  E_emergent  E_Young", fmt="%.6e")
print("✅ Datos guardados en E_vs_Young.txt")
