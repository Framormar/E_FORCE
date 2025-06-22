"""
Barre T y P en un cristal FCC Lennard-Jones (fix npt) y mide
E(S,X) = F_obs / Σ |F_int|  para cada combinación (T,P).
"""
from lammps import lammps
import numpy as np

# ------------- ajustes globales -------------
rho       = 0.8442      # densidad reducida
n_cells   = 6           # 6×6×6 celdas → 864 átomos
dt        = 0.002
eq_steps  = 15000       # equilibración
prod_steps= 60000       # producción
dump      = 3000        # stride medición
f_atom    = 0.05        # fuerza externa por átomo (+x) ε/σ
T_list    = [10, 50, 100, 150, 200]   # K (unid. LJ)
P_list    = [0, 1, 2, 3, 5]           # ε σ⁻³

results = []

for T in T_list:
    for P in P_list:
        lmp = lammps(); c = lmp.command
        # -- cristal L-J --
        c("units lj"); c("atom_style atomic"); c("boundary p p p")
        c(f"lattice fcc {rho}")
        c(f"region box block 0 {n_cells} 0 {n_cells} 0 {n_cells}")
        c("create_box 1 box"); c("create_atoms 1 box")
        c("pair_style lj/cut 2.5"); c("pair_coeff 1 1 1 1 2.5")
        c("mass 1 1"); c("velocity all create 1.0 12345")
        c(f"timestep {dt}")
        c("minimize 1e-4 1e-6 100 1000")
        c(f"fix npt all npt temp {T} {T} 2.0 iso {P} {P} 5.0")
        # fuerza externa homogénea
        c(f"variable f_atom equal {f_atom}")
        c("fix pull all addforce ${f_atom} 0 0")
        # fuerzas internas (restar f_atom en +x)
        c('variable fx_int atom "fx - v_f_atom"')
        c('variable fmag_i atom "abs(v_fx_int)"')
        c('compute  Fsum all reduce sum v_fmag_i')
        c('variable f_sum equal c_Fsum')
        # fuerza observable = f_atom · Natoms
        c("variable Nat equal count(all)")
        c('variable f_obs equal "v_f_atom * v_Nat"')
        # termo y corridas
        c("thermo 1000")
        c(f"run {eq_steps}")
        samples = prod_steps // dump
        E_vals = np.zeros(samples)
        for i in range(samples):
            c(f"run {dump}")
            f_sum = lmp.extract_variable("f_sum",None,0)
            f_obs = lmp.extract_variable("f_obs",None,0)
            E_vals[i] = f_obs / f_sum
        results.append((T, P, E_vals.mean(), E_vals.std()))
        lmp.close()

np.savetxt("E_vs_TP.txt", results,
           header="T    P    E_mean    E_std", fmt="%.6g")
print("✔️  Escaneo T–P completado ⇒ archivo E_vs_TP.txt")
