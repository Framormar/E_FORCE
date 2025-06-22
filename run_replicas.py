from lammps import lammps
import numpy as np, random, sys

rho = 0.8442; n = 6; dt = 0.002; eq = 15000
f_atom = 0.05; T = 100.0
seeds = [123, 456, 789, 321, 654]          # 5 réplicas

E_vals = []
for seed in seeds:
    l = lammps(); c = l.command
    c("units lj"); c("atom_style atomic"); c("boundary p p p")
    c(f"lattice fcc {rho}")
    c(f"region box block 0 {n} 0 {n} 0 {n}")
    c("create_box 1 box"); c("create_atoms 1 box")
    c("pair_style lj/cut 2.5"); c("pair_coeff 1 1 1 1 2.5")
    c("mass 1 1"); c(f"velocity all create 1.0 {seed}")
    c(f"timestep {dt}")
    c("minimize 1e-4 1e-6 100 1000")
    c(f"fix nvt all nvt temp 1.0 {T} 2.0")
    c(f"variable f_atom equal {f_atom}")
    c("fix pull all addforce ${f_atom} 0 0")
    c('variable fx_int atom "fx - v_f_atom"')
    c('variable fmag atom "abs(v_fx_int)"')
    c("compute Fsum all reduce sum v_fmag")
    c("variable f_sum equal c_Fsum")
    c('variable f_obs equal "v_f_atom * count(all)"')
    c("thermo 0"); c(f"run {eq}")
    E = l.extract_variable("f_obs",None,0)/l.extract_variable("f_sum",None,0)
    E_vals.append(E); l.close()

E_vals = np.array(E_vals)
print("E =", E_vals.mean(), "±", E_vals.std())

