from lammps import lammps
import numpy as np

N = 8000            # nº átomos
L = (N/0.70)**(1/3) # caja cúbica para ρ=0.70
lmp = lammps(); c = lmp.command

c("units lj")
c("atom_style atomic")
c("boundary p p p")
c(f"region box block 0 {L} 0 {L} 0 {L}")
c("create_box 1 box")
c(f"create_atoms 1 random {N} 12345 box")
c("pair_style lj/cut 2.5")
c("pair_coeff 1 1 1 1 2.5")
c("mass 1 1")

# Minimización suave
c("pair_modify shift yes")          # evita energías ∞
c("minimize 1e-5 1e-6 1000 10000")
c("pair_modify shift no")

# Termostatado
c("velocity all create 0.5 9876 mom yes rot no")
c("neighbor 0.3 bin")
c("timestep 0.001")
c("fix nvt all nvt temp 0.5 1.2 2.0")

# Carga externa débil en +x
c("variable f_atom equal 0.01")
c("fix pull all addforce ${f_atom} 0 0")

# Operador E
c('variable fx_int atom "fx - v_f_atom"')
c('variable fmag atom "abs(v_fx_int)"')
c("compute F all reduce sum v_fmag")
c('variable f_sum equal c_F')
c('variable f_obs equal "v_f_atom * count(all)"')

c("thermo 0")
c("run 20000")      # 20 000 pasos con Δt=0.001 → 20 τ

E = lmp.extract_variable("f_obs", None, 0) / \
    lmp.extract_variable("f_sum", None, 0)
print("E_fluido =", E)
