"""
Mide el operador emergente Ey bajo cizalla (fuerza en +y) en un sólido FCC
Lennard-Jones.  Devuelve Ey.

ρ = 0.8442, T = 100 K (unidades LJ)
"""
from lammps import lammps

rho, n = 0.8442, 6         # densidad, celdas
dt, eq = 0.002, 15000
tau   = 0.05               # fuerza/átomo en +y

l = lammps(); c = l.command
c("units lj")
c("atom_style atomic")
c("boundary p p p")
c(f"lattice fcc {rho}")
c(f"region box block 0 {n} 0 {n} 0 {n}")
c("create_box 1 box")
c("create_atoms 1 box")
c("pair_style lj/cut 2.5")
c("pair_coeff 1 1 1 1 2.5")
c("mass 1 1")
c("velocity all create 1.0 12345")
c(f"timestep {dt}")
c("minimize 1e-4 1e-6 100 1000")
c(f"fix nvt all nvt temp 1.0 100.0 2.0")

# Cizalla: esfuerzo tau en +y
c(f"variable tau equal {tau}")
c("fix pull all addforce 0 ${tau} 0")

# Opera­dor E en componente y
c('variable fy_int atom "fy - v_tau"')
c('variable fmag atom "abs(v_fy_int)"')
c("compute Fsum all reduce sum v_fmag")
c('variable f_sum equal c_Fsum')
c('variable f_obs equal "v_tau * count(all)"')

c("thermo 0")
c(f"run {eq}")

E_y = l.extract_variable("f_obs",None,0) / l.extract_variable("f_sum",None,0)
print("Ey (cizalla) =", E_y)
