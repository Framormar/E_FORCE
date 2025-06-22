from lammps import lammps; import numpy as np
T,rho,n,dt,eq,pr,dump,fatom = 100.0,0.8442,6,0.002,20000,400000,4000,0.05
l=lammps(); c=l.command
c("units lj"); c("atom_style atomic"); c("boundary p p p")
c(f"lattice fcc {rho}"); c(f"region box block 0 {n} 0 {n} 0 {n}")
c("create_box 1 box"); c("create_atoms 1 box")
c("pair_style lj/cut 2.5"); c("pair_coeff 1 1 1 1 2.5"); c("mass 1 1")
c("velocity all create 1.0 12345"); c(f"timestep {dt}")
c("minimize 1e-4 1e-6 100 1000"); c(f"fix nvt all nvt temp 1.0 {T} 2.0")
c(f"variable f_atom equal {fatom}")
c("group face region box"); c("fix pull face addforce ${f_atom} 0 0")
c("variable Nface equal count(face)")
c('variable fx_int atom "fx - v_f_atom*(gmask(face)>0)"')
c('variable fmag_i atom "abs(v_fx_int)"')
c("compute Fsum all reduce sum v_fmag_i"); c("variable f_sum equal c_Fsum")
c('variable f_obs equal "v_f_atom * v_Nface"')
c(f"dump d1 all custom {dump} dump.lj id type x y z fx")
c("thermo 1000"); c(f"run {eq}")
E=np.zeros(pr//dump)
for i in range(len(E)):
    c(f"run {dump}")
    f_sum=float(l.extract_variable("f_sum",None,0))
    f_obs=float(l.extract_variable("f_obs",None,0))
    E[i]=f_obs/f_sum
np.savetxt("E_vs_t.txt",E)
print("E medio =",E.mean(),"σ =",E.std())
