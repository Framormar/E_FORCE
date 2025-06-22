from lammps import lammps

a0      = 4.08        # parámetro de red físico Au (Å)
cells   = 6           # 6×6×6 → 864 átomos
dt      = 0.001
eq      = 20000
f_atom  = 0.05        # eV/Å

l = lammps(); c = l.command
c("units metal")
c("atom_style atomic")
c("boundary p p p")
c(f"lattice fcc {a0}")                       # ← parámetro real
c(f"region box block 0 {cells} 0 {cells} 0 {cells}")
c("create_box 1 box")
c("create_atoms 1 box")

c("pair_style eam")
c("pair_coeff * * Au_u3.eam")                # SETFL, sin nombre
c("mass 1 196.97")

# vecinos: cutoff (5.5 Å) + skin (1 Å) < L/2 ⇒ one 4000 suficiente
c("neighbor 7.0 bin")
c("neigh_modify delay 0 every 1 check yes one 4000")

c("velocity all create 300 12345 mom yes rot no")
c(f"timestep {dt}")
c("minimize 1e-4 1e-6 100 1000")
c("fix nvt all nvt temp 300 300 0.2")

c(f"variable f_atom equal {f_atom}")
c("fix pull all addforce ${f_atom} 0 0")

c('variable fx_int atom "fx - v_f_atom"')
c('variable fmag atom "abs(v_fx_int)"')
c("compute Fsum all reduce sum v_fmag")
c('variable f_sum equal c_Fsum')
c('variable f_obs equal "v_f_atom * count(all)"')

c("thermo 0")
c(f"run {eq}")

E_x = l.extract_variable("f_obs", None, 0) / l.extract_variable("f_sum", None, 0)
print(f"Ex (Au, a0={a0} Å) = {E_x:.6e}")

l.close()
