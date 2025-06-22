"""
Calcula las componentes del operador E (Ex, Ey, Ez)
para tres configuraciones de un cristal FCC Lennard-Jones:

  tag 0 → perfecto
  tag 1 → 1 % de vacancias aleatorias
  tag 2 → dislocación de borde creada manualmente
           (se desplaza la mitad superior media celda en +x)

Genera defects_E_manual.txt con columnas:
tag  Ex  Ey  Ez
"""
import numpy as np, random
from lammps import lammps

# Parámetros globales
T, rho, n_cells   = 100.0, 0.8442, 6
dt, eq, prod, st  = 0.002, 15000, 60000, 3000
f_atom            = 0.05   # fuerza/átomo (+x,+y,+z)

def base_system():
    L = lammps(); c = L.command
    c("units lj"); c("atom_style atomic"); c("boundary p p p")
    c(f"lattice fcc {rho}")
    c(f"region box block 0 {n_cells} 0 {n_cells} 0 {n_cells}")
    c("create_box 1 box"); c("create_atoms 1 box")
    c("pair_style lj/cut 2.5"); c("pair_coeff 1 1 1 1 2.5")
    c("mass 1 1"); c("velocity all create 1.0 12345")
    return L, c

def add_common(c):
    c(f"timestep {dt}")
    c("minimize 1e-4 1e-6 100 1000")
    c(f"fix nvt all nvt temp 1.0 {T} 2.0")
    # fuerza homogénea en x,y,z
    c(f"variable f_atom equal {f_atom}")
    c("fix pull all addforce ${f_atom} ${f_atom} ${f_atom}")
    c("variable Nat equal count(all)")
    # define E en cada componente
    for ax, fx in zip("xyz", ("fx", "fy", "fz")):
        c(f'variable f{ax}_int atom "{fx} - v_f_atom"')
        c(f'variable fmag_{ax} atom "abs(v_f{ax}_int)"')
        c(f'compute  F{ax} all reduce sum v_fmag_{ax}')
        c(f'variable f_sum_{ax} equal c_F{ax}')
        c(f'variable f_obs_{ax} equal "v_f_atom * v_Nat"')

def measure_E(L):
    c = L.command
    c("thermo 1000"); c(f"run {eq}")
    E = {ax: [] for ax in "xyz"}
    for _ in range(prod // st):
        c(f"run {st}")
        for ax in "xyz":
            obs = L.extract_variable(f"f_obs_{ax}", None, 0)
            denom = L.extract_variable(f"f_sum_{ax}", None, 0)
            E[ax].append(obs / denom)
    return tuple(np.mean(E[ax]) for ax in "xyz")

# Almacena resultados
records = []

# 0: perfecto
L, c = base_system()
add_common(c)
records.append((0, *measure_E(L)))
L.close()

# 1: vacancias 1 %
L, c = base_system()
nat = n_cells**3 * 4
vac_ids = random.sample(range(1, nat + 1), int(0.01 * nat))
c("group vac id " + " ".join(map(str, vac_ids)))
c("delete_atoms group vac")
add_common(c)
records.append((1, *measure_E(L)))
L.close()

# 2: dislocación manual
L, c = base_system()
# definimos región superior
c("variable Ly equal ly")
c("variable halfLy equal v_Ly/2.0")
c("region up block INF INF ${halfLy} INF INF INF units box")
c("group up region up")
# desplazamos media celda en +x
c("variable shift equal 0.5")
c("displace_atoms up move ${shift} 0 0 units box")
add_common(c)
records.append((2, *measure_E(L)))
L.close()

# Guardar
np.savetxt("defects_E_manual.txt", records,
           header="tag  Ex  Ey  Ez", fmt="%.6e")
print("✔️  Archivo defects_E_manual.txt generado")
