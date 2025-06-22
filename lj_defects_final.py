"""
Calcula las tres componentes del operador E (Ex, Ey, Ez) para:

 0 – Cristal perfecto
 1 – Cristal con 1 % de vacancias (aleatorias)
 2 – Cristal con dislocación de borde (requiere USER-DISLOC)

Guarda los resultados en defects_E.txt con columnas:
tag  Ex  Ey  Ez
"""
from lammps import lammps
import numpy as np, random

# ---------- parámetros comunes ----------
T       = 100.0
rho     = 0.8442
n_cells = 6            # 6×6×6 celdas FCC → 864 átomos
dt      = 0.002
eq      = 15000        # pasos de equilibración
pr      = 60000        # pasos de producción
dump    = 3000         # stride medición
f_atom  = 0.05         # fuerza homogénea (+x,+y,+z) ε/σ

def build_system(tag):
    L = lammps(); c = L.command
    c("units lj"); c("atom_style atomic"); c("boundary p p p")
    c(f"lattice fcc {rho}")
    c(f"region box block 0 {n_cells} 0 {n_cells} 0 {n_cells}")
    c("create_box 1 box"); c("create_atoms 1 box")

    # ---- defectos según 'tag' ----
    if tag == 1:                             # 1 % vacancias
        nat = n_cells**3 * 4
        n_vac = int(0.01 * nat)
        ids = random.sample(range(1, nat+1), n_vac)
        c("group vac id " + " ".join(map(str, ids)))
        c("delete_atoms group vac")
    if tag == 2:                             # dislocación de borde (requiere USER-DISLOC)
        c("dislocation edge 0.5  0 0 1  1 0 0  1")   # Burgers b=0.5[100]

    # ---- potencial y dinámica ----
    c("pair_style lj/cut 2.5")
    c("pair_coeff 1 1 1 1 2.5")
    c("mass 1 1")
    c("velocity all create 1.0 12345")
    c(f"timestep {dt}")
    c("minimize 1e-4 1e-6 100 1000")
    c(f"fix nvt all nvt temp 1.0 {T} 2.0")

    # fuerza externa en x,y,z
    c(f"variable f_atom equal {f_atom}")
    c("fix pull all addforce ${f_atom} ${f_atom} ${f_atom}")
    c("variable Nat equal count(all)")

    # --- variables y computes para Ex,Ey,Ez ---
    for ax, fx in zip(("x","y","z"), ("fx","fy","fz")):
        # fuerza interna = total - fuerza externa aplicada
        c(f'variable f{ax}_int atom "{fx} - v_f_atom"')
        c(f'variable fmag_{ax} atom "abs(v_f{ax}_int)"')
        c(f'compute Fsum_{ax} all reduce sum v_fmag_{ax}')
        c(f'variable f_sum_{ax} equal c_Fsum_{ax}')
        c(f'variable f_obs_{ax} equal "v_f_atom * v_Nat"')
    return L

def measure_E(L):
    c = L.command
    c("thermo 1000"); c(f"run {eq}")
    samples = pr // dump
    E = {ax: [] for ax in ("x","y","z")}
    for _ in range(samples):
        c(f"run {dump}")
        for ax in ("x","y","z"):
            obs = L.extract_variable(f"f_obs_{ax}", None, 0)
            s   = L.extract_variable(f"f_sum_{ax}", None, 0)
            E[ax].append(obs / s)
    return (np.mean(E["x"]), np.mean(E["y"]), np.mean(E["z"]))

results = []
for tag, label in [(0,"perfecto"), (1,"vacancias"), (2,"disloc")]:
    print(f"▶️  Simulación: {label}")
    L = build_system(tag)
    Ex, Ey, Ez = measure_E(L)
    L.close()
    print(f"   E = ({Ex:.3e}, {Ey:.3e}, {Ez:.3e})")
    results.append((tag, Ex, Ey, Ez))

np.savetxt("defects_E.txt", results,
           header="tag  Ex  Ey  Ez", fmt="%.6e")
print("\n✔️  Archivo defects_E.txt guardado")
