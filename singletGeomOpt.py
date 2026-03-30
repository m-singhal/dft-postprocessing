import subprocess as sp
import os
import numpy as np
import time
import shutil


# =========================
# OPTIMIZER
# =========================
def Optimizer(R, F, alpha=0.015):
    R = R + alpha*F
    return R

# =========================
# READ FORCES FROM OUTCAR
# =========================
def ReadOutcar(pathOutcar):
    with open(pathOutcar, "r") as f:
        lines = f.readlines()

    start = end = None

    for i, line in enumerate(lines):
        if "TOTAL-FORCE" in line:
            start = i + 2
        if "total drift" in line:
            end = i - 1

    if start is None or end is None:
        raise ValueError(f"Force data not found in OUTCAR: {pathOutcar}")

    RF = np.loadtxt(lines[start:end])

    return RF[:, :3], RF[:, 3:]

# =========================
# UPDATE GEOMETRY
# =========================
def Update(pathOutcarMixed, pathOutcarTriplet, fTol=0.002):
    R, F = ReadOutcar(os.path.expanduser(pathOutcarMixed))
    rT, fT = ReadOutcar(os.path.expanduser(pathOutcarTriplet))

    # Effective excited-state force
    fS = 2 * F - fT

    maxForce = np.abs(fS).max()

    if maxForce < fTol:
        return fS, R, maxForce, "Converged"
    else:
        Rnew = Optimizer(R, fS)
        return fS, Rnew, maxForce, "Ongoing"

# =========================
# WRITE POSCAR
# =========================
def WritePOSCAR(pathPrevPoscar, pathOutcarMixed, pathOutcarTriplet, pathSavePOSCAR):
    fS, R, maxForce, status = Update(pathOutcarMixed, pathOutcarTriplet)

    with open(os.path.expanduser(pathPrevPoscar), "r") as f:
        header = f.readlines()[:7]

    for path in pathSavePOSCAR:
        path = os.path.expanduser(path)
        with open(path, "w") as f:
            f.writelines(header)
            f.write("Cartesian\n")
            np.savetxt(f, R)

    return fS, maxForce, status

# =========================
# ENERGY + ZPL
# =========================
def ZPL(outcarGs, outcarEsMixed, outcarEsTriplet):

    def get_energy(path):
        with open(path, "r") as f:
            lines = f.readlines()
        for line in reversed(lines):
            if "free  energy   TOTEN" in line:
                return float(line.strip().split()[-2])
        raise ValueError(f"TOTEN not found in {path}")

    Eg = get_energy(outcarGs)
    Em = get_energy(outcarEsMixed)
    Et = get_energy(outcarEsTriplet)

    Es = 2 * Em - Et

    return Es, Es - Eg


# =========================
# Main
# =========================
path = os.path.dirname(os.path.abspath(__file__))
pathGS = f"{path}/gs_outcar"
status = "NA"
totIter = 200
os.makedirs(f"{path}/poscar", exist_ok=True)
shutil.copyfile(
    f"{path}/mixed/POSCAR",
    f"{path}/poscar/POSCAR0"
)
j1 = sp.run(["sbatch", "run.sh"], cwd=f"{path}/mixed", capture_output=True, text=True)
j2 = sp.run(["sbatch", "run.sh"], cwd=f"{path}/triplet", capture_output=True, text=True)

for iter in range(1,totIter+1):
    stat = sp.run(["squeue", "-u", "uqmsin17"],
    cwd=f"{path}",
    capture_output=True,
    text=True) 
    print(stat.stdout)
    jobID = [j1.stdout.strip().split()[-1],j2.stdout.strip().split()[-1]]
    while jobID[0] in stat.stdout or jobID[1] in stat.stdout:
        print(f"Iteration={iter} | Job IDs={jobID} | VASP jobs are still running, waiting 30 seconds...")
        time.sleep(30)
        stat = sp.run(["squeue", "-u", "uqmsin17"],
        cwd=f"{path}",
        capture_output=True,
        text=True)
        

    
    fS,maxForce,status = WritePOSCAR(
        pathPrevPoscar=f"{path}/poscar/POSCAR{iter-1}",
        pathOutcarMixed=f"{path}/mixed/OUTCAR",
        pathOutcarTriplet=f"{path}/triplet/OUTCAR",
        pathSavePOSCAR=[f"{path}/mixed/POSCAR",f"{path}/triplet/POSCAR"]
    )

    shutil.copyfile(
    f"{path}/mixed/POSCAR",
    f"{path}/poscar/POSCAR{iter}"
)

    os.makedirs(f"{path}/forces", exist_ok=True)
    np.savetxt(f"{path}/forces/FORCES{iter}", fS)

    print(f"Iteration={iter} | Status={status} | MaximumForce={maxForce}")
    total_energy, zpl = ZPL(f'{pathGS}/OUTCAR', f'{path}/mixed/OUTCAR', f'{path}/triplet/OUTCAR')
    print(f"Total energy = {total_energy} eV, zpl = {zpl} eV")
    if status == "Ongoing":
        j1 = sp.run(["sbatch", "run.sh"], cwd=f"{path}/mixed", capture_output=True, text=True)
        j2 = sp.run(["sbatch", "run.sh"], cwd=f"{path}/triplet", capture_output=True, text=True)
        continue
    elif status == "Converged":
        print("Convergence Achieved, excited state geometry obtained.")
        break
    elif iter == totIter:
        print("Maximum iterations reached without convergence.")
        break