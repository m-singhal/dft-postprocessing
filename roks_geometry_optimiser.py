import subprocess as sp
import os
import numpy as np
import time
import shutil

def optimizer(R, F, alpha=0.015):
    return (R + alpha*F)

def read_outcar(path_outcar):
    with open(path_outcar, "r") as f:
        lines = f.readlines()
        start = end = None
        for index,line in enumerate(lines):
            if "TOTAL-FORCE" in line:
                start = index + 2
            if "total drift" in line:
                end = index - 1
        if start is None or end is None:
            raise ValueError(f"Force data not found in OUTCAR: {path_outcar}")
        RF = np.loadtxt(lines[start:end])
    return (RF[:,:3], RF[:,3:])

def update_structure(path_roks_outcar, f_tol = 0.01):
    R,F = read_outcar(os.path.expanduser(path_roks_outcar))
    status = "NA"
    maxForce = np.abs(F).max()
    if  maxForce < f_tol:
        status = "Converged"
        return (F,R,maxForce,status)
    else:
        status = "Ongoing"
        return (F, optimizer(R, F),maxForce,status)

def write_poscar(path_prev_poscar, path_roks_outcar, path_save_poscar):
    F, R, maxForce, status = update_structure(path_roks_outcar)
    with open(os.path.expanduser(path_prev_poscar), "r") as f:
        lines = f.readlines()
        initial_text = lines[:7]
    
    for path in path_save_poscar:
        path = os.path.expanduser(path)
        with open(path, "w") as f:
            f.writelines(initial_text)
            f.write("Cartesian\n")
            np.savetxt(f, R)
    return (F, maxForce, status)

def ZPL(path_gs_outcar, path_roks_outcar):
    with open(path_gs_outcar, "r") as f:
        lines = f.readlines()
        for line in lines:
            if "TOTEN" in line:
                Eg = float(line.strip().split()[-2])
    
    with open(path_roks_outcar, "r") as f:
        lines = f.readlines()
        for line in lines:
            if "TOTEN" in line:
                Eroks = float(line.strip().split()[-2])
        
    return (Eroks - Eg) 


path_calculation_roks = "/scratch/pawsey1141/uqmsin17/cbcn_defect/rBN/bulk/pbe/defect/es/true_relaxation"
path_gs_outcar = "/scratch/pawsey1141/uqmsin17/cbcn_defect/rBN/bulk/pbe/defect/gs/scf"
status = "NA"
totIter = 150
os.makedirs(f"{path_calculation_roks}/poscar", exist_ok=True)
shutil.copyfile(
    f"{path_calculation_roks}/01/POSCAR",
    f"{path_calculation_roks}/poscar/POSCAR0"
)
j0 = sp.run(["sbatch", "run.sh"], cwd=f"{path_calculation_roks}", capture_output=True, text=True)

for iter in range(1,totIter+1):
    stat = sp.run(["squeue", "-u", "uqmsin17"],
    cwd=f"{path_calculation_roks}",
    capture_output=True,
    text=True) 
    print(stat.stdout)
    job_id = j0.stdout.strip().split()[-1]
    while job_id in stat.stdout:
        print(f"Iteration={iter} | Job IDs={job_id} | VASP jobs are still running, waiting 30 seconds...")
        time.sleep(30)
        stat = sp.run(["squeue", "-u", "uqmsin17"],
        cwd=f"{path_calculation_roks}",
        capture_output=True,
        text=True)
        
    F,maxForce,status = write_poscar(
        path_prev_poscar=f"{path_calculation_roks}/poscar/POSCAR{iter-1}",
        path_roks_outcar=f"{path_calculation_roks}/01/OUTCAR",
        path_save_poscar=[f"{path_calculation_roks}/01/POSCAR",f"{path_calculation_roks}/02/POSCAR"]
    )

    shutil.copyfile(
    f"{path_calculation_roks}/01/POSCAR",
    f"{path_calculation_roks}/poscar/POSCAR{iter}"
)

    os.makedirs(f"{path_calculation_roks}/forces", exist_ok=True)
    np.savetxt(f"{path_calculation_roks}/forces/FORCES{iter}", F)

    print(f"Iteration={iter} | Status={status} | MaximumForce={maxForce}")
    print(f"zpl = {ZPL(f'{path_gs_outcar}/OUTCAR', f'{path_calculation_roks}/01/OUTCAR')} eV")
    if status == "Ongoing":
        j0 = sp.run(["sbatch", "run.sh"], cwd=f"{path_calculation_roks}", capture_output=True, text=True)
        continue
    elif status == "Converged":
        print("Convergence Achieved, excited state geometry obtained.")
        break
    elif iter == totIter:
        print("Maximum iterations reached without convergence.")
        break

    