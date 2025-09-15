import subprocess as sp
import os
import numpy as np
import time
import shutil

def Optimizer(R, F, alpha=0.015):
    return (R + alpha*F)

def ReadOutcar(pathOutcar):
    with open(pathOutcar, "r") as f:
        lines = f.readlines()
        start = end = None
        for index,line in enumerate(lines):
            if "TOTAL-FORCE" in line:
                start = index + 2
            if "total drift" in line:
                end = index - 1
        if start is None or end is None:
            raise ValueError(f"Force data not found in OUTCAR: {pathOutcar}")
        RF = np.loadtxt(lines[start:end])
    return (RF[:,:3], RF[:,3:])

def Update(pathOutcarSinglet, pathOutcarTriplet, fTol = 0.01):
    R,F = ReadOutcar(os.path.expanduser(pathOutcarSinglet))
    rT,fT = ReadOutcar(os.path.expanduser(pathOutcarTriplet))
    fS = 2*F - fT
    status = "NA"
    maxForce = np.abs(fS).max()
    if  maxForce < fTol:
        status = "Converged"
        return (fS,R,maxForce,status)
    else:
        status = "Ongoing"
        return (fS,Optimizer(R, fS),maxForce,status)

def WritePOSCAR(pathPrevPoscar, pathOutcarSinglet, pathOutcarTriplet, pathSavePOSCAR):
    fS, R, maxForce, status = Update(pathOutcarSinglet, pathOutcarTriplet)
    with open(os.path.expanduser(pathPrevPoscar), "r") as f:
        lines = f.readlines()
        initialText = lines[:7]
    
    for path in pathSavePOSCAR:
        path = os.path.expanduser(path)
        with open(path, "w") as f:
            f.writelines(initialText)
            f.write("Cartesian\n")
            np.savetxt(f, R)
    return (fS, maxForce, status)

def ZPL(outcarGs, outcarEsSinglet, outcarEsTriplet):
    with open(outcarGs, "r") as f:
        lines = f.readlines()
        for line in lines:
            if "TOTEN" in line:
                Eg = float(line.strip().split()[-2])
    
    with open(outcarEsSinglet, "r") as f:
        lines = f.readlines()
        for line in lines:
            if "TOTEN" in line:
                Es = float(line.strip().split()[-2])
    
    with open(outcarEsTriplet, "r") as f:
        lines = f.readlines()
        for line in lines:
            if "TOTEN" in line:
                Et = float(line.strip().split()[-2])
    
    return (2*Es - Et - Eg)


path = "/scratch/pawsey1141/uqmsin17/cbcn_defect/rBN/bulk/pbe/defect/es/true_relaxation"
pathGS = "/scratch/pawsey1141/uqmsin17/cbcn_defect/rBN/bulk/pbe/defect/gs/scf"
status = "NA"
totIter = 150
shutil.copyfile(
    f"{path}/singlet/POSCAR",
    f"{path}/poscar/POSCAR0"
)
j1 = sp.run(["sbatch", "run.sh"], cwd=f"{path}/singlet", capture_output=True, text=True)
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
        pathOutcarSinglet=f"{path}/singlet/OUTCAR",
        pathOutcarTriplet=f"{path}/triplet/OUTCAR",
        pathSavePOSCAR=[f"{path}/singlet/POSCAR",f"{path}/triplet/POSCAR"]
    )

    shutil.copyfile(
    f"{path}/singlet/POSCAR",
    f"{path}/poscar/POSCAR{iter}"
)

    np.savetxt(f"{path}/forces/FORCES{iter}", fS)

    print(f"Iteration={iter} | Status={status} | MaximumForce={maxForce}")
    print(f"zpl = {ZPL(f'{pathGS}/OUTCAR', f'{path}/singlet/OUTCAR', f'{path}/triplet/OUTCAR')} eV")
    if status == "Ongoing":
        j1 = sp.run(["sbatch", "run.sh"], cwd=f"{path}/singlet", capture_output=True, text=True)
        j2 = sp.run(["sbatch", "run.sh"], cwd=f"{path}/triplet", capture_output=True, text=True)
        continue
    elif status == "Converged":
        print("Convergence Achieved, excited state geometry obtained.")
        break
    elif iter == totIter:
        print("Maximum iterations reached without convergence.")
        break

    