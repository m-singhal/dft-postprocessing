import subprocess as sp
import os
import numpy as np
import time
import shutil

def Optimizer(R, F, alpha=0.001):
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

def Update(pathOutcar, fTol = 0.01):
    R,F = ReadOutcar(os.path.expanduser(pathOutcar))
    status = "NA"
    maxForce = np.abs(F).max()
    if  maxForce < fTol:
        status = "Converged"
        return (R,maxForce,status)
    else:
        status = "Ongoing"
        return (Optimizer(R, F),maxForce,status)

def WritePOSCAR(pathPrevPoscar, pathOutcar, pathSavePOSCAR):
    R, maxForce, status = Update(pathOutcar)
    with open(os.path.expanduser(pathPrevPoscar), "r") as f:
        lines = f.readlines()
        initialText = lines[:7]
    
    for path in pathSavePOSCAR:
        path = os.path.expanduser(path)
        with open(path, "w") as f:
            f.writelines(initialText)
            f.write("Cartesian\n")
            np.savetxt(f, R)
    return (maxForce, status)

status = "NA"
totIter = 100
shutil.copyfile(
    "/scratch/pawsey1141/uqmsin17/cbcn_defect/trilayer/gs/my_relax/init/POSCAR",
    f"/scratch/pawsey1141/uqmsin17/cbcn_defect/trilayer/gs/my_relax/poscar/POSCAR0"
)
j1 = sp.run(["sbatch", "run.sh"], cwd="/scratch/pawsey1141/uqmsin17/cbcn_defect/trilayer/gs/my_relax/init", capture_output=True, text=True)

for iter in range(1,totIter+1):
    stat = sp.run(["squeue", "-u", "uqmsin17"],
    cwd="/scratch/pawsey1141/uqmsin17/cbcn_defect/trilayer/gs/my_relax/init",
    capture_output=True,
    text=True) 
    print(stat.stdout)
    jobID = j1.stdout.strip().split()[-1]
    while jobID in stat.stdout:
        print(f"Iteration={iter} | Job IDs={jobID} | VASP jobs are still running, waiting 30 seconds...")
        time.sleep(30)
        stat = sp.run(["squeue", "-u", "uqmsin17"],
        cwd="/scratch/pawsey1141/uqmsin17/cbcn_defect/trilayer/gs/my_relax/init",
        capture_output=True,
        text=True)
        

    
    maxForce,status = WritePOSCAR(
        pathPrevPoscar=f"/scratch/pawsey1141/uqmsin17/cbcn_defect/trilayer/gs/my_relax/poscar/POSCAR{iter-1}",
        pathOutcar="/scratch/pawsey1141/uqmsin17/cbcn_defect/trilayer/gs/my_relax/init/OUTCAR",
        pathSavePOSCAR=["/scratch/pawsey1141/uqmsin17/cbcn_defect/trilayer/gs/my_relax/init/POSCAR"]
    )
    shutil.copyfile(
    "/scratch/pawsey1141/uqmsin17/cbcn_defect/trilayer/gs/my_relax/init/POSCAR",
    f"/scratch/pawsey1141/uqmsin17/cbcn_defect/trilayer/gs/my_relax/poscar/POSCAR{iter}"
)
    
    print(f"Iteration={iter} | Status={status} | MaximumForce={maxForce}")
    if status == "Ongoing":
        j1 = sp.run(["sbatch", "run.sh"], cwd="/scratch/pawsey1141/uqmsin17/cbcn_defect/trilayer/gs/my_relax/init", capture_output=True, text=True)
        continue
    elif status == "Converged":
        print("Convergence Achieved, excited state geometry obtained.")
        break
    elif iter == totIter - 1:
        print("Maximum iterations reached without convergence.")
        break