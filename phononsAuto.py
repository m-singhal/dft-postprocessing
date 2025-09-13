import subprocess as sp
import os
import time 


pathPhonons = os.path.expanduser("/scratch/pawsey1141/uqmsin17/cbcn_defect/rBN/bulk/pbe/defect/phonons/PHONOPY")
pathWavecar = os.path.expanduser("/scratch/pawsey1141/uqmsin17/cbcn_defect/rBN/bulk/pbe/defect/gs/scf/WAVECAR")
pathPhonopy = os.path.expanduser("/software/projects/pawsey1141/uqmsin17/myenv/bin/phonopy")
supercellDim = "1 1 1"

createDisps = sp.run([pathPhonopy,"-d",f"--dim={supercellDim}"], cwd = pathPhonons, capture_output=True, text=True)
createDisps = createDisps.stdout.strip().split()

numDisps = int(createDisps[createDisps.index("supercells:")+1])

if numDisps >= 1000: 
    sp.run(
        "for i in $(seq -w 1 999); do mv POSCAR-$i POSCAR-0$i; done",
        shell=True, cwd=pathPhonons
    )

sp.run(
    f"for i in $(seq -w 1 {numDisps}); do mkdir $i; mv POSCAR-$i ./$i/POSCAR; cp {pathPhonons}/INCAR {pathPhonons}/KPOINTS {pathPhonons}/POTCAR {pathPhonons}/run.sh ./$i/; ln -sf {pathWavecar} ./$i/; done",
    shell=True, cwd=pathPhonons
)

width = len(str(numDisps))
dollar_i = [f"{i:0{width}d}" for i in range(1, numDisps + 1)] 

if numDisps <= 500:
    for i in dollar_i:
        sp.run(
                "sbatch run.sh",
                shell=True, cwd=f"{pathPhonons}/{i}"
            )
        print(f"Submitted {i}!")
    print("\n")
    print(f"Jobs completed.")
        
else:
    nbatch = (numDisps // 500)
    dollar_iter = [] 
    for i in range(nbatch):
        dollar_iter.append(dollar_i[i*500:(i+1)*500])
    if nbatch < (numDisps/500):
        dollar_iter.append(dollar_i[nbatch*500:numDisps])
    
    for i in dollar_iter:
        JobID = []
        for j in i:
            submit = sp.run(
                    "sbatch run.sh",
                    shell=True, cwd=f"{pathPhonons}/{j}", capture_output=True, text=True
                )
            print(f"Submitted {j}!")
            JobID.append(submit.stdout.strip().split()[-1])
            stat = sp.run(["squeue", "-u", "uqmsin17"], capture_output=True, text=True)
        while any(job in stat.stdout.split() for job in JobID):
            time.sleep(30)
            stat = sp.run(["squeue", "-u", "uqmsin17"], capture_output=True, text=True)
        print("\n")
        print(f"Jobs completed till {i[-1]}.")
       
sp.run(f"{pathPhonopy} -f */vasprun.xml", shell=True, cwd=pathPhonons)
sp.run(f"{pathPhonopy} -p -s band.conf", shell=True, cwd=pathPhonons)