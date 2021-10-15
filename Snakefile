# vi: syntax=python
# vi: ft=python

import os

os.environ["OMP_NUM_THREADS"] = "32"


rule all:
  input:
    "prd.log"

##singularity: "docker://gromacs/gromacs"


rule convert_pdb:
  input: "protein.pdb"
  output: gro="protein.gro", top="topol.top"
  #singularity: "docker://gromacs/gromacs"
  shell:
    "echo '5\n3' | gmx pdb2gmx -f {input} -ignh -o {output.gro}"

rule setup_box:
  input: "protein.gro"
  output: "newbox.gro"
  #singularity: "docker://gromacs/gromacs"
  shell:
    "echo Protein | gmx editconf -f {input} -o {output} -princ -d 2"


rule solvate:
  input: gro="newbox.gro"
  output: gro="solv.gro"
  #singularity: "docker://gromacs/gromacs"
  shell: "gmx solvate -cp newbox.gro -cs tip4p.gro -o solv.gro -p topol.top"

#wget http://www.mdtutorials.com/gmx/umbrella/Files/ions.mdp

rule setup_add_ions:
  input: "solv.gro"
  output: "ions.tpr"
  #singularity: "docker://gromacs/gromacs"
  shell: "gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1"


rule add_ions:
  input: "ions.tpr"
  output: "solv_ions.gro"
  #singularity: "docker://gromacs/gromacs"
  shell: "echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.1"

#wget http://www.mdtutorials.com/gmx/umbrella/Files/minim.mdp
#wget http://www.mdtutorials.com/gmx/umbrella/Files/npt.mdp
#
rule setup_minimize:
  input: "solv_ions.gro"
  output: "em.tpr"
  #singularity: "docker://gromacs/gromacs"
  shell: "gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 10"


rule minimize:
  input: "em.tpr"
  output: protected("em.gro"), protected("em.log"), protected("em.edr")
  #singularity: "docker://gromacs/gromacs"
  threads: 8
  shell: "gmx mdrun -v -deffnm $(basename {input} .tpr) -ntmpi 1 -ntomp {threads}"


rule setup_equilibrium1:
  input: mdp="eq1.mdp", gro="em.gro"
  output: "eq1.tpr"
  #singularity: "docker://gromacs/gromacs"
  shell: "gmx grompp -f {input.mdp} -c {input.gro} -p topol.top -r {input.gro} -o {output} -maxwarn 10"


rule equilibrium1:
  input: "eq1.tpr"
  output: gro=protected("eq1.gro"),
          cpt=protected("eq1.cpt"),
          trr=protected("eq1.trr"),
          log=protected("eq1.log")
  threads: 8
  #singularity: "docker://gromacs/gromacs"
  shell: "gmx mdrun -deffnm $(basename {input} .tpr) -ntmpi 1 -ntomp {threads}"

rule setup_equilibrium2:
  input: mdp="eq2.mdp", gro="eq1.gro"
  output: "eq2.tpr"
  #singularity: "docker://gromacs/gromacs"
  shell: "gmx grompp -f {input.mdp} -c {input.gro} -p topol.top -r {input.gro} -o {output} -maxwarn 10"

rule equilibrium2:
  input: tpr="eq2.tpr", prev="eq1.cpt"
  output: gro=protected("eq2.gro"),
          cpt=protected("eq2.cpt"),
          trr=protected("eq2.trr"),
          log=protected("eq2.log")
  threads: 8
  #singularity: "docker://gromacs/gromacs"
  shell: "gmx mdrun -deffnm $(basename {input.tpr} .tpr) -ntmpi 1 -ntomp {threads}"

rule setup_production:
  input: mdp="prd.mdp", gro="eq2.gro"
  output: "prd.tpr"
  #singularity: "docker://gromacs/gromacs"
  shell: "gmx grompp -f {input.mdp} -c {input.gro} -p topol.top -r {input.gro} -o {output} -maxwarn 10"


rule production:
  input: "prd.tpr"
  output: gro=protected("prd.gro"), cpt=protected("prd.cpt"), trr=protected("prd.trr"), log=protected("prd.log")
  threads: 8
  #singularity: "docker://gromacs/gromacs"
  shell: "gmx mdrun -deffnm prd "


