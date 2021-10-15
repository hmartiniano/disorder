all:
	nohup snakemake -s Snakefile --cores 32 --use-singularity > snakemake.log &

clean:
	rm *.trr *.edr *.gro *.log *.tpr *.cpt *.top *.itp *.mdp

