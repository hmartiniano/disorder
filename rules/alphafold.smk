rule install_alphafold:
	output:
	shell:
		"wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabfold_linux.sh -O resources/ &&"
		"cd resources &&"
		"bash install_colabfold_linux.sh"
