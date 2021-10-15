rule do_dssp:
	input: trr="prd.trr", tpr="prd.tpr"
	output: protected("scount.xvg"), protected("area.xpm")
	shell:
		"echo 1 | nohup gmx do_dssp -f {input.trr} -s {input.tpr} -sc scount.xvg -a area.xpm "
