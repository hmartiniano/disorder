rule gyration:
	input: trr="prd.trr", tpr="prd.tpr"
	output: protected("gyration.xvg")
	shell:
		"echo 1 | nohup gmx gyration -f {input.trr} -s {input.tpr} -o gyration.xvg"
