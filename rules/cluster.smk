# vim: set ts=4 retab
rule run_cluster:
	input: trr="prd.trr", tpr="prd.tpr"
	output: protected("clust-id.xvg"), protected("clust-size.xvg"), protected("cluster.log"), protected("clusters.pdb"), protected("rmsd-clust.xpm"), protected("rmsd-dist.xvg")
	shell:
		"gmx cluster -f prd.trr -s prd.tpr -dist rmsd-dist.xvg -sz clust-size.xvg -clid clust-id.xvg -cl clusters.pdb"

"""
rule extract_cluster:
	input: trr="prd.trr", tpr="prd.tpr", clid="clust-id.xvg"
	output:
    shell:
		"gmx extract-cluster -f prd.trr -s prd.tpr -clusters {input.clid} -o clusters.pdb"
"""
