from Bio import Phylo
import io
import glob


outfile="all_scogs_concatenated.tree"
with open(outfile, "w") as out:
	files=glob.glob("hybrid_orthofinder/21_genetrees/*.contree")
	for file in files:
		tree=Phylo.read("hybrid_orthofinder/21_genetrees/OG0006220.contree", "newick")
		for branch in tree.get_terminals():
			rename_branch=str(branch).split("_")[1]
			if rename_branch.startswith("Pp"):

				rename_branch=rename_branch.replace("Pp2192","CC2192")
				rename_branch=rename_branch.replace("Pp2298","PK2298")
			branch.name=rename_branch
		
		out.write(tree.format("newick"))