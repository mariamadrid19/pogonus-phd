
git clone https://github.com/lh3/minigraph
cd minigraph && make

#Build pangenome using minigraph
minigraph -cxggs -L 100 -o pogonus.gfa sorted_prim_dud.fasta sorted_prim_nieu.fasta

#Validate the graph
gfatools pogonus.gfa
