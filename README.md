# DateFBD
Divergence time distribution from fossil ages and topologies

Software includes
 - 'dist'
	computes the divergence time distibution associated to a given clade from a set of possible trees, the fossil ages, and the diversification rates
 - 'tree'
	computes the divergence time distribution associated to a given set of nodes of a single tree,  the fossil ages, and the diversification rates
 - 'esti'
	estimates the speciation, extinction and fossilization rates from a set of possible trees, the fossil ages, and the diversification rates
 - 'draw'
	draws a single tree with the fossil ages

type
	'> make all'
	in a console opened on the src directory for compiling software.

Directory "src" contains the C sources of software

Directory "data" contains the dataset studied in "Exact distribution of divergence times from fossil ages and topologies":
	'CotylosauriaTree.phy' contains 100 equiparsimonious trees of Cotylosauria
	'CotylosauriaAges.csv' contains the fossil ages
	'Optimisation_Parameters' are the parameters of the numerical optimisation used for determining the speciation, extinction and fossilization rates


A complete description of the options of the programs is given below.

------------
 dist 
------------

--------------------------
REQUIREMENT

	'dist' requires the gsl and pthread libraries.

--------------------------
COMPILING

	Just type
	> make dist
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'dist' computes the divergence time distibution associated to a given clade from a set of possible trees, the fossil ages, and the diversification rates


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	dist - Computation of the divergence time distibution associated to a given clade from a set of possible trees, the fossil ages, and the diversification rates
	
SYNOPSIS
	dist [OPTIONS] <input Tree File> <input Fossil File> <input List Clade> [output File]

DESCRIPTION
	Compute the distribution of the divergence time associated to the clade corresponding to the list of tips given in <input List Clade> by sampling into the trees contianed in <input Tree File> (in Newick format) into the fossil ranges provided in <input Fossil File> (in csv format) and output the distribution as a .csv table <output File>.csv

	Options are
	-z <input Tree File>
		output the tree in text debug format in the console and exit 
	-o <origin bound inf> [origin bound sup]
		set the origin range; if origin bound sup is not given it is set to the most ancient fossil age
	-e <end bound inf> <end bound sup>
		set the end time range
	-p <speciation rate> <extinction rate> <fossilization rate>
		set the speciation, extinction and fossilization rates
	-s <number>
		set the number of samples
	-d
		return the distribution (otherwise the density is returned by default)
	-u <value>
		set the step discretizing the time distribution to <value>
	-s <number>
		set the number of thread running in parallell
	-h
		display help

--------------------------


-------
 tree 
-------

--------------------------
REQUIREMENT

	'tree' requires the gsl, pthread and cairo libraries.

--------------------------
COMPILING

	Just type
	> make tree
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'tree' computes the divergence time distibutions associated to a given set of nodes from a single tree, the fossil ages, and the diversification rates and draw a tree figure


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	tree - Computation of the divergence time distibutions associated to a given set of nodes from a single tree, the fossil ages, and the diversification rates and draw a tree figure
SYNOPSIS
	tree [OPTIONS] <input Tree File> <input Fossil File> [output File]

DESCRIPTION
	Compute the divergence time distibutions associated to a given set of nodes from the tree of <input Tree File>, the fossil ages of <input Fossil File>, and the diversification rates and draw a tree figure

	Options are
	-z <input Tree File>
		output the tree in text debug format in the console and exit 
	-o <origin bound inf> [origin bound sup]
		set the origin range; if origin bound sup is not given it is set to the most ancient fossil age
	-e <end bound inf> <end bound sup>
		set the end time range
	-p <speciation rate> <extinction rate> <fossilization rate>
		set the speciation, extinction and fossilization rates
	-s <number>
		set the number of samples
	-d
		return the distribution (otherwise the density is returned by default)
	-u <value>
		set the step discretizing the time distribution to <value>
	-t <number>
		set the number of thread running in parallell
	-n <number>
		compute the divergence distribution associated to node <number>; option can be used several times; if it is not used all the divergence times are computed
	-x <number>
		set the graphic format of the output (option is required if one wants a graphic output)
			-x 1 -> pdf
			-x 2 -> postscript
			-x 3 -> png
			-x 4 -> svg
			-x 5 -> LaTeX (psTricks)
			-x 6 -> LaTeX (TikZ)
	-h
		display help

--------------------------


-------
 esti 
-------

--------------------------
REQUIREMENT

	'esti' requires the gsl, nlopt and pthread libraries.

--------------------------
COMPILING

	Just type
	> make esti
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'esti' estimates the speciation, extinction and fossilization rates from a set of trees, the fossil ages


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	esti - estimates the speciation, extinction and fossilization rates 
	
SYNOPSIS
	esti [OPTIONS] <input Tree File> <input Fossil File> <input List Clade> [output File]

DESCRIPTION
	Estimate the speciation, extinction and fossilization rates by sampling into the trees contained in <input Tree File> (in Newick format), into the fossil ranges provided in <input Fossil File> (in csv format) and output the result as a text file <output File>

	Options are
	-o <origin bound inf> [origin bound sup]
		set the origin range; if origin bound sup is not given it is set to the most ancient fossil age
	-e <end bound inf> <end bound sup>
		set the end time range
	-s <number>
		set the number of samples
	-t <number>
		set the number of thread running in parallell
	-n <option File>
		set the options of the numerical optimization to that provided in <option File>
	-h
		display help

--------------------------

-------
 draw 
-------

--------------------------
REQUIREMENT

	'draw' requires the cairo library.

--------------------------
COMPILING

	Just type
	> make draw
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'draw' draws a single tree with the fossil ages


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	draw - Drawing a single tree with the fossil ages
	
SYNOPSIS
	draw [OPTIONS] <input Tree File> <input Fossil File>

DESCRIPTION
	Output a figure in various graphic formats of the tree in <input Tree File> with the fossil ages of  <input Fossil File>

	Options are
	-z <input Tree File>
		output the tree in text debug format in the console and exit 
	-o <origin> 
		set the origin time
	-e <end> 
		set the end time
	-x <number>
		set the graphic format of the output (option is required if one wants a graphic output)
			-x 1 -> pdf
			-x 2 -> postscript
			-x 3 -> png
			-x 4 -> svg
			-x 5 -> LaTeX (psTricks)
			-x 6 -> LaTeX (TikZ)
	-h
		display help

--------------------------
