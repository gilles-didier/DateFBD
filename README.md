# DateFBD
Divergence time distribution from fossil ages and topologies

Software includes
 - 'dspe'
	computes the divergence time distibution associated to a given clade from a set of possible trees, the fossil stratigraphic intervals and provides coda files for MCMC convergence diagnostics
 - 'tspe'
	computes the divergence time distribution associated to a given set of nodes of a single tree and the fossil stratigraphic intervals
 - 'dext'
	computes the extinction time distibutions associated to a given set of clades from a set of possible trees, the fossil stratigraphic intervals and provides coda files for MCMC convergence diagnostics
 - 'text'
	computes the extinction time distribution associated to all the terminal nodes of a single tree with the fossil stratigraphic intervals
 - 'qext'
	computes the quantile at a given order, of extinction time distibution associated to a given set of clades from a set of possible trees, the fossil stratigraphic intervals and provides coda files for MCMC convergence diagnostics
 - 'cext'
	computes the probability that a set of extinct taxa goes extinct before another one from a set of possible trees, the fossil stratigraphic intervals and provides coda files for MCMC convergence diagnostics
 - 'diag'
	same input as 'dspe', provides coda files for MCMC convergence diagnostics
 - 'esti'
	estimates (by maximum likelihood) the speciation, extinction and fossilization rates from a set of possible trees, the fossil stratigraphic intervals
 - 'draw'
	draws a single tree with the fossil stratigraphic intervals
 - 'samp'
	sample observable trees under the FBD model 
 - R package 'rfbd'
	implements 'dspe', 'dext', 'qext', 'cext' and 'diag' in the R environment (some features of the standalone software such as the multithreading are not available yet)

type
	> make all
	in a console opened on the src directory for compiling standalone software.
	
type	
	> install.packages("<path to>/rfbd_0.0.tar.gz", repos = NULL, type = "source")

	> library(rfbd)
	in a R console for loading the R package (where <path to> is the path to the file 'rfbd_0.0.tar.gz')

Directory "src" contains the C sources of the standalone software

Archive file "rfbd_0.0.tar.gz" contains the R package 'rfdb', including source files and manual

Directory "data" contains two subdirectories
	'Speciation' contains the dataset studied in "Exact distribution of divergence times from fossil ages and topologies":
		'CotylosauriaTree.phy' contains 100 equiparsimonious trees of Cotylosauria
		'CotylosauriaAges.csv' contains the fossil ages
		'Optimisation_Parameters' are the parameters of the numerical optimisation used for determining the speciation, extinction and fossilization rates
	'Extinction' contains the dataset studied in "Distributions of extinction times from fossil ages and tree topologies: the example of some mid-Permian synapsid extinctions":
		'TreeEupely100.phy' contains 100 equiparsimonious trees of Eupelicosauria
		'FossilAges.csv' contains the fossil ages
		'Edaphosauridae', 'Ophiacodontidae' and 'Sphenacodontidae' contain the so-called clades


A complete description of the options of the programs is given below.

------------
 dspe 
------------

--------------------------
REQUIREMENT

	'dspe' requires the gsl and pthread libraries.

--------------------------
COMPILING

	Just type
	> make dspe
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'dspe' computes the divergence time distibution associated to a given clade from a set of possible trees, the fossil ages, and the diversification rates


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	dspe - computes the divergence time distibution associated to a given clade from a set of possible trees, the fossil stratigraphic intervals and provides coda files for MCMC convergence diagnostics
	
SYNOPSIS
	dspe [OPTIONS] <input List Clade> <input Tree File> <input Fossil File> [output File]

DESCRIPTION
	Compute the distribution of the divergence time associated to the clade corresponding to the list of tips given in <input List Clade> by sampling into the trees contained in <input Tree File> (in Newick format) into the fossil ranges provided in <input Fossil File> (in csv format) and output the distribution as a .csv table <output File>.csv and two files <output File>.out and <output File>.ind which can be read by the R package coda for MCMC diagnostics

	Options are
	-z <input Tree File>
		output the tree in text debug format in the console and exit 
	-o <origin bound inf> [origin bound sup]
		set the origin range; if origin bound sup is not given it is set to the most ancient fossil age
	-w <speciation rate width> <extinction rate width> <fossilization rate width>
		set the width of the sliding windows used for samplint the speciation, extinction and fossilization rates during the MCMC
	-i <speciation rate upper bound> <extinction rate  upper bound> <fossilization rate upper bound>
		set the upper bounds of the interval in which the speciation, extinction and fossilization rates are uniformly drawn at the beginning of the MCMC (lower bounds are all 0)
	-b <number>
		set the number of iterations for the burning stage
	-g <number>
		set the thinning parameter (only one iteration in <number> is considered)  
	-f <proportion>
		set the proportion of moves in the rates in the MCMC proposal (other moves are drawn uniformly amont the fossil ages)
	-a <speciation proportion> <extinction proportion>
		set the relation proportion of moves in the speciation and the extinction rates (thus in the fossilization rate)
	-s <number>
		set the number of samples required
	-r <number>
		set the random seed
	-d
		return the distribution (otherwise the density is returned by default)
	-u <value>
		set the step discretizing the time distribution to <value>
	-t <number>
		set the number of thread running in parallell
	-h
		display help

EXAMPLE

./dspe -o -400 -s 2000 -u 0.1 -b 10000 -g 300 -r 5555 -f 0.25 -a 0.33 0.33 -w 0.1 0.1 0.01 -i 5. 5. 5. -t 40 Amniota.txt CotylosauriaTree.phy CotylosauriaAges.csv amniota2000_400


--------------------------

-------
 tspe 
-------

--------------------------
REQUIREMENT

	'tspe' requires the gsl, pthread and cairo libraries.

--------------------------
COMPILING

	Just type
	> make tspe
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'tspe' computes the divergence time distributions associated to a given set of nodes from a single tree, the fossil ages, and the diversification rates and draw a tree figure


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	tspe - Computation of the divergence time distributions associated to a given set of nodes from a single tree, the fossil ages, and the diversification rates and draw a tree figure
SYNOPSIS
	tspe [OPTIONS] <input Tree File> <input Fossil File> [output File]

DESCRIPTION
	Compute the divergence time distributions associated to a given set of nodes from the tree of <input Tree File>, the fossil ages of <input Fossil File>, and the diversification rates and draw a tree figure

	Options are
	-z <input Tree File>
		output the tree in text debug format in the console and exit 
	-o <origin bound inf> [origin bound sup]
		set the origin range; if origin bound sup is not given it is set to the most ancient fossil age
	-w <speciation rate width> <extinction rate width> <fossilization rate width>
		set the width of the sliding windows used for samplint the speciation, extinction and fossilization rates during the MCMC
	-i <speciation rate upper bound> <extinction rate  upper bound> <fossilization rate upper bound>
		set the upper bounds of the interval in which the speciation, extinction and fossilization rates are uniformly drawn at the beginning of the MCMC (lower bounds are all 0)
	-b <number>
		set the number of iterations for the burning stage
	-g <number>
		set the thinning parameter (only one iteration in <number> is considered)  
	-f <proportion>
		set the proportion of moves in the rates in the MCMC proposal (other moves are drawn uniformly amont the fossil ages)
	-a <speciation proportion> <extinction proportion>
		set the relation proportion of moves in the speciation and the extinction rates (thus in the fossilization rate)
	-s <number>
		set the number of samples required
	-d
		return the distribution (otherwise the density is returned by default)
	-u <value>
		set the step discretizing the time distribution to <value>
	-n <number>
		compute the divergence distribution associated to node <number>; option can be used several times; if it is not used all the divergence times are computed; the number of each node can be found by using the option -z or the prgram 'draw' with option -d
	-x <number>
		set the graphic format of the output (option is required if one wants a graphic output)
			-f 1 -> pdf
			-f 2 -> postscript
			-f 3 -> png
			-f 4 -> svg
			-f 5 -> LaTeX (psTricks)
			-f 6 -> LaTeX (TikZ)
	-h
		display help

EXAMPLE

./tspe -o -400 -s 2000 -b 10000 -g 300 -f 0.25 -a 0.33 0.33 -w 0.1 0.1 0.01 -i 5. 5. 5. -t 40 -x 6 -n 0 -n 1 -n 2  -n 5 -n 7 -n 24 -n 25 -n 26 -n 27 -n 32 -n 67 -n 68 -n 103 -n 104 -n 129 -n 130 -n 147 -n 202 -n 203 -n 204 -n 208 -n 209 -n 236 -n 237 -n 276 CotylosauriaTree.phy CotylosauriaAges.csv distribution

--------------------------

------------
 dext 
------------

--------------------------
REQUIREMENT

	'dext' requires the gsl and pthread libraries.

--------------------------
COMPILING

	Just type
	> make dext
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'dext' computes the extinction time distibution associated to given sets of clades from a set of possible trees, the fossil ages


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	dext - computes the extinction time distibution associated to a given clade from a set of possible trees, the fossil stratigraphic intervals and provides coda files for MCMC convergence diagnostics
	
SYNOPSIS
	dext [OPTIONS] <input Tree File> <input Fossil File> [output File]

DESCRIPTION
	Compute the distribution of the divergence time associated to the clade corresponding to the list of tips given in <input List Clade> by sampling into the trees contained in <input Tree File> (in Newick format) into the fossil ranges provided in <input Fossil File> (in csv format) and output the distribution as a .csv table <output File>.csv and two files <output File>.out and <output File>.ind which can be read by the R package coda for MCMC diagnostics

	Options are
	-z <input Tree File>
		output the tree in text debug format in the console and exit 
	-o <origin bound inf> [origin bound sup]
		set the origin range; if origin bound sup is not given it is set to the most ancient fossil age
	-w <speciation rate width> <extinction rate width> <fossilization rate width>
		set the width of the sliding windows used for samplint the speciation, extinction and fossilization rates during the MCMC
	-i <speciation rate upper bound> <extinction rate  upper bound> <fossilization rate upper bound>
		set the upper bounds of the interval in which the speciation, extinction and fossilization rates are uniformly drawn at the beginning of the MCMC (lower bounds are all 0)
	-b <number>
		set the number of iterations for the burning stage
	-g <number>
		set the thinning parameter (only one iteration in <number> is considered)  
	-f <proportion>
		set the proportion of moves in the rates in the MCMC proposal (other moves are drawn uniformly amont the fossil ages)
	-a <speciation proportion> <extinction proportion>
		set the relation proportion of moves in the speciation and the extinction rates (thus in the fossilization rate)
	-s <number>
		set the number of samples required
	-r <number>
		set the random seed
	-n <taxa name>
		add the taxa to the list of taxa whose extinction distribution to be computed
	-c <file name>
		add the clade made of the taxa listed in the file to the list of clades whose extinction distribution to be computed
	-d
		return the distribution (otherwise the density is returned by default)
	-u <value>
		set the step discretizing the time distribution to <value>
	-t <number>
		set the number of thread running in parallell
	-h
		display help

EXAMPLE

./dext -o -400 -s 300 -b 10000 -g 200 -f 0.2 -a 0.33 0.33 -w 0.5 0.5 0.5 -i 5. 5. 5. -t 1 -m -230 -u 0.05 -c ../data/Ophiacodontidae -c ../data/Edaphosauridae -c ../data/Sphenacodontidae ~/Dropbox/Extinction/dev/data/Eupely100.phy ~/Dropbox/Extinction/dev/data/CotylosauriaAges.csv dist


--------------------------

-------
 text 
-------

--------------------------
REQUIREMENT

	'text' requires the gsl, pthread and cairo libraries.

--------------------------
COMPILING

	Just type
	> make text
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'text' computes the extinction time distributions associated to all the terminal nodes of a single tree with the fossil ages and draw a tree figure


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	text - Computation of the divergence time distributions associated to a given set of nodes from a single tree, the fossil ages, and the diversification rates and draw a tree figure
SYNOPSIS
	text [OPTIONS] <input Tree File> <input Fossil File> [output File]

DESCRIPTION
	Compute the divergence time distributions associated to a given set of nodes from the tree of <input Tree File>, the fossil ages of <input Fossil File>, and the diversification rates and draw a tree figure

	Options are
	-z <input Tree File>
		output the tree in text debug format in the console and exit 
	-o <origin bound inf> [origin bound sup]
		set the origin range; if origin bound sup is not given it is set to the most ancient fossil age
	-w <speciation rate width> <extinction rate width> <fossilization rate width>
		set the width of the sliding windows used for samplint the speciation, extinction and fossilization rates during the MCMC
	-i <speciation rate upper bound> <extinction rate  upper bound> <fossilization rate upper bound>
		set the upper bounds of the interval in which the speciation, extinction and fossilization rates are uniformly drawn at the beginning of the MCMC (lower bounds are all 0)
	-b <number>
		set the number of iterations for the burning stage
	-g <number>
		set the thinning parameter (only one iteration in <number> is considered)  
	-f <proportion>
		set the proportion of moves in the rates in the MCMC proposal (other moves are drawn uniformly amont the fossil ages)
	-a <speciation proportion> <extinction proportion>
		set the relation proportion of moves in the speciation and the extinction rates (thus in the fossilization rate)
	-s <number>
		set the number of samples required
	-d
		return the distribution (otherwise the density is returned by default)
	-u <value>
		set the step discretizing the time distribution to <value>
	-n <number>
		compute the divergence distribution associated to node <number>; option can be used several times; if it is not used all the divergence times are computed; the number of each node can be found by using the option -z or the prgram 'draw' with option -d
	-x <number>
		set the graphic format of the output (option is required if one wants a graphic output)
			-f 1 -> pdf
			-f 2 -> postscript
			-f 3 -> png
			-f 4 -> svg
			-f 5 -> LaTeX (psTricks)
			-f 6 -> LaTeX (TikZ)
	-h
		display help

EXAMPLE

./text -o -400 -s 100 -b 100 -g 1 -f 0.5 -a 0.33 0.33 -w 0.5 0.5 0.5 -i 5. 5. 5. -t 1 -m -230 -u 1 -y 500 -x 6 ~/Dropbox/Extinction/dev/data/Eupely10.phy ~/Dropbox/Extinction/dev/data/CotylosauriaAgesNew.csv tree

--------------------------

-------
 qext 
-------

--------------------------
REQUIREMENT

	'qext' requires the gsl, pthread and cairo libraries.

--------------------------
COMPILING

	Just type
	> make qext
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'qext' computes the quantile at a given order, of extinction time distibution associated to a given set of clades from a set of possible trees, the fossil stratigraphic intervals and provides coda files for MCMC convergence diagnostics


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	qext - Computation of the divergence time distributions associated to a given set of nodes from a single tree, the fossil ages, and the diversification rates and draw a tree figure
SYNOPSIS
	qext [OPTIONS] <input Tree File> <input Fossil File> [output File]

DESCRIPTION
	Compute the divergence time distributions associated to a given set of nodes from the tree of <input Tree File>, the fossil ages of <input Fossil File>, and the diversification rates and draw a tree figure

	Options are
	-z <input Tree File>
		output the tree in text debug format in the console and exit 
	-o <origin bound inf> [origin bound sup]
		set the origin range; if origin bound sup is not given it is set to the most ancient fossil age
	-w <speciation rate width> <extinction rate width> <fossilization rate width>
		set the width of the sliding windows used for samplint the speciation, extinction and fossilization rates during the MCMC
	-i <speciation rate upper bound> <extinction rate  upper bound> <fossilization rate upper bound>
		set the upper bounds of the interval in which the speciation, extinction and fossilization rates are uniformly drawn at the beginning of the MCMC (lower bounds are all 0)
	-b <number>
		set the number of iterations for the burning stage
	-g <number>
		set the thinning parameter (only one iteration in <number> is considered)  
	-f <proportion>
		set the proportion of moves in the rates in the MCMC proposal (other moves are drawn uniformly amont the fossil ages)
	-a <speciation proportion> <extinction proportion>
		set the relation proportion of moves in the speciation and the extinction rates (thus in the fossilization rate)
	-s <number>
		set the number of samples required
	-d
		return the distribution (otherwise the density is returned by default)
	-u <value>
		set the step discretizing the time distribution to <value>
	-n <number>
		compute the divergence distribution associated to node <number>; option can be used several times; if it is not used all the divergence times are computed; the number of each node can be found by using the option -z or the prgram 'draw' with option -d
	-h
		display help

EXAMPLE

./cext -o -400 -s 1000 -b 10000 -g 200 -f 0.2 -a 0.33 0.33 -w 0.5 0.5 0.5 -i 5. 5. 5. -t 1 -m -230 -u 0.05 ~/Dropbox/Extinction/dev/data/Sphenacodontidae ~/Dropbox/Extinction/dev/data/Ophiacodontidae ~/Dropbox/Extinction/dev/data/Eupely.phy ~/Dropbox/Extinction/dev/data/CotylosauriaAges.csv

--------------------------

-------
 cext 
-------

--------------------------
REQUIREMENT

	'cext' requires the gsl, pthread and cairo libraries.

--------------------------
COMPILING

	Just type
	> make cext
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'cext' computes the probability that a set of extinct taxa goes extinct before another one from a set of possible trees, the fossil stratigraphic intervals and provides coda files for MCMC convergence diagnostics


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	cext - Computation of the divergence time distributions associated to a given set of nodes from a single tree, the fossil ages, and the diversification rates and draw a tree figure
SYNOPSIS
	cext [OPTIONS] <input Tree File> <input Fossil File> [output File]

DESCRIPTION
	Compute the divergence time distributions associated to a given set of nodes from the tree of <input Tree File>, the fossil ages of <input Fossil File>, and the diversification rates and draw a tree figure

	Options are
	-z <input Tree File>
		output the tree in text debug format in the console and exit 
	-o <origin bound inf> [origin bound sup]
		set the origin range; if origin bound sup is not given it is set to the most ancient fossil age
	-w <speciation rate width> <extinction rate width> <fossilization rate width>
		set the width of the sliding windows used for samplint the speciation, extinction and fossilization rates during the MCMC
	-i <speciation rate upper bound> <extinction rate  upper bound> <fossilization rate upper bound>
		set the upper bounds of the interval in which the speciation, extinction and fossilization rates are uniformly drawn at the beginning of the MCMC (lower bounds are all 0)
	-b <number>
		set the number of iterations for the burning stage
	-g <number>
		set the thinning parameter (only one iteration in <number> is considered)  
	-f <proportion>
		set the proportion of moves in the rates in the MCMC proposal (other moves are drawn uniformly amont the fossil ages)
	-a <speciation proportion> <extinction proportion>
		set the relation proportion of moves in the speciation and the extinction rates (thus in the fossilization rate)
	-s <number>
		set the number of samples required
	-d
		return the distribution (otherwise the density is returned by default)
	-u <value>
		set the step discretizing the time distribution to <value>
	-n <number>
		compute the divergence distribution associated to node <number>; option can be used several times; if it is not used all the divergence times are computed; the number of each node can be found by using the option -z or the prgram 'draw' with option -d
	-h
		display help

EXAMPLE

./cext -o -400 -s 1000 -b 10000 -g 200 -f 0.2 -a 0.33 0.33 -w 0.5 0.5 0.5 -i 5. 5. 5. -t 1 -m -230 -u 0.05 ~/Dropbox/Extinction/dev/data/Sphenacodontidae ~/Dropbox/Extinction/dev/data/Ophiacodontidae ~/Dropbox/Extinction/dev/data/Eupely.phy ~/Dropbox/Extinction/dev/data/CotylosauriaAges.csv

--------------------------

------------
 diag 
------------

--------------------------
REQUIREMENT

	'diag' requires the gsl and pthread libraries.

--------------------------
COMPILING

	Just type
	> make diag
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'diag' provides coda files for MCMC convergence diagnostics from a set of possible trees, the fossil stratigraphic intervals


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	diag - provides coda files for MCMC convergence diagnostics from a set of possible trees, the fossil stratigraphic intervals
	
SYNOPSIS
	diag [OPTIONS] <input Tree File> <input Fossil File> [output File]

DESCRIPTION
	Provide coda files for MCMC convergence diagnostics by sampling into the trees contained in <input Tree File> (in Newick format) into the fossil ranges provided in <input Fossil File> (in csv format). Output two files <output File>.out and <output File>.ind which can be read by the R package coda for MCMC diagnostics

	Options are
	-z <input Tree File>
		output the tree in text debug format in the console and exit 
	-o <origin bound inf> [origin bound sup]
		set the origin range; if origin bound sup is not given it is set to the most ancient fossil age
	-w <speciation rate width> <extinction rate width> <fossilization rate width>
		set the width of the sliding windows used for samplint the speciation, extinction and fossilization rates during the MCMC
	-i <speciation rate upper bound> <extinction rate  upper bound> <fossilization rate upper bound>
		set the upper bounds of the interval in which the speciation, extinction and fossilization rates are uniformly drawn at the beginning of the MCMC (lower bounds are all 0)
	-b <number>
		set the number of iterations for the burning stage
	-g <number>
		set the thinning parameter (only one iteration in <number> is considered)  
	-f <proportion>
		set the proportion of moves in the rates in the MCMC proposal (other moves are drawn uniformly amont the fossil ages)
	-a <speciation proportion> <extinction proportion>
		set the relation proportion of moves in the speciation and the extinction rates (thus in the fossilization rate)
	-s <number>
		set the number of samples required
	-r <number>
		set the random seed
	-t <number>
		set the number of thread running in parallell
	-h
		display help

EXAMPLE

./diag -o -400 -s 2000 -b 10000 -g 300 -r 5555 -f 0.25 -a 0.33 0.33 -w 0.1 0.1 0.01 -i 5. 5. 5. -t 40 CotylosauriaTree.phy CotylosauriaAges.csv amniota2000_400

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
	-d
		print the numbers associated to the nodes
	-x <number>
		set the graphic format of the output (option is required if one wants a graphic output)
			-f 1 -> pdf
			-f 2 -> postscript
			-f 3 -> png
			-f 4 -> svg
			-f 5 -> LaTeX (psTricks)
			-f 6 -> LaTeX (TikZ)
	-h
		display help

--------------------------

-------
 samp 
-------

--------------------------
REQUIREMENT

	'samp' requires the cairo library.

--------------------------
COMPILING

	Just type
	> make samp
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'samp' draws a single tree with the fossil ages


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	samp - Drawing a single tree with the fossil ages
	
SYNOPSIS
	samp [OPTIONS] <input Tree File> <input Fossil File>

DESCRIPTION
	Output a figure in various graphic formats of the tree in <input Tree File> with the fossil ages of  <input Fossil File>

	Options are
	-z <input Tree File>
		output the tree in text debug format in the console and exit 
	-o <origin> 
		set the origin time
	-e <end> 
		set the end time
	-d
		print the numbers associated to the nodes
	-x <number>
		set the graphic format of the output (option is required if one wants a graphic output)
			-f 1 -> pdf
			-f 2 -> postscript
			-f 3 -> png
			-f 4 -> svg
			-f 5 -> LaTeX (psTricks)
			-f 6 -> LaTeX (TikZ)
	-h
		display help

--------------------------
