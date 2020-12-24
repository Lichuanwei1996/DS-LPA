To run the program, under MS-DOS, type:



benchmark [FLAG] [P]


[FLAG] [P]



-N		number of nodes

-k		average degree

-maxk		maximum degree

-mu		mixing parameter

-t1		minus exponent for the degree sequence

-t2		minus exponent for the community size distribution

-minc		minimum for the community sizes

-maxc		maximum for the community sizes

-on		number of overlapping nodes

-om		number of memberships of the overlapping nodes

-C              average clustering coefficient

Example1:

 benchmark -N 1000 -k 15 -maxk 50 -mu 0.1 -minc 20 -maxc 50

Example2:

 benchmark -f flags.dat -t1 3


If you want to produce a kind of Girvan-Newman benchmark, you can type:

 benchmark -N 128 -k 16 -maxk 16 -mu 0.1 -minc 32 -maxc 32

Output:
1) network.dat contains the list of edges (nodes are labelled from 1 to the number of nodes; the edges are ordered and repeated twice, i.e. source-target and target-source).


2) community.dat contains a list of the nodes and their membership (memberships are labelled by integer numbers >=1).


3) statistics.dat contains the degree distribution (in logarithmic bins), the community size distribution, and the distribution of the mixing parameter