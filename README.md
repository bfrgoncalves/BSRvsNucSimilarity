# BSRvsNucSimilarity

#####This program plots the values of Blast Score Ratio and nucleotide sequence similarity in an interval of n generations, under a specific mutation rate.

python BSRvsNucSimilarity.py -i <Inputfile> -n <Generations> -c <Ncalc> -m <MutationRate>

-i "Input query file with nucleotide sequences"

-n 'Number of generations passed between each calculation of BSR and nucleotide similarity.'

-c 'Number of BSR and nucleotide similarity calculations.'

-m "Mutation Rate"

Example:

python BSRvsNucSimilarity.py -i query.fasta -n 50 -c 50 -m 0.06