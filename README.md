# BSRvsNucSimilarity

#####This program plots the values of Blast Score Ratio and nucleotide sequence similarity in an interval of n generations, under a specific mutation rate.

python BSRvsNucSimilarity.py -i <Inputfile> -n <Generations> -m <MutationRate>

-i "Input query file with nucleotide sequences"

-n 'Number of generations passed and calculations of BSR and nucleotide similarity.'

-m "Mutation Rate"

Example:

python BSRvsNucSimilarity.py -i query.fasta -n 20 -m 0.06