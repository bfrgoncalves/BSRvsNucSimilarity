from Bio import SeqIO
from BCBio import GFF
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
import sys
import os
import re
import HTSeq
import subprocess
from Bio.Seq import Seq
from os import listdir
from os.path import isfile, join
import time
from CommonFastaFunctions import runBlastParser
import random
import matplotlib.pyplot as plt
import argparse


#Criar a Blast DB:
def Create_Blastdb( questionDB, overwrite, dbtypeProt,dbName ):
    isProt=dbtypeProt

    if not os.path.isfile(dbName + ".nin") and not os.path.isfile(dbName + ".nhr") and not os.path.isfile(dbName + ".nsq"):
        
        if not isProt:
            os.system( "makeblastdb -in " + questionDB + " -out " + dbName + " -dbtype nucl -logfile " + dbName + "_blast.log" )
        else:
            os.system( "makeblastdb -in " + questionDB + " -out " + dbName + " -dbtype prot -logfile " + dbName + "_blast.log" )

    elif overwrite:
        if not isProt:
            os.system( "makeblastdb -in " + questionDB + " -out " + dbName + " -dbtype nucl -logfile " + dbName + "_blast.log" )
        else:
            os.system( "makeblastdb -in " + questionDB + " -out " + dbName + " -dbtype prot -logfile " + dbName + "_blast.log" )

    else:
        print "BLAST DB files found. Using existing DBs.."  
    return( dbName )


def Create_FASTAquery(geneName,geneSequence):
        of_handle=open(nameFASTA, 'w')
        os.chmod(nameFASTA, 0755)
        of_handle.write(">gn|"  + geneName + "|\n" + geneSequence +"\n")
        of_handle.close()

def Align_sort_key(Result):
    #x=len(Result.split(";")[0].split("--"))
    AlignPos=Result.split(";")[0].split("--")[0]
    return int(AlignPos)
####################################################################################################################################


def reverseComplement(strDNA):

    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    strDNArevC = ''
    for l in strDNA:

        strDNArevC += basecomplement[l]

    return strDNArevC[::-1]

def translateSeq(DNASeq):
    seq=DNASeq
    try:
        #print seq
        myseq= Seq(seq)
        #print "aqui"
        #print myseq
        protseq=Seq.translate(myseq, table=11,cds=True)
        #print protseq

    except:
        try:
            seq=reverseComplement(seq)
            myseq= Seq(seq)
            #print myseq
            protseq=Seq.translate(myseq, table=11,cds=True)
                        
        except:
            try:
                seq=seq[::-1]
                myseq= Seq(seq)
                #print myseq
                protseq=Seq.translate(myseq, table=11,cds=True)
            except:
                try:
                    seq=seq[::-1]                           
                    seq=reverseComplement(seq)
                    myseq= Seq(seq)
                    #print myseq
                    protseq=Seq.translate(myseq, table=11,cds=True)
                except:
                    raise
    return protseq



def getOwnBlastScore(FASTAfile):
    gene_fp = HTSeq.FastaReader(FASTAfile)
    #alleleI=0
    names=""
    alleleProt=''
    proteome=""
    bestmatch=0
    for allele in gene_fp: #new db for each allele to blast it against himself
        #print allele.seq
        #print type(allele.seq)
        x = str(translateSeq(allele.seq))
        #print str(allele.name)
        #names=allele.name.split("|")[3]
        #print allele.seq
        alleleProt+=">"+str(allele.name)+"\n"+x+"\n"
        proteome+=">"+str(allele.name)+"\n"+x+"\n"
        #print alleleProt
    with open('referenceAA.fasta', "wb") as f:
        f.write(alleleProt)
    with open('proteome.fasta', "wb") as v:
        v.write(proteome)
    name = "Databases/reference_db"
    Gene_Blast_DB_name = Create_Blastdb('referenceAA.fasta',1,True, name)
        # --- get BLAST score ratio --- #
    cline = NcbiblastpCommandline(query='proteome.fasta', db=name, out=blast_out_file, outfmt=5)
        #print cline
    allelescore=0
    blast_records = runBlastParser(cline,blast_out_file, "")
    allelescores={}
    for blast_record in blast_records:
        found=False 
        for alignment in blast_record.alignments:
            if found is False:
                #print blast_record.query, alignment.hit_def
                for match in alignment.hsps:
                    #print alignment.hit_def
                        #print "---------------------"
                    #print alignment.hit_def
                    #print blast_record.query
                    #print alignment.hit_def
                    try:
                        if allelescores[str(alignment.hit_def)] < match.score:
                            allelescores[str(alignment.hit_def)] = int(match.score)
                            break
                    except KeyError:
                        allelescores[str(alignment.hit_def)] = int(match.score)
                        break
                bestmatch=allelescores[str(alignment.hit_def)]
            else:
                break
    #print allelescores
    #for i in allelescores:
        #hitsName.append(str(i)+";"+str(allelescores[i])+";")
    #hitsName.sort(key=Align_sort_key)
    #print hitsName
    #return alleleI,allelescores,Gene_Blast_DB_name
    #print alleleI
    #print len(allelescores)
    print bestmatch
    return bestmatch



def getBlastScoreRatios(pathQuery,allelescores):
    
    g_fp = HTSeq.FastaReader(pathQuery)
    countGenes=0
    proteome=""
    for contig in g_fp:
        #IdCDS=str(j)+"--"+str(countCDS)
        #print type(contig.seq)
        protseq=translateSeq(contig.seq)
        proteome+=">"+str(contig.name)+"\n"+str(protseq)+"\n"
        #print proteome
        if os.path.isfile('proteome.fasta'):
            os.remove('proteome.fasta')
        with open('proteome.fasta', "wb") as v:
            v.write(proteome)

    name = "Databases/reference_db"
    cline = NcbiblastpCommandline(query='proteome.fasta', db=name, out=blast_out_file, outfmt=5)
        #print cline

    allelescore=0
    blast_records = runBlastParser(cline,blast_out_file, "")

    #os.remove('proteome.fasta')
    
    blastScoreRatio=0
    countT=0
    matchR=[]
    m=[]
    prevName=""
    bestmatch=""
    #bestmatches={}
    for blast_record in blast_records:
        found=False 
        for alignment in blast_record.alignments:
            if found is False:
                #print blast_record.query, alignment.hit_def
                scoreToUse=0
                for match in alignment.hsps:
                    if float(scoreToUse)<float(match.score):
                        scoreToUse=str(match.score)
                    #print alignment.hit_def
                        #print "---------------------"
                    #print alignment.hit_def
                    #print blast_record.query
                    #print alignment.hit_def
                    blastScoreRatio = float(match.score) / float(allelescores)
                    #print blastScoreRatio

                    cdsStrName=blast_record.query
                    if(blastScoreRatio == 1 and bestmatches[str(alignment.hit_def)][2]=="No"):
                        bestmatches[str(alignment.hit_def)]=[str(match.score),str(blastScoreRatio),"Yes"]
                        #print alignment
                        #print match
                    elif(blastScoreRatio == 1 and match.score>float(bestmatches[str(alignment.hit_def)][0])):
                        bestmatches[str(alignment.hit_def)]=[str(match.score),str(blastScoreRatio),"Yes"]
                        #print match
                    elif(match.score>float(bestmatches[str(alignment.hit_def)][0]) and blastScoreRatio>0.1 and blastScoreRatio>float(bestmatches[str(alignment.hit_def)][1])):
                        #print match.query
                        #print match.sbjct
                        #print allelescores
                        bestmatches[str(alignment.hit_def)]=[str(match.score),str(blastScoreRatio),"Yes"]

                bestmatch=bestmatches[str(alignment.hit_def)][1]
                print 

            else:
                break


    return bestmatch

def getNucSimilarity(pathQuery, IsCreateDB, BestScore, IsFirstTime):
    DbName = "Databases/NucSimi_db"
    if IsCreateDB == "yes":
        Gene_Blast_DB_name = Create_Blastdb('reference.fasta',1,True,DbName)
        # --- get BLAST score ratio --- #
    cline = NcbiblastpCommandline(query=pathQuery, db=DbName, out=blast_out_file, outfmt=5, evalue=0.001)
        #print cline
    allelescore=0
    blast_records = runBlastParser(cline,blast_out_file, "")
    nucScores={}
    #print blast_records
    for blast_record in blast_records:
        found=False 
        for alignment in blast_record.alignments:
            if found is False:
                #print blast_record.query, alignment.hit_def
                for match in alignment.hsps:
                    #print alignment.hit_def
                        #print "---------------------"
                    #print alignment.hit_def
                    #print blast_record.query
                    #print alignment.hit_def
                    print 
                    try:
                        if nucScores[str(alignment.hit_def)] < match.score:
                            nucScores[str(alignment.hit_def)] = float(match.score)
                            #print nucScores
                            #print str(alignment.hit_def)
                            break
                    except KeyError:
                        nucScores[str(alignment.hit_def)] = float(match.score)
                        break
            else:
                break

    if IsFirstTime:
        BestScore = nucScores[str(1)]

    return BestScore, (nucScores[str(1)]/BestScore)

def checkSimilarity(querySeq,refSeq):
    countmatches=0
    for i in range(0,len(querySeq)):
        if querySeq[i] == refSeq[i]:
            countmatches+=1

    return float(float(countmatches)/float(len(refSeq)))

def mutate(orig_string, mutation_rate=0.00066):
    bases="ACTG"
    result = []
    mutations = []
    for base in orig_string:
        if random.random() < mutation_rate:
            new_base = bases[bases.index(base) - random.randint(1, 3)] # negatives are OK
            result.append(new_base)
            mutations.append((base, new_base))
        else:
            result.append(base)
    #return "".join(result), mutations
    return "".join(result)

def checkStops(sequence):
    sequence=sequence
    bases="ACG"
    stopCodons = ["TAG","TAA","TGA"]
    newCodons=[]
    codons = [sequence[i:i+3] for i in range(0,len(sequence),3)]
    countCodons=0
    for j in codons:
        for i in stopCodons:
            if i == j:
                return True
            #print x
    return False

def getFASTAfeatures(pathF):
    g_fp = HTSeq.FastaReader(pathF)
    countGenes=0
    #print "aqui"
    for contig in g_fp:
        countGenes+=1
        bestmatches[str(contig.name)] = ["0","0","No"] # aleleID : score, score ratio, perfectmatch, key name of the DNA sequence string, Isfound, "NewAllele"

####################################################################################################################################################


#python database_search.py isBSR

def main():
    
    
    
    parser = argparse.ArgumentParser(description="This program plots the values of Blast Score Ratio and nucleotide sequence similarity in an interval of n generations, under a specific mutation rate.")
    parser.add_argument('-i', nargs='?', type=str, help="Input fasta file with nucleotide sequences", required=True)
    parser.add_argument('-n', nargs='?', type=int, help='Number of generations passed between each calculation of BSR and nucleotide similarity.', required=True)
    parser.add_argument('-c', nargs='?', type=int, help='Number of BSR and nucleotide similarity calculations.', required=True)
    parser.add_argument('-m', nargs='?', type=float, help="Mutation Rate", required=False)

    args = parser.parse_args()

    mutRate=args.m
    inputFileQ = args.i
    gen=args.n
    nPlot=args.c

    pathRef= "reference_seq123.fasta"
    pathQuery = "query_seq123.fasta"

    BSR = []
    NuScore= []
    generations = []
    print 'Running ' + str(gen*nPlot) +" generations."

    inputF = HTSeq.FastaReader(inputFileQ)

    for contig in inputF:
        print 'Running ' + str(gen*nPlot) +" generations for sequence " + contig.name + "."
        ref = open(pathRef,"w")
        ref.write(">"+contig.name+"\n"+contig.seq+"\n")
        ref.close()
        query = open(pathQuery,"w")
        query.write(">"+contig.name+"\n"+contig.seq+"\n")
        query.close()

        allelescores={}
        bestmatches={}
        getFASTAfeatures(pathRef)
        refScore = getOwnBlastScore(pathRef)
        bestBSR = getBlastScoreRatios(pathQuery,refScore)
        #print bestBSR
        queryFile = HTSeq.FastaReader(pathQuery)
        for contig in queryFile:
            name=contig.name
            refSeq = contig.seq

        Identity = checkSimilarity(refSeq,refSeq)
        #bestNuScore, Identity = getNucSimilarity("query.fasta","yes",0,True)
        Iden=Identity
        BSR.append(float(bestBSR))
        NuScore.append(float(Identity))
        generations.append([float(bestBSR),float(Identity)])

        for i in range(0,nPlot):
            newSeq = ""
            queryFile = HTSeq.FastaReader(pathQuery)
            for contig in queryFile:
                name=contig.name
                for j in range(0,gen):
                    start = contig.seq[0:3]
                    #print start
                    end = contig.seq[len(contig.seq)-3:]
                    seqToUse = contig.seq[3:len(contig.seq)-3]
                    #print seqToUse
                    newSeq = mutate(seqToUse,mutRate)
                    IsStop = checkStops(newSeq)
                    if IsStop:
                        newSeq = contig.seq
                    else:
                        newSeq = start + newSeq + end
                newfile = open(pathQuery,"w")
                newfile.write(">"+str(name)+"\n"+newSeq+"\n")
                newfile.close()

            bestmatches={}
            allelescores={}
            getFASTAfeatures(pathRef)
            refScore = getOwnBlastScore(pathRef)
            bestBSR = getBlastScoreRatios(pathQuery,refScore)
            Identity = checkSimilarity(newSeq,refSeq)
            #bestNuScore, Identity = getNucSimilarity("query.fasta","no", bestNuScore, False)
            BSR.append(float(bestBSR))
            NuScore.append(float(Identity))

        generations.append([BSR,NuScore])


        os.remove("proteome.fasta")
        os.remove("referenceAA.fasta")
        os.remove(pathQuery)
        os.remove(pathRef)

    print generations
    plt.plot(BSR,NuScore, 'ro')
    plt.axis([0, 1, 0, 1])
    plt.xlabel('BSR')
    plt.ylabel('Nuc Similarity')
    plt.show()

    print "DONE"

if __name__ == "__main__":
    blast_out_file = "Databases/results_subj.xml"
    allelescores={}
    bestmatches={}
    main()

