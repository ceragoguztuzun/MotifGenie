import requests
from bs4 import BeautifulSoup
import numpy  as np 
import pandas as pd 
import argparse
import sys
import matplotlib.pyplot as plt
plt.ion()
import logomaker as lm
import math


jasparfileName = None
topN = None
thresholdRatioForConsensus = 0.66
coreMotifSearchFlag = False
seqLogoFlag = False
cellLine = None
genome_selected = "hg38"
upperScoreThr = 9999
lowerScoreThr = 0

def getData(genome_selected, chr_no, start_pos, end_pos):
    url = "api.genome.ucsc.edu/getData/sequence?genome={genome_name};chrom={chr};start={start};end={end}".format(genome_name=genome_selected, chr=chr_no, start=start_pos, end=end_pos)
    resp = requests.get('https://'+url)
    if resp.status_code != 200:
        # This means something went wrong.
        if resp.json()['error']:
            print("ERROR:",resp.json()['error'])
        else:
            print("ERROR while connecting")
        sys.exit(1)

    dna_seq = resp.json()['dna'].upper()
    return dna_seq

def getChrData(genome_selected, chr_no):
    url = "api.genome.ucsc.edu/getData/sequence?genome={genome_name};chrom={chr}".format(genome_name=genome_selected, chr=chr_no)
    resp = requests.get('https://'+url)
    if resp.status_code != 200:
        # This means something went wrong.
        if resp.json()['error']:
            print("ERROR:",resp.json()['error'])
        else:
            print("ERROR while connecting")
        sys.exit(1)
    dna_seq = resp.json()['dna'].upper()
    start_pos = int(resp.json()['start'])
    end_pos = int(resp.json()['end'])
    return dna_seq, start_pos, end_pos

# converts positive DNA strand to negative DNA strand and vica versa
def convertDNAStrand(seq):
    seq_negative = seq[::-1]
    seq_negative = seq_negative.lower()
    seq_negative = seq_negative.replace("a","T")
    seq_negative = seq_negative.replace("t","A")
    seq_negative = seq_negative.replace("g","C")
    seq_negative = seq_negative.replace("c","G")
    return seq_negative


# this method displays the consensus of the inputted motif_df. 
# if ratio of position is below threshold_ratio_for_n it is displayed as "n"
# threshold_ratio_for_consensus is 0.66 which is referenced by ESR1's pfm
# consensus gets trimmed if threshold_ratio_for_consensus is not equal to 0

def getConsensus( motif_df, threshold_ratio_for_consensus = 0.66, threshold_ratio_for_n = 0.6):
    
    ratios_list = [ motif_df.loc[i].max()/motif_df.loc[i].sum() for i in range( motif_df.shape[0])]
    indices = []
    indices = [m for m in range(len(ratios_list)) if ratios_list[m] >= threshold_ratio_for_consensus]
    
    if threshold_ratio_for_consensus != 0:
        motif_df = motif_df[ indices[0]: indices[-1]+1]
        ratios_list = ratios_list[ indices[0]: indices[-1]+1]
    
    consensus_str = [ motif_df.loc[ indices[0] + m].idxmax() if ratios_list[m] >= threshold_ratio_for_n else "n" for  m in range(len(ratios_list))]
    consensus_str = ''.join( consensus_str)

    return motif_df, consensus_str

# method to generate and display sequence logo for the motifs of a TF
'''
    W(b,i) = log2(P(b,i)/P(b))
    P(b,i) = frequency from table
    P(b) = expected frequency of nucleotides (taken as 0.25)
    
    The sequence score gives an indication of how different the sequence is from a random sequence.
    e. The score is 0 if the sequence has the same probability of being a functional site and of being a random site. 
    The score is greater than 0 if it is more likely to be a functional site than a random site, 
    and less than 0 if it is more likely to be a random site than a functional site.
'''
def displaySeqLogo( list_of_motifs, search_type, jasparfileName):

    store_frequencies = []
    motif_length = len(min(list_of_motifs, key = len))
    # iterate for each position considering the minimum motif length -------------- (might modify minimum part later) 
    for column_of_motif in range(motif_length):
        (motif, count) = np.unique([m[column_of_motif] for m in list_of_motifs], return_counts=True)
        frequencies = np.asarray( (motif, count) )
        frequencies_percentage = [float(freq)/len([m[column_of_motif] for m in list_of_motifs]) for freq in frequencies[1]]
    
        # organizing data
        zeros = np.zeros(4)
        for i in range(len(frequencies[0])):
            nucleotide = frequencies[0][i]
            if nucleotide == "A":
                zeros[0] = frequencies_percentage[i]
            elif nucleotide == "C":
                zeros[1] = frequencies_percentage[i]
            elif nucleotide == "G":
                zeros[2] = frequencies_percentage[i]
            elif nucleotide == "T":
                zeros[3] = frequencies_percentage[i]
            
        store_frequencies.append(zeros)

    store_frequencies = np.asarray(store_frequencies)

    for x in range(store_frequencies.shape[0]):
        hu = (store_frequencies[x][0]+0.01)*np.log2(store_frequencies[x][0]+0.01)+(store_frequencies[x][1]+0.01)*np.log2(store_frequencies[x][1]+0.01)+(store_frequencies[x][2]+0.01)*np.log2(store_frequencies[x][2]+0.01)+(store_frequencies[x][3]+0.01)*np.log2(store_frequencies[x][3]+0.01)
        store_frequencies[x][0] = store_frequencies[x][0]*(2.0+hu)
        store_frequencies[x][1] = store_frequencies[x][1]*(2.0+hu)
        store_frequencies[x][2] = store_frequencies[x][2]*(2.0+hu)
        store_frequencies[x][3] = store_frequencies[x][3]*(2.0+hu)

    table_of_frequencies = pd.DataFrame(store_frequencies, columns=['A', 'C', 'G', 'T']) # ---> table of frequencies obtained
    
    name = jasparfileName + "-" + search_type
    logo = lm.Logo(table_of_frequencies, stack_order = "big_on_top")
    logo.ax.set_xlabel('Position\n('+name+")",fontsize=10)
    logo.ax.set_ylabel("Entropy (bits)", labelpad=-1,fontsize=10)
    #logo.savefig(name + '.png', bbox_inches='tight')
    logo.ax.figure.canvas.start_event_loop(sys.float_info.min) #workaround for Exception in Tkinter callback
    logo.ax.figure.savefig(name+".png", dpi=300, bbox_inches='tight') 
    
    #logo.ax.figure.show()


# searches motifs from JASPAR in inputted LOCUS:
def searchinSeq(sequence, conversionFlag, startIndex, coreMotifSearchFlag, scores_list, antiparallelFlag, peak_name, entry_name = "", topN=15, fileSelection=0): 
    
    # read pfm file from JASPAR
    file_name = open(fileSelection,"r")
    contents = file_name.read()
    
    # form a matrix from read file
    contents = contents.split("\n")
    tf_name = contents[0].split(" ")[1]

    motif_df = {"A":[],
                "C":[],
                "G":[],
                "T":[]}
    keys = ["A","C","G","T"]

    for i in range(4):
        position_line = contents[i+1].split(" ")
        position_line = [float(x) for x in position_line if x != ""]
        motif_df[keys[i]] = position_line


    motif_df = pd.DataFrame(motif_df)
    # ...
    
    if antiparallelFlag == 1:
        motif_df.columns = ['T', 'G', 'C', 'A']
        motif_df = motif_df.iloc[::-1].reset_index(drop=True)
    
    # form score matrix for profile scoring
    match_df = np.zeros((motif_df.shape[0], 4))

    # core motif
    if (thresholdRatioForConsensus):
        thr_consensus = thresholdRatioForConsensus
    else:
        thr_consensus = 0.66
    thr_n = 0.6
    
    if coreMotifSearchFlag == 0:
        motif_df, consensus = getConsensus( motif_df, 0, thr_n)
    elif coreMotifSearchFlag == 1:
        motif_df, consensus = getConsensus( motif_df, thr_consensus, thr_n)

    for i in range( motif_df.shape[0]):
        sum = 0
        for j in range(4):
            sum += motif_df.iat[i,j]
        for j in range(4):
            motif_df.iat[i,j] /= (sum+1)
        total = np.log( 1/ (sum+1))

        match_df[i,0] = 3*( np.log( 1-motif_df.iat[i,0])/ total)-2*( np.log(1-motif_df.iat[i,1]) / total + np.log( 1- motif_df.iat[i,2]) / total + np.log( 1-motif_df.iat[i,3])/total )

        match_df[i,1] = 3*( np.log( 1-motif_df.iat[i,1])/ total)-2*( np.log(1-motif_df.iat[i,0]) / total + np.log( 1- motif_df.iat[i,2]) / total + np.log( 1-motif_df.iat[i,3])/total )
        
        match_df[i,2] = 3*( np.log( 1-motif_df.iat[i,2])/ total)-2*( np.log(1-motif_df.iat[i,1]) / total + np.log( 1- motif_df.iat[i,0]) / total + np.log( 1-motif_df.iat[i,3])/total )
        
        match_df[i,3] = 3*( np.log( 1-motif_df.iat[i,3])/ total)-2*( np.log(1-motif_df.iat[i,1]) / total + np.log( 1- motif_df.iat[i,2]) / total + np.log( 1-motif_df.iat[i,0])/total )
    
    # calculate scores
    i = 0
    while i <= (len(sequence) - motif_df.shape[0]):
        subseq = sequence[i:i+motif_df.shape[0]]
        
        score = 0
        if antiparallelFlag == 0:
            for m in range(len(subseq)):
                if subseq[m] == "A":
                    score += match_df[m,0]
                elif subseq[m] == "C":
                    score += match_df[m,1]
                elif subseq[m] == "G":
                    score += match_df[m,2]
                elif subseq[m] == "T":
                    score += match_df[m,3]
                else:
                    score -= 3
        elif antiparallelFlag == 1:
            for m in range(len(subseq)):
                if subseq[m] == "T":
                    score += match_df[m,0]
                elif subseq[m] == "G":
                    score += match_df[m,1]
                elif subseq[m] == "C":
                    score += match_df[m,2]
                elif subseq[m] == "A":
                    score += match_df[m,3]
                else:
                    score -= 3
                    
        if conversionFlag == 1:
            pos = len(sequence) - i - len(subseq)
        elif conversionFlag == 0:
            pos = i
        
        if antiparallelFlag == 1:
            # Antiparallel Search
            stat = "-" 
        else:
            # Parallel Search 
            stat = "+"

        antiparallelConsensus = convertDNAStrand(consensus)
        startIndex = int(startIndex)
        scores_list.append( (subseq, score, pos + startIndex, pos + startIndex + len(subseq) - 1, stat, thr_n, thr_consensus, peak_name, entry_name))
        
        i+=1
    return scores_list

def outputResult(scores_list, coreMotifSearchFlag, search_type, topN):
    
    list_of_motifs =  []
    dash = '-' * 122
    global upperScoreThr
    global lowerScoreThr

    # Display general header
    print(dash)
    print(search_type)
    #print("Threshold for n =", scores_list[0][5])
    if coreMotifSearchFlag == True:
        print("Threshold to trim consensus =", scores_list[0][6])

    # Display output
    scores_list = sorted(scores_list, key = lambda x: x[1], reverse = True)

    if upperScoreThr != None:
        upperScoreThr = float(upperScoreThr)
    if lowerScoreThr != None:
        lowerScoreThr = float(lowerScoreThr)

    if upperScoreThr != None and lowerScoreThr != None: 
        if upperScoreThr < lowerScoreThr:
            print("Upper Score Threshold cannot be less than the Lower Score Threshold.")
            sys.exit(0)
        scores_list = list(filter(lambda x:float(x[1]) >= lowerScoreThr and upperScoreThr >= float(x[1]), scores_list))
    
    elif upperScoreThr != None:
        scores_list = list(filter(lambda x: float(x[1]) <= upperScoreThr , scores_list))
        
    elif lowerScoreThr != None: 
        scores_list = list(filter(lambda x:float(x[1]) >= lowerScoreThr, scores_list))

    result = scores_list[:topN]
    
    print(dash)
    if result[0][8] == "":
        print('{:<30s}{:<20s}{:<30s}{:<15s}{:<15s}'.format("SEQUENCE","SCORE","POSITION","PARALLELITY","PEAK")) 
    else:
        print('{:<30s}{:<20s}{:<30s}{:<15s}{:<15s}{:<30s}'.format("SEQUENCE","SCORE","POSITION","PARALLELITY","PEAK","ENTRY")) 
    print(dash)
    
    for i in range(len(result)):
        if result[i][8] == "":
            print('{:<30s}{:<20f}{:<30s}{:<15s}{:<15s}'.format(result[i][0], result[i][1], str(result[i][2]) + " - " + str(result[i][3]), result[i][4],result[i][7]))
        else:
            print('{:<30s}{:<20f}{:<30s}{:<15s}{:<15s}{:<30s}'.format(result[i][0], result[i][1], str(result[i][2]) + " - " + str(result[i][3]), result[i][4],result[i][7],result[i][8]))
            
        list_of_motifs.append(result[i][0])
        
    # Display seqlogo
    if seqLogoFlag == True:
        displaySeqLogo( list_of_motifs, search_type, jasparfileName)
    print("\n")

def searchCallHelper( geneName, conversionFlag, scores_list, peaks, topN=15, fileSelection=0):
    
    if geneName != None:
        chr_dict = displayChrMenu(genome_selected,0)
        
    
    for peak_interval in peaks:
        chr_no = peak_interval[0]
        startLocation = int(peak_interval[1])
        endLocation = int(peak_interval[2])
        peak_name = peak_interval[3]
        peak_prc = peak_interval[4]

        
        # gene specific search
        if geneName != None:
            ex_list = chr_dict[chr_no]
            for entry in ex_list:
                if entry['name2'] == geneName and entry['name'][:2] == "NM":
                    if startLocation >= entry['txStart'] and endLocation <= entry['txEnd']:
                    
                        entry_name = entry['name']
                        seq = getData(genome_selected, chr_no, startLocation, endLocation)
                        seq = convertDNAStrand(seq) #seq_negative
                        start_index = startLocation
                    
                        scores_list = searchinSeq( seq, conversionFlag, start_index, coreMotifSearchFlag, scores_list, 0, peak_name, entry_name, topN, fileSelection) #parallel search
                        scores_list = searchinSeq( seq, conversionFlag, start_index, coreMotifSearchFlag, scores_list, 1, peak_name, entry_name, topN, fileSelection) #antiparallel seach
    
        #####
        else:
            seq = getData(genome_selected, chr_no, startLocation, endLocation)
            seq = convertDNAStrand(seq) #seq_negative
            start_index = startLocation
                
            scores_list = searchinSeq( seq, conversionFlag, start_index, coreMotifSearchFlag, scores_list, 0, peak_name, "", topN, fileSelection) #parallel search
            scores_list = searchinSeq( seq, conversionFlag, start_index, coreMotifSearchFlag, scores_list, 1, peak_name, "", topN, fileSelection) #antiparallel seach
            
    if geneName != None and coreMotifSearchFlag == True:
        search_type = "MOTIF SEARCH IN PEAK SEQUENCES OF "+geneName+" GENE with CORE MOTIF"
    elif geneName == None and coreMotifSearchFlag == True:
        search_type = "MOTIF SEARCH IN PEAK SEQUENCES with CORE MOTIF"
    elif geneName != None and coreMotifSearchFlag == False:
        search_type = "MOTIF SEARCH IN PEAK SEQUENCES OF "+geneName
    elif geneName == None and coreMotifSearchFlag == False:
        search_type = "MOTIF SEARCH IN PEAK SEQUENCES"

    outputResult(scores_list, coreMotifSearchFlag, search_type, topN)


def navigateSearch(peaks, gene_name):
    """**************************************************************SEARCH IN PEAK SEQUENCES******************************************************************"""
    #if searchInWholeSeqFlag == "on":
    # searched in negative strand whole sequence, position found wrt positive strand
    result_df = []
    
    # search motifs in whole sequence 
    if len(jasparfileName) > 0:
        if (topN):
            searchCallHelper(gene_name,1,result_df, peaks, topN, jasparfileName)
        else:
            searchCallHelper(gene_name,1,result_df, peaks, fileSelection = jasparfileName)

def displayGenomeMenu():
    resp = requests.get('https://api.genome.ucsc.edu/list/ucscGenomes')
    if resp.status_code != 200:
        # This means something went wrong.
        if resp.json()['error']:
            print("ERROR:",resp.json()['error'])
        else:
            print("ERROR while connecting")
        sys.exit(1)
    genomes_list = list(resp.json()['ucscGenomes'].keys())
    print( genomes_list)
    return genomes_list

def displayChrMenu(genome_selected,displayFlag):
    url = "api.genome.ucsc.edu/getData/track?genome={genome_name};track=refGene".format(genome_name=genome_selected) 
    resp = requests.get('https://'+url)
    if resp.status_code != 200:
        # This means something went wrong.
        if resp.json()['error']:
            print("ERROR:",resp.json()['error'])
        else:
            print("ERROR while connecting")
        sys.exit(1)

    chr_dict = resp.json()['refGene']
    if displayFlag == 1:
        print(chr_dict.keys())
    return chr_dict

def findPeakLocations(peakIntervalThr, genome_selected, cellLine, chr_no, startLocation, endLocation, tf_name):
    # read chipseq peak file
    url = "motifgenie-merge-peaks.ue.r.appspot.com/?tf={tf}&cellline={cl}&locus={chr}:{start}-{end}&pctg={thr}".format(tf=tf_name, cl=cellLine, chr=chr_no, start=startLocation, end=endLocation,thr=peakIntervalThr)

    resp = requests.get('https://'+url)
    if resp.status_code != 200:
        # This means something went wrong.
        if resp.json()['error']:
            print("ERROR:",resp.json()['error'])
        else:
            print("ERROR while connecting")
        sys.exit(1)
    
    html_file = BeautifulSoup(resp.text, "html.parser")

    peaks = []
    
    for interval in html_file.prettify().split("<pre>")[1].split("</pre>")[0].split("\n")[:-1]:
        items = interval.split("\t")
        peaks.append((items[0], items[1], items[2],items[3],items[4]))

    return peaks

def getTFname(pfm_file):
    file_name = open(pfm_file,"r")
    contents = file_name.read()

    # form a matrix from read file
    contents = contents.split("\n")
    tf_name = contents[0].split(" ")[1]

    return tf_name
 

def run(args):
    # ---CONFIGURE INPUTS---
    geneName = args.geneName

    global jasparfileName
    global topN
    global thresholdRatioForConsensus
    global coreMotifSearchFlag
    global seqLogoFlag
    global cellLine
    global genome_selected #as cistrome uses hg38, but can be changed later TODO
    global upperScoreThr
    global lowerScoreThr
    jasparfileName = args.jasparFile

    topN = args.topN
    thresholdRatioForConsensus = args.thresholdRatioForConsensus
    coreMotifSearchFlag = args.coreMotifSearchFlag
    seqLogoFlag = args.seqLogoFlag

    cellLine = args.cellLine
    genome_selected = "hg38"
    upperScoreThr = args.upperScoreThreshold
    lowerScoreThr = args.lowerScoreThreshold

    peakIntervalThr = args.peakIntervalThreshold
 
    locus = args.locus
    # get chr_no, startLocation, endLocation
    chr_no = locus.split(':')[0]
    locus = locus.replace(',','')
    startLocation = int(locus.split(':')[1].split('-')[0])
    endLocation = int(locus.split(':')[1].split('-')[1]) 

    #1. if genome isn't inputted, alert, display genome menu, get user input until its ok
    '''
    if genome_selected == None:
        print("Genome selection is not made. Please write one genome name from this menu:")
        genomes_list = displayGenomeMenu()
        genome_selected = input ("Enter genome:") 
        while genome_selected not in genomes_list:
            genome_selected = input ("Enter genome:") 
        #genome_selected OK
    '''
    
    #2. upperscorethr > lowerscorethr
    if upperScoreThr != None and lowerScoreThr != None:
        while upperScoreThr < lowerScoreThr:
            print("Lower score threshold cannot exceed the upper score threshold. Please rewrite them:")
            lowerScoreThr = input ("Enter the lower score threshold:")  
            upperScoreThr = input ("Enter the upper score threshold:")  
    
    #3. peak theshold must be within 0-100
    while peakIntervalThr < 0 or peakIntervalThr > 100:
        print("Peak interval threshold must be valid. Please rewrite it:")
        peakIntervalThr = input ("Enter the peak interval threshold:")  
    
    tf_name = getTFname(jasparfileName)
    peaks = findPeakLocations(peakIntervalThr, genome_selected, cellLine, chr_no, startLocation, endLocation, tf_name)
    if len(peaks) == 0:
        print("No Peaks Found, cannot do analysis.")
        sys.exit(0)

    navigateSearch(peaks, geneName)

def main():
    parser=argparse.ArgumentParser(description="MotifGenie is a Python command line tool for searching transcription factor binding motif sequences in merged ChIP-Seq binding sites. First, MotifGenie analyzes multiple ChIP-Seq samples simultaneously to identify shard binding regions for a locus of interest, then, searches the binding motif sequence of a given transcription factor in these shared binding regions.")

    parser.add_argument("-jf", "--jasparFile", help="pfm file name from JASPAR.", dest="jasparFile", type=str, required=True)
    parser.add_argument("-cl", "--cellLine", help="Cell line input.",dest="cellLine",type=str,required=True)
    parser.add_argument("-locus", "--locus", help="Write as: chr5:1,234,234-1,245,345",dest="locus",type= str,required=True)
    parser.add_argument("-pit", "--peakIntervalThreshold", help="Threshold of DNA binding sequence To be used in Peak Locations File to find Peak Sequence intervals.", dest="peakIntervalThreshold", type=int, required=True)
    
    parser.add_argument("-ust", "--upperScoreThreshold", help="Upper score threshold for displaying output.", dest="upperScoreThreshold", type=float)
    parser.add_argument("-lst", "--lowerScoreThreshold", help="Lower score threshold for displaying output.",dest="lowerScoreThreshold",type=float)
    parser.add_argument("-gene", "--geneName", help="Gene name to do search on.",dest="geneName",type=str)

    # user input numerical values and thresholds
    parser.add_argument("-n", "--topN", help="Number of motifs to be displayed which have the highest scores. (15 by default)", dest="topN", type=int)
    parser.add_argument("-trc", "--thresholdRatioForConsensus", help="Positions of consensus get trimmed if frequency ratio is greater than the thresholdRatioForConsensus. Used with coreMotifSearchFlag. (0.66 by default)", dest="thresholdRatioForConsensus", type=float)

    group1 = parser.add_mutually_exclusive_group()
    group1.add_argument("-cms","--coreMotifSearchFlag",action="store_true",help="Runs the search with Core Motifs", dest="coreMotifSearchFlag")

    group2 = parser.add_mutually_exclusive_group()
    group2.add_argument("-sl","--seqLogoFlag",action="store_true",help="Displays the Sequence Logo after each motif search when the argument is used.", dest="seqLogoFlag")

    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main()
