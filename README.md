# MotifGenie
MotifGenie is a Python command line tool for searching transcription factor binding motif sequences in merged ChIP-Seq binding sites. First, MotifGenie analyzes multiple ChIP-Seq samples simultaneously to identify shard binding regions for a locus of interest, then, searches the binding motif sequence of a given transcription factor in these shared binding regions.

*By:* Çerağ Oğuztüzün, Pelin Yaşar, Kerim Yavuz, Mesut Muyan, Tolga Can

## Requirements 
beautifulsoup4-4.9.3
`$ pip3 install bs3` 

logomaker-0.8 
`$ pip3 install logomaker`
## Installation  
Download the mg.py script in a working directory.
## Input Files 
MotifGenie requires a pfm file from the JASPAR database as input, hence the mg.py script must be in the same directory as the pfm file. Pfm file contains a Position Frequency Matrix of the motif to be searched.

An example pfm file of the ESR1 transcription factor from the JASPAR database, where MA0112.1 refers to the ID:
```
>MA0112.1 ESR1
  1.00   1.00   7.00   2.00   0.00   0.00   0.00   6.00   1.00   2.00   3.00   1.00   1.00   5.00   0.00   0.00   2.00   3.00
  5.00   5.00   1.00   0.00   0.00   0.00   7.00   0.00   7.00   5.00   2.00   1.00   0.00   1.00   8.00   9.00   4.00   4.00
  1.00   1.00   1.00   7.00   9.00   0.00   2.00   2.00   1.00   1.00   4.00   1.00   8.00   3.00   0.00   0.00   0.00   1.00
  2.00   2.00   0.00   0.00   0.00   9.00   0.00   1.00   0.00   1.00   0.00   6.00   0.00   0.00   1.00   0.00   3.00   1.00
```
## How To Run
### Parameters
MotifGenie is used by a set of required and optional arguments through the command line interface. MotifGenie does search on the hg38 genome. The cell line information to be inputted must be present in CellLine_list.txt found in the repository. Also, the name of the Transcription Factor pfm file to be inputted must be present in the TF_list.txt file which is found in the repository. If the TF name is not in the list, peaks will not be found.

The **required** arguments include:
- `-jf JASPARFILE` or `--jasparFile JASPARFILE`: Pfm file name from JASPAR. (for example MA0112.1.pfm)
- `-cl CELLLINE` or `--cellLine CELLLINE`: Cell line input. (for example MCF-7)
- `-locus LOCUS` or `--locus LOCUS`: Written as: chr5:1,234,234-1,245,345 or chr5:1234234-1245345
- `-pit PEAKINTERVALTHRESHOLD` or `--peakIntervalThreshold PEAKINTERVALTHRESHOLD`: Threshold of DNA binding sequence to be used in peak locations to find peak sequence intervals. This means filtering the peak sequences where binding percentage passes the inputted threshold value.

The **optional** arguments include:
- `-h` or `--help`: displays the documentation
- `-gene GENENAME` or `--geneName GENENAME`: Gene name to search using UCSC Genome Browser.
- `-ust UPPERSCORETHRESHOLD` or `--upperScoreThreshold UPPERSCORETHRESHOLD`: Upper score threshold for displaying output.
- `-lst LOWERSCORETHRESHOLD` or `--lowerScoreThreshold LOWERSCORETHRESHOLD`: Lower score threshold for displaying output.
- `-n TOPN` or `--topN TOPN`: Number of motifs to be displayed that have the highest scores. (15 by default)
- `-cms` or `--coreMotifSearchFlag`: Runs the search with Core Motifs.
- `-trc THRESHOLDRATIOFORCONSENSUS` or `--thresholdRatioForConsensus THRESHOLDRATIOFORCONSENSUS`: Positions of consensus get trimmed if the frequency ratio is greater than the inputted threshold. Used with coreMotifSearchFlag, when searching using Core Motifs. (0.66 by default)
- `-sl` or `--seqLogoFlag`: Displays and saves the Sequence Logo after motif search when the argument is used.

The command line interface structure: 
```
mg.py [-h] -jf JASPARFILE -cl CELLLINE -locus LOCUS -pit
             PEAKINTERVALTHRESHOLD [-ust UPPERSCORETHRESHOLD]
             [-lst LOWERSCORETHRESHOLD] [-gene GENENAME] [-n TOPN]
             [-trc THRESHOLDRATIOFORCONSENSUS] [-cms] [-sl]
```
### Sample Search Commands
- This searches the motif given through the pfm file name, on the given locus using the cell line and the peak interval threshold (which are the must-have arguments).

`python mg.py -jf MA0112.1.pfm -cl MCF-7 -locus chr5:139638349-139693882 -pit 20`
- This search is similar to the above one but does a core motif search by inputting the threshold ratio for consensus (-trc) value.

`python mg.py -jf MA0112.1.pfm -cl MCF-7 -locus chr5:139638349-139693882 -pit 20 -trc 0.7 -cms`
- This does search on the overlapping sequences of CXXC5 gene's sequences obtained from the UCSC database, and the peak sequences. The search also uses core motif search and specifies the number of outputs (-n) to be 35 motifs.

`python mg.py -jf MA0112.1.pfm -cl MCF-7 -locus chr5:139638349-139693882 -pit 20 -n 35 -gene CXXC5 -trc 0.7 -cms`
- This search is similar to the above one, but specifies the upper score threshold (-ust) and lower score threshold (-lst) to filter the motifs to be outputted. The sequence logo flag (-sl) can be included to each of these commands for viewing the Sequence Logo of the outputted set of motifs.

`python mg.py -jf MA0112.1.pfm -cl MCF-7 -locus chr5:139638349-139693882 -pit 20 -n 30 -gene CXXC5 -ust 8.3 -lst 6.9`
## Output Description 
MotifGenie outputs the search results as a table in the command line interface. If the Sequence Logo flag is used, the generated Sequence Logo will be saved in the directory of the script as a PNG file.

A sample output:
```
--------------------------------------------------------------------------------------------------------------------------
MOTIF SEARCH IN PEAK SEQUENCES OF CXXC5 GENE with CORE MOTIF
Threshold to trim consensus = 0.7
--------------------------------------------------------------------------------------------------------------------------
SEQUENCE                      SCORE               POSITION                      PARALLELITY    PEAK           ENTRY
--------------------------------------------------------------------------------------------------------------------------
GGCCAGGGTCACCT                18.064151           139651263 - 139651276         -              peak5          NM_001317200
GGCCAGGGTCACCT                18.064151           139651263 - 139651276         -              peak5          NM_001317201
GGCCAGGGTCACCT                18.064151           139651263 - 139651276         -              peak5          NM_001317205
GGCCAGGGTCACCT                18.064151           139651263 - 139651276         -              peak5          NM_001317204
GGCCAGGGTCACCT                18.064151           139651263 - 139651276         -              peak5          NM_001317202
```

- `HEADER`: the summary of the search parameters such as, whether a specific gene name was specified, core motif search was used and the inputted consensus trimming threshold.
- `SEQUENCE`: the sequence of motifs found sorted with respect to the scores.
- `SCORE`: the log-based score given to the motif which is correlated with the similarity of the inputted motif that was searched on the reference sequence.
- `POSITION`: the start and end locations where the motif was found on the reference genome.
- `PARALLELITY`: whether the motif was found in a parallel (+) or antiparallel (-) search.
- `PEAK`: the number of peak sequence where the motif was found on.
- `ENTRY`: the entry number referring to the UCSC database, in which the sequence for the inputted gene was found. This column is only presented for gene-specific search.

Sample Sequence Logo output (not related to the above result):
![MA0112 1 pfm-MOTIF SEARCH IN PEAK SEQUENCES OF CXXC5](https://user-images.githubusercontent.com/38559757/109000332-9fac6b80-76b4-11eb-97a3-3d7e2a7e1623.png)

