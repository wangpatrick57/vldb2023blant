## How to Run
We have provided files containing intermediate output between each stage so that they may all be run independently. For example, stage 3 can be tested even if stages 1 and 2 are not yet working. Example output of stage 1 is in indexes/ and example output of stage 2 is in alignments/.

### Stage 1: Indexing
This stage is written in C, with the source code files being in src/. The relevant sections to read are the function SampleGraphletIndexAndPrint() in blant-sampling.c and the section in blant.c which calls SampleGraphletIndexAndPrint(). Other source files include various helper functions which allow the BLANT project to run, such as reading in edge list files, deduplicating graphlets, etc.

First, run "PAUSE=0 NO8=0 make base" in order to make the ./blant executable. Additional instructions about making ./blant can be found in https://github.com/waynebhayes/BLANT. After making, run "./blant -k8 -lDEG2 -mi -sINDEX [network]" and redirect the output to an index file. After doing so, you should deduplicate the lines in the index file without changing the order of the lines using "./dedup.sh [index]" We have put example index files in the indexes/ directory. The networks used to generate those index files are in networks/.

### Stage 2: Alignment
This stage is written in Python, and the source files are in the root directory of the repository. To run this stage, simply run "python3 alignment.py [index1] [network1] [index2] [network2]" and redirect the output to an alignments file. Example alignment files are in alignments/.

### Stage 2: Merging
This stage is also written in Python, and the source files are in the root directory of the repository. To run this stage, simply run "python3 merging.py [alignments] [network1] [network2]", and redirect the output to an output file. Example alignment files are in merged/. 10000 iterations is a good starting amount for any network pair that should run in a few minutes
