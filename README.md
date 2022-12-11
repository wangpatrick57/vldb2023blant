# Table of Contents
 1. Outside Links
 2. How to Run
 3. Source Code Overview

# Outside Links
arXiv link to paper: https://arxiv.org/abs/2212.02771<br>
Active repo of indexing stage: https://github.com/waynebhayes/BLANT/blob/master/src/blant.c<br>
Active repo of alignment and merging stages: https://github.com/wangpatrick57/plant<br>

_This repo is a frozen version of the previous two repos which represents the relevant code for the version of the paper in the arXiv link. Significants parts of the code (old experiments or unrelated functionality) are removed to aid readability_

# How to Run
### Docker
We have provided an Ubuntu Docker environment to run the code in. Use `./run_docker.sh` to start an interactive Docker container. This docker container has a volume to the root directory of the repository, so source code may be easily modified.

### Intermediate Files
We have provided files containing intermediate output between each stage so that they may all be run independently. Example output of stage 1 is in indexes/ and example output of stage 2 is in alignments/.

### Stage 1: Indexing
This stage is written in C and will need to be compiled with a Makefile we have provided. Running `make` for the first time will sleep for 100 seconds and show some instructions. `make clean` will clean code-related files and `make pristine` will remove all other files generated in the make process, resulting in a pristine repo. `PAUSE=0 make base` will make everything, but will take up to an hour or more to run as the k=8 files in canon_maps/ take a long time to generate. If you want to just try indexing with k<=7, you may run `PAUSE=0 NO8=1 make base`.

After making, an executable called "blant" will be created. For historical reasons, this name is slightly confusing: the exectuable "blant" represents the first stage in the algorithm called "BLANT" in the paper. This executable has two parameters, k and D, described in the paper. All the .ind files in indexes/ were created with k=8 and D=2, as described in the paper. blant performs a lot of deduplication but does not do it perfectly, so deduplication must be done on the final output file as well with `dedup.sh`. An example run of this stage looks like `./blant -k8 -D2 networks/mouse.el >indexes/mouse.ind; dedup.sh indexes/mouse.ind`. Note that this command will take about 30 minutes to run, so a faster way to test blant out would be to run `./blant -k6 -D1 networks/syeast.el`.

### Stage 2: Alignment
The code to run this stage is in alignment.py. alignment.py takes in four parameters: two indexes and two networks. An example run of this stage looks like `python3 alignment.py indexes/mouse.ind networks/mouse.el indexes/rat.ind networks/rat.el >alignments/mouse-rat.alns`.

### Stage 3: Merging
The code to run this stage is in merging.py. merging.py takes in three parameters: an alignments file and two networks. An example run of this stage looks like `python3 merging.py alignments/mouse-rat.alns networks/mouse.el networks/rat.el`.

# Source Code Overview
### C code for indexing
All the code used for this stage is in the src/ directory. Most files here are helper files used for various utility functions such as reading in a graph file, storing a graph object, utilizing canon_maps/, etc. The "core" logic of the indexing algorithm is in the `SampleGraphletIndexAndPrint()` function in blant.c.

### Python code for alignment and merging
All of the code used for these two stages are in .py files in the root directory of the repo. Many files are used for various utility functions like reading, printing, etc. The "core" logic for the alignment step is in `find_seeds()` function in seeding_algorithm_core.py. The "core" logic for the merging step is in the `SimAnnealGrow` class in simanneal_grow_core.py.
