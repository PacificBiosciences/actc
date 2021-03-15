<h1 align="center"><img width="300px" src="actc.png"/></h1>
<h1 align="center">actc</h1>
<p align="center">Align clr to ccs reads.</p>

# What is this about?
If you want to have one-click solution to align subreads to the respective CCS reads.

# How to compile

    mkdir build && cd build
    meson --prefix ~/mytools
    ninja
    ninja install

# How to run

    actc movie.subreads.bam movie.ccs.bam aligned.bam --log-level INFO

# Output files
Main output file `aligned.bam` contains all alignments,
in the same order as they occur in the inputs.

Auxilliary file `aligned.fasta` contains all references of the alignment file.

# Pre-conditions
The CCS file must be a subset of the ZMWs of the subread input file.
