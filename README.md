<img width="200px" src="actc.png" align="right"/>
<h1 align="center">actc</h1>
<p align="center">Align clr to ccs reads.</p>

***

# What is this about?
If you want to have one-click solution to align subreads to the respective CCS reads.

# How to install

## Use bioconda

```
conda install -c bioconda actc
```

## Download binary directly

  Get binary from [releases](https://github.com/PacificBiosciences/actc/releases) and make it executable.

## Compile manually

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
 * The CCS file must be a subset of the ZMWs of the subread input file.
 * The CLR file must be indexed, a `subreads.bam.pbi` can be generated via `pbindex`
 * The CLR file must be sorted by ZMW, default for `subreads.bam`

# Chunking
You can parallelize by chunking, via `--chunk` is possible. For this, the
`ccs.bam` file must accompanied by a `.pbi` file.

An example workflow, all `actc` invocations can run simultaneously:

    actc movie.subreads.bam movie.ccs.bam aligned.1.bam --chunk 1/10 -j <THREADS>
    actc movie.subreads.bam movie.ccs.bam aligned.2.bam --chunk 2/10 -j <THREADS>
    ...
    actc movie.subreads.bam movie.ccs.bam aligned.10.bam --chunk 10/10 -j <THREADS>

# How to index BAM files
To generate the BAM index
`.bam.pbi`, use `pbindex`, which can be installed with `conda install pbbam`.

    pbindex movie.ccs.bam
    pbindex movie.subreads.bam

# Changelog
  * 0.4.0
    * Accept REVIO files
  * 0.3.0
    * Properly handle flanks
  * 0.2.0
    * Improve runtime and reduce memory footprint for low-complexity molecules
  * 0.1.0
    * Add _actc_ `@PG` to BAM header
    * Auto generate version from `meson.build`
    * Improve (de)compression, ergo throughput, by threading
  * 0.0.3
    * Use subreads read groups
    * Fix FASTA output for `--chunk`
    * Initial attempt to increase throughput
  * 0.0.2
    * Add `--chunk` capabilities
  * 0.0.1
    * Intial release
