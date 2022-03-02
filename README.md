[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/jermp/sshash.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/jermp/sshash/context:cpp)

Notes on mindex development
===========================

Currently, the reference-based indexing code is not well-integrated and there are a few cumbersome steps to get the full index.  We'll provide an example below, assuming we start with a reference genome sequence in a fasta file called `ref.fa`.  Further, we will assume for the time being that `ref.fa` contains *no invalid characters* (in practice this is achieved by running a genome through the `fixFasta` functionality of pufferfish).

First we need to build the compacted reference de Bruijn Graph.  We do this using [`cuttlefish`](https://github.com/COMBINE-lab/cuttlefish).

```
$ cuttlefish build -s ref.fa -k 31 -t 16 -o ref_dbg -m 8 -f 3
```

This creates 2 files called `ref_dbg.cf_seq` and `ref_dbg.cf_seg`.  Next we will build sshash on the `ref_dbg.cf_seg` file.

```
$ ./build <( awk '{ print ">"$1"\n"$2; }' ref_dbg.cf_seg ) 31 20 --canonical-parsing -o ref_idx.sshash
```
We are using input redirection here to convert the `cf_seg` file to a fasta file on the fly.  Next we will build the contig table:

```
$ ./build_contig_table ref_dbg 31 -o ref_idx.ctab
```

Now, we have a complete index consisting of the sequence index (the `.sshash` file) and the contig index (the `.ctab` file).  Right now,
the reference index loader assumed these specific suffixes.  Of course, this is fragile and error prone, so we'll make it better.  Finally
you can test out the index with:

```
$ ./test_load ref_idx ref.fa
```

This will iterate over all of the positions of `ref.fa` and ensure that the observed k-mer has an entry in the reference index that points
back to the current reference sequence, at the current position, in the provided (forward) orientation.


SSHash
======

This is a compressed dictionary data structure for k-mers
(strings of length k over the DNA alphabet {A,C,G,T}), based on **S**parse and **S**kew **Hash**ing.

**A (pre-print) paper describing the data structure can be found [here](https://www.biorxiv.org/content/10.1101/2022.01.15.476199).**

For a dictionary of n k-mers,
two basic queries are supported:

- i = Lookup(g), where i is in [0,n) if the k-mer g is found in the dictionary or i = -1 otherwise;
- g = Access(i), where g is the k-mer associated to the identifier i.

If also the abundances of the k-mers (their frequency counts) are stored in the dictionary, then the dictionary is said to be *weighted* and it also supports:

- c = Abundance(i), where i is a given k-mer identifier.

A membership query (determine if a given k-mer is present in the dictionary or not) is, therefore, supported by means of the lookup query.
The dictionary can also stream through all k-mers of a given DNA file
(.fasta or .fastq formats) to determine their membership to the dictionary.

**NOTE**: The Lookup query assumes that two k-mers being the *reverse complement* of each other are the same.

#### Table of contents
* [Compiling the Code](#compiling-the-code)
* [Dependencies](#dependencies)
* [Build a Dictionary](#build-a-dictionary)
* [Examples](#Examples)
* [Input Files](#input-files)
* [Large-scale Benchmark](#large-scale-benchmark)
* [Author](#author)

Compiling the Code
------------------

The code is tested on Linux with `gcc` and on Mac with `clang`.
To build the code, [`CMake`](https://cmake.org/) is required.

Clone the repository with

	git clone --recursive https://github.com/jermp/sshash.git

If you have cloned the repository without `--recursive`, you will need to perform the following commands before
compiling:

    git submodule init
    git submodule update

To compile the code for a release environment (see file `CMakeLists.txt` for the used compilation flags), it is sufficient to do the following:

    mkdir build
    cd build
    cmake ..
    make -j

For a testing environment, use the following instead:

    mkdir debug_build
    cd debug_build
    cmake .. -D CMAKE_BUILD_TYPE=Debug -D SSHASH_USE_SANITIZERS=On
    make

Dependencies
------------

The repository has minimal dependencies: it only uses the [PTHash](https://github.com/jermp/pthash) library (for minimal perfect hashing), and `zlib` to read gzip-compressed streams.

To automatically pull the PTHash dependency, just clone the repo with
`--recursive` as explained in [Compiling the Code](#compiling-the-code).

If you do not have `zlib` installed, you can do

	sudo apt-get install zlib1g

if you are on Linux/Ubuntu, or

	brew install zlib

if you have a Mac.

Build a Dictionary
------------------

The driver program
called `build` can be used to build a dictionary.

From within the directory
where the code was compiled (see the section [Compiling the Code](#compiling-the-code)), run the command:

	./build --help

to show the usage of the driver program (reported below for convenience).

	Usage: ./build [-h,--help] input_filename k m [-s seed] [-l l] [-c c] [--canonical-parsing] [--abundances] [-o output_filename] [--check] [--bench] [--verbose]

	 input_filename
		Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not:
		- without duplicate nor invalid kmers
		- one DNA sequence per line.
		For example, it could be the de Bruijn graph topology output by BCALM.

	 k
		K-mer length (must be <= 31).

	 m
		Minimizer length (must be < k).

	 [-s seed]
		Seed for construction (default is 1).

	 [-l l]
		A (integer) constant that controls the space/time trade-off of the dictionary. A reasonable values lies between 2 and 12 (default is 6).

	 [-c c]
		A (floating point) constant that trades construction speed for space effectiveness of minimal perfect hashing. A reasonable value lies between 3.0 and 10.0 (default is 3.000000).

	 [--canonical-parsing]
		Canonical parsing of k-mers. This option changes the parsing and results in a trade-off between index space and lookup time.

	 [--abundances]
		Also store the abundances in compressed format.

	 [-o output_filename]
		Output file name where the data structure will be serialized.

	 [--check]
		Check correctness after construction.

	 [--bench]
		Run benchmark after construction.

	 [--verbose]
		Verbose output during construction.

	 [-h,--help]
		Print this help text and silently exits.
		

Examples
--------

For the examples, we are going to use some collections
of *stitched unitigs* from the directory `../data/unitigs_stitched`.
These collections were built for k = 31, so dictionaries should be built with k = 31 as well to ensure correctness.

In the section [Input Files](#input-files), we explain how
such collections of stitched unitigs can be obtained from raw FASTA files.

### Example 1

	./build ../data/unitigs_stitched/salmonella_enterica_k31_ust.fa.gz 31 13 --check --bench -o salmonella_enterica.index

This example builds a dictionary for the k-mers read from the file `../data/unitigs_stitched/salmonella_enterica_k31_ust.fa.gz`,
with k = 31 and m = 13. It also check the correctness of the dictionary (`--check` option), run a performance benchmark (`--bench` option), and serializes the index on disk to the file `salmonella_enterica.index`.

To run a performance benchmark after construction of the index,
use:

	./bench salmonella_enterica.index

To also store the abundances, use the option `--abundances`:

	./build ../data/unitigs_stitched/with_abundances/salmonella_enterica_k31_ust.abundances.fa.gz 31 13 --abundances --check --verbose
	
### Example 2

	./build ../data/unitigs_stitched/salmonella_100_k31_ust.fa.gz 31 15 -l 2 -o salmonella_100.index

This example builds a dictionary from the input file `../data/unitigs_stitched/salmonella_100_k31_ust.fa.gz` (a pangenome consisting in 100 genomes of *Salmonella Enterica*), with k = 31, m = 15, and l = 2. It also serializes the index on disk to the file `salmonella_100.index`.

To perform some streaming membership queries, use:

	./query salmonella_100.index ../data/queries/SRR5833294.10K.fastq.gz

if your queries are meant to be read from a FASTQ file, or

	./query salmonella_100.index ../data/queries/salmonella_enterica.fasta.gz --multiline
	
if your queries are to be read from a (multi-line) FASTA file.

### Example 3

	./build ../data/unitigs_stitched/salmonella_100_k31_ust.fa.gz 31 13 -l 4 -s 347692 --canonical-parsing -o salmonella_100.canon.index

This example builds a dictionary from the input file `../data/unitigs_stitched/salmonella_100_k31_ust.fa.gz` (same used in Example 2), with k = 31, m = 13, l = 4, using a seed 347692 for construction (`-s 347692`), and with the canonical parsing modality (option `--canonical-parsing`). The dictionary is serialized on disk to the file `salmonella_100.canon.index`.

The	 "canonical" version of the dictionary offers more speed for only a little space increase (for a suitable choice of parameters m and l), especially under low-hit workloads -- when the majority of k-mers are not found in the dictionary. (For all details, refer to the paper.)
	
Below a comparison between the dictionary built in Example 2 (not canonical)
and the one just built (Example 3, canonical).

	./query salmonella_100.index ../data/queries/SRR5833294.10K.fastq.gz
	index size: 10.3981 [MB] (6.36232 [bits/kmer])
	==== query report:
	num_kmers = 460000
	num_valid_kmers = 459143 (99.8137% of kmers)
	num_positive_kmers = 46 (0.0100187% of valid kmers)
	num_searches = 42/46 (91.3043%)
	num_extensions = 4/46 (8.69565%)
	elapsed = 229.159 millisec / 0.229159 sec / 0.00381932 min / 498.172 ns/kmer
	
	./query salmonella_100.canon.index ../data/queries/SRR5833294.10K.fastq.gz
	index size: 11.0657 [MB] (6.77083 [bits/kmer])
	==== query report:
	num_kmers = 460000
	num_valid_kmers = 459143 (99.8137% of kmers)
	num_positive_kmers = 46 (0.0100187% of valid kmers)
	num_searches = 42/46 (91.3043%)
	num_extensions = 4/46 (8.69565%)
	elapsed = 107.911 millisec / 0.107911 sec / 0.00179852 min / 234.589 ns/kmer

We see that the canonical dictionary is twice as fast as the regular dictionary
for low-hit workloads,
even on this tiny example, for only +0.4 bits/k-mer.

Input Files
-----------

SSHash is meant to index k-mers from collections that do not contain duplicates
nor invalid k-mers (strings containing symbols different from {A,C,G,T}).
These collections can be obtained, for example, by extracting the maximal unitigs of a de Bruijn graph.

Doing so is easy to do using the tool [BCALM2](https://github.com/GATB/bcalm).
This tool builds a compacted de Bruijn graph and outputs its maximal unitigs.
From the output of BCALM2, we can then *stitch* (i.e., glue) some unitigs to reduce the number of nucleotides. The stitiching process is carried out using the [UST](https://github.com/medvedevgroup/UST) tool.

Below we provide a complete example (assuming both BCALM2 and UST are installed correctly) that downloads the Human (GRCh38) Chromosome 13 and extracts the maximal stitiched unitigs for k = 31.

	mkdir DNA_datasets
	wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz -O DNA_datasets/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz
	~/bcalm/build/bcalm -in ~/DNA_datasets/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz -kmer-size 31 -abundance-min 1 -nb-cores 8
	~/UST/ust -k 31 -i ~/Homo_sapiens.GRCh38.dna.chromosome.13.fa.unitigs.fa
	gzip Homo_sapiens.GRCh38.dna.chromosome.13.fa.unitigs.fa.ust.fa
	rm ~/Homo_sapiens.GRCh38.dna.chromosome.13.fa.unitigs.fa

See also the script `scripts/download_and_preprocess_datasets.sh`
for precise arguments.

Large-scale Benchmark
---------------------

*Pinus Taeda* ("pine", [GCA_000404065.3](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/404/065/GCA_000404065.3_Ptaeda2.0/GCA_000404065.3_Ptaeda2.0_genomic.fna.gz)) and *Ambystoma Mexicanum* ("axolotl", [GCA_002915635.2](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/915/635/GCA_002915635.3_AmbMex60DD/GCA_002915635.3_AmbMex60DD_genomic.fna.gz))
are some of the largest genome assemblies, respectively counting
10,508,232,575 and 17,987,935,180 distinct k-mers for k = 31.

After running BCALM2 and UST, we build the indexes as follows.

    ./build ~/DNA_datasets.larger/GCA_000404065.3_Ptaeda2.0_genomic.ust_k31.fa.gz 31 20 -l 6 -c 7 -o pinus.m20.index
    ./build ~/DNA_datasets.larger/GCA_000404065.3_Ptaeda2.0_genomic.ust_k31.fa.gz 31 19 -l 6 -c 7 --canonical-parsing -o pinus.m19.canon.index
    ./build ~/DNA_datasets.larger/GCA_002915635.3_AmbMex60DD_genomic.ust_k31.fa.gz 31 21 -l 6 -c 7 -o axolotl.m21.index
    ./build ~/DNA_datasets.larger/GCA_002915635.3_AmbMex60DD_genomic.ust_k31.fa.gz 31 20 -l 6 -c 7 --canonical-parsing -o axolotl.m20.canon.index
    
The following table summarizes the space of the dictionaries.

| Dictionary        |Pine       || Axolotl  ||
|:------------------|:---:|:----:|:---:|:---:|
|                   | GB     | bits/k-mer  | GB    | bits/k-mer |
| SSHash, regular   | 13.21  | 10.06       | 22.28 | 9.91       |
| SSHash, canonical | 14.94  | 11.37       | 25.03 | 11.13      |

    

To query the dictionaries, we use [SRR17023415](https://www.ebi.ac.uk/ena/browser/view/SRR17023415) fastq reads
(23,891,117 reads, each of 150 bases) for the pine,
and [GSM5747680](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5747680) multi-line fasta (15,548,160 lines) for the axolotl.

Timings have been collected on an Intel Xeon Platinum 8276L CPU @ 2.20GHz,
using a single thread.

| Dictionary        |Pine       || Axolotl  ||
|:------------------|:---:|:----:|:---:|:---:|
|                   |(>75% hits)||(>86% hits)|
|                   | tot (min) | avg (ns/k-mer) | tot (min) | avg (ns/k-mer) |
| SSHash, regular   | 19.2      | 400            | 4.2       | 269            |
| SSHash, canonical | 14.8      | 310            | 3.2       | 208            |

Below the complete query reports.

    ./query pinus.m20.index ~/DNA_datasets.larger/queries/SRR17023415_1.fastq.gz
    ==== query report:
    num_kmers = 2866934040
    num_valid_kmers = 2866783488 (99.9947% of kmers)
    num_positive_kmers = 2151937575 (75.0645% of valid kmers)
    num_searches = 418897117/2151937575 (19.466%)
    num_extensions = 1733040458/2151937575 (80.534%)
    elapsed = 1146.58 sec / 19.1097 min / 399.933 ns/kmer

    ./query pinus.m19.canon.index ~/DNA_datasets.larger/queries/SRR17023415_1.fastq.gz
    ==== query report:
    num_kmers = 2866934040
    num_valid_kmers = 2866783488 (99.9947% of kmers)
    num_positive_kmers = 2151937575 (75.0645% of valid kmers)
    num_searches = 359426304/2151937575 (16.7025%)
    num_extensions = 1792511271/2151937575 (83.2975%)
    elapsed = 889.779 sec / 14.8297 min / 310.359 ns/kmer

    ./query axolotl.m21.index ~/DNA_datasets.larger/queries/Axolotl.Trinity.CellReports2017.fasta.gz --multiline
    ==== query report:
    num_kmers = 931366757
    num_valid_kmers = 748445346 (80.3599% of kmers)
    num_positive_kmers = 650467884 (86.9092% of valid kmers)
    num_searches = 124008258/650467884 (19.0645%)
    num_extensions = 526459626/650467884 (80.9355%)
    elapsed = 250.173 sec / 4.16955 min / 268.608 ns/kmer

    ./query axolotl.m20.canon.index ~/DNA_datasets.larger/queries/Axolotl.Trinity.CellReports2017.fasta.gz --multiline
    ==== query report:
    num_kmers = 931366757
    num_valid_kmers = 748445346 (80.3599% of kmers)
    num_positive_kmers = 650467884 (86.9092% of valid kmers)
    num_searches = 106220473/650467884 (16.3299%)
    num_extensions = 544247411/650467884 (83.6701%)
    elapsed = 193.871 sec / 3.23119 min / 208.158 ns/kmer
    
Author
------

Giulio Ermanno Pibiri - <giulio.ermanno.pibiri@isti.cnr.it>
