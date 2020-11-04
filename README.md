# Fingerprint Aligner
Fingerprint Aligner (FA) is an All vs One multiple threaded Multiple Sequence Alignment program.

## About the program
FA extracts a fingerprint from a sequence given as input, the reference sequence. The fingerprint is then utilised to align all other query sequences.
Therefore, to utilise this program you should have two files, one with the sequences you want to align, and a second one with the sequence you want to align those sequences to.

The program was first developed to rapidly align SARS-CoV-2 genome sequences, but its range of applications can go far beyond.
By using NCBI's SARS-CoV-2 sequence as the reference genome, it was possible to align more than 700 sequences in under 25 minutes.

# Installation
The installation of a binary in the system is still a work in progress, but for now you can build your own binary or utilise the ones which come with the project.
In the folders `cmake-build-debug` and `cmake-build-release` you can find the pre-built binaries of the program. If you have no intention of debugging the program, utilise the binary inside the release folder.

## Building
Clone the project on your machine, then build it by running `make` inside one of the build folders.
```
git clone https://github.com/ovaltriangle/fingerprint-aligner.git
cd fa/cmake-build-release
make
```

# Usage
```
Aligns genomes using a reference fingerprint.
Usage:
  fa [OPTION...]

  -r, --reference arg           The file containing the reference genomes to
                                be used.
  -q, --query arg               The file containing the query genomes(s).
  -o, --output arg              Sequence file output. If not specified,
                                'out.fa' will be used. (default: out.fa)
  -t, --threads arg             The number of threads to use (>1 for
                                multithreading). (default: 1)
  -p, --protein-threshold arg   The threshold for the number of nucleotides
                                before a protein is considered so. (default:
                                60)
  -s, --similarity-threshold arg
                                Minimum percentage of identity when comparing
                                two sequences. (default: 0.9)
  -w, --weights arg             The weights to be used in the Hirschberg's NW
                                score calculation(insertion, deletion,
                                substitution and match scores). Separate the values
                                with a comma. (default: -2,-2,-1,1)
  -h, --help                    Print this message and exit.
```
Fingerprint Aligner expects two files; a file containing the reference genome and a file containing all queries genomes to be aligned to the reference.
The program currently accepts only FASTA files and works only with DNA and RNA sequences.

You can specify the number of threads to be used by using the `-t` option. Be cautious in using large number of threads, as it can freeze or crash your machine.

## Example
All shown examples assume you have two files, `reference.fasta` and `queries.fasta`
```
# Basic usage
fa -r reference.fasta -q queries.fasta

# Extended options
fa --reference reference.fasta --query queries.fasta

# Align the sequences utilising 6 threads
fa -r reference.fasta -q queries.fasta -t 6

# Align the sequences utilising 6 threads and custom settings
fa -r reference.fasta -q queries.fasta -o result.fasta -t 6 -p 90 -s 0.95 -w -5,-5,-3,2
```
