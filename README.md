# FastDup

Identifies duplicate reads. This tool locates and tags duplicate reads in a coordinate ordered SAM or BAM file.

Use the same algorithm as picard MarkDuplicates and output identical results.
Use spdlog as log tool and the default level is 'info'.

### Features

* Fast - with the same number of threads `FastDup` is ~8X faster than GATK MarkDuplicatesSpark.
  And `FastDup` achives ~20X performance improvement than Picard MarkDuplicates.
* Generate identical outputs compared to Picard MarkDuplicates.
* The same detailed metrics data witch Picard MarkDuplicates.
* All data processed in memory and low-memory footprint even for large input files. 

### Limitations

* Although `FastDup` can detecte all the same duplicates as Picard MarkDuplicates. They may mark 
  different reads as duplicates because the reads sort algorithm in Picard MarkDuplicates is unstable.
  Considering there are 2 reads(A, B and A is in front of B in file) in a duplicate group and they
  have the same score, Picard Markduplicates may mark A as duplicate because B may be in front of A
  after sorting. While `FastDup` use stable sort algorithm and always mark B as duplicate.
* In optical duplicates detection, Picard Markduplicates use short (int16_t) as data type in parsing 
  tile/region, x coordinate and y coordinate from a read name, which may data overflow as these integers
  may exceed the range of short type. `FastDup` fixes this bug. But for consistency with Picard Markduplicates,
  we keep this bug in source codes. Just change the data type in PhysicalLocation struct in read_ends.h file
  to fix this bug.
* `FastDup` use the data characteristics in coordinate ordered SAM/BAM files to improve the performance of
  detecting duplicates, thus the input should be ordered by coordinate in advance.

## Requirements

### Build tools

* autoconf (for htslib)
* cmake
* c++17 (gcc >= 8.1 or clang >= 7 should work.)

### Libraries needed

* zlib
* libbz2
* liblzma
* libcurl
* libdeflate (optional)

## Install

Download a distribution tarball `FastDup.tar.gz` or clone the source codes from github.

```bash
# build htslib
cd FastDup/ext/htslib
autoreconf -i
./configure
make

# build FastDup
cd FastDup
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

The generated binary fastdup will be in the build/bin folder.

## Usage

Get help

```bash
./fastdup --help
```

Mark duplicates on an input BAM file using 8 threads

```bash
./fastdup --input in_test.bam --output out_md.bam --metrics stats.txt --num-threads 8
```