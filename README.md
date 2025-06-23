# FastDup: A Scalable Duplicate Marking Tool using Speculation-and-Test Mechanism

FastDup is a tool designed to locate and tag duplicate reads in a coordinate-sorted SAM or BAM file. It uses the same core algorithm as Picard MarkDuplicates to produce identical results and utilizes `spdlog` for logging, with the default level set to 'info'.

## ✨ Key Features

*   **🚀 Blazing Fast**: With the same number of threads, `FastDup` is approximately **8X** faster than GATK MarkDuplicatesSpark and **20X** faster than Picard MarkDuplicates.
*   **✅ Identical Results**: Generates outputs that are identical to those of Picard MarkDuplicates.
*   **📊 Detailed Metrics**: Provides the same detailed metrics data as Picard MarkDuplicates.
*   **🧠 Memory Efficient**: All data is processed in memory, maintaining a low memory footprint even with very large input files.

## ⚠️ Limitations

*   **Marking Stability**: Although `FastDup` detects the exact same set of duplicates as Picard MarkDuplicates, it may mark a different read as the "duplicate" within a group. This is because the sorting algorithm in Picard MarkDuplicates is unstable.
    *   **For example**: Consider a duplicate group with two reads, A and B, having the same score, where A appears before B in the file. Picard might reorder them, placing B before A, and consequently mark A as the duplicate. In contrast, `FastDup` uses a stable sort algorithm, which preserves the original relative order and thus always marks B as the duplicate.
*   **Data Overflow in Optical Duplicate Detection**: Picard MarkDuplicates uses the `short` (int16_t) data type when parsing tile/region and x/y coordinates from a read name. This can lead to a data overflow if these integer values exceed the range of a `short`. `FastDup` identifies this bug but retains this behavior in the source code to maintain consistency with Picard.
    *   **To fix this**: You can resolve this issue by changing the relevant data types in the `PhysicalLocation` struct within the `read_ends.h` file.
*   **Input File Requirement**: `FastDup`'s performance improvements rely on the data characteristics of coordinate-sorted files. Therefore, the **input SAM/BAM file must be sorted by coordinate beforehand**.

## 🛠️ Requirements

Before you begin, ensure you have the following tools and libraries installed.

```bash
# Install autoconf (for htslib), cmake, a C++17 compiler (GCC >= 8.1 or Clang >= 7 should work),
# zlib, libbz2, liblzma, libcurl, and libdeflate (optional).
sudo apt update
sudo apt install -y autoconf cmake g++-8 zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libdeflate-dev
```

## 🚀 Installation

You can clone the source code directly from GitHub.

```bash
# 1. Clone the repository
git clone https://github.com/zzhofict/FastDup.git
cd FastDup

# 2. Build the bundled htslib
cd ext/htslib
autoreconf -i
./configure
make
cd ../..

# 3. Build FastDup
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

## 💡 Usage

### Get help
```bash
# Navigate to the project root
cd FastDup

# Run the executable with the --help flag
./build/bin/fastdup --help
```

### Example Command
Mark duplicates on an input BAM file using 8 threads and generates a metrics file.
```bash
# Navigate to the project root
cd FastDup

# Run the command
./build/bin/fastdup \
    --input ./test/input/in_test.bam \
    --output ./test/output/out_md.bam \
    --metrics stats.txt \
    --num-threads 8
```
