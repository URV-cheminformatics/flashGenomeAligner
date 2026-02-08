# flashGenomeAligner 🧬

**flashGenomeAligner** is an ultra-efficient, parallel C++ tool designed for rapid gene-by-gene alignment and comprehensive mutation analysis of viral genomes. Specifically optimized for Zika and Dengue viruses, it utilizes SIMD instructions (via the Parasail library) and OpenMP multithreading to achieve high-performance processing.

## 🚀 Features

*   **Ultra-Fast Alignment:** Leverages Parasail for efficient Smith-Waterman/Needleman-Wunsch alignment with SIMD support.
*   **Parallel Processing:** Uses OpenMP to utilize all available CPU cores.
*   **Detailed Mutation Analysis:**
    *   Detects SNPs, Insertions, and Deletions.
    *   Translates codons to identify Missense, Nonsense, and Synonymous mutations.
    *   Tracks Frameshifts and intergenic regions.
*   **Structured Output:** Generates comprehensive TSV reports for downstream analysis.

## 📂 Project Structure

```text
flashGenomeAligner/
├── src/                # C++ Source code
├── external/           # External libraries (Parasail)
├── data/               # Example data and references
```

## 🛠 Prerequisites

*   **C++ Compiler:** with C++17 support (GCC, Clang, MSVC).
*   **CMake:** Version 3.10 or higher.
*   **OpenMP:** Usually included with GCC/Clang (install `libomp-dev` on Linux/macOS).
*   **Parasail:** The Parasail alignment library must be available in `external/parasail`.

## 📦 Compilation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/URV-cheminformatics/flashGenomeAligner.git
    cd flashGenomeAligner
    ```

2.  **Ensure Parasail is built/available:**
    This project uses Parasail for fast SIMD-optimized sequence alignment. The repository does not vendor Parasail by default, so you must clone and build it locally under `external/parasail`.

    ```bash
    # From the project root
    git clone https://github.com/jeffdaily/parasail.git external/parasail

    # Build Parasail (Release)
    cmake -S external/parasail -B external/parasail/build -DCMAKE_BUILD_TYPE=Release
    cmake --build external/parasail/build -j
    ```
    On Linux/macOS the build will produce `external/parasail/build/libparasail.so` (or `.a`).
    
    On Windows it will produce `parasail.lib` and/or `parasail.dll` depending on the generator/toolchain.

4.  **Build with CMake (Recommended for Performance):**
    
    It is crucial to build in **Release** mode to enable compiler optimizations (O3/AVX2) and SIMD instructions. Debug builds will be significantly slower.

    **Linux / macOS (single-config generators: Makefiles, Ninja)**

    ```bash
    # Configure (Release)
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

    # Build
    cmake --build build -j
    ```

    **Windows (Visual Studio / multi-config generators)**
    
    ```bash
    # Configure (Visual Studio generator creates multi-config builds)
    cmake -S . -B build

    # Build Release
    cmake --build build --config Release
    ```

## � Workflow: How to Analyze Your Data

To analyze your own viral sequences (e.g., Dengue, Zika, SARS-CoV-2), you need three input files:

1.  **Reference Genome (`.fasta` / `.fna`)**: One FASTA file containing the reference sequence.
2.  **Gene Definitions (`.tsv`)**: A tab-separated file defining the start/end positions of each gene in the reference.
    *   Format: `GeneName` `Start` `End` `Description`
3.  **Query Sequences (`.fasta`)**: The FASTA file containing all the viral sequences you want to analyze.

### Running the Analysis

Once you have your files, run the tool from the command line:

```bash
# General Syntax
./flashGenomeAligner <genes_layout.tsv> <reference.fasta> <query.fasta> [output_results.tsv]

# Example: Analyzing Zika Virus
./flashGenomeAligner data/Zika/genes.tsv data/Zika/ref.fasta data/Zika/samples.fasta results/zika_analysis.tsv
```

### Running the Default Test (SARS-CoV-2)
If you run the program without arguments, it will execute the built-in SARS-CoV-2 pilot test located in `data/examples/tests/SARS_CoV2_test/`:

```bash
./flashGenomeAligner
```

## 📝 License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

This project uses the **Parasail** library for high-performance sequence alignment:
*   Daily, Jeff. (2016). Parasail: SIMD C library for global, semi-global, and local pairwise sequence alignments. *BMC Bioinformatics*, 17(1), 1-11. [doi:10.1186/s12859-016-0930-z](https://doi.org/10.1186/s12859-016-0930-z)

## 👥 Author

**Diego Arcos**
*   Email: **diegoarcos33@gmail.com**


*   **Final Degree Project in Biotechnology**
*   **Universitat Rovira i Virgili (URV)** - Tarragona, Spain
*   **Cheminformatics and Nutrition Research Group**
*   Supervisor: **Santiago Garcia Vallvé**

---
*Created as part of the Final Degree Project (TFG) - 2025.*
