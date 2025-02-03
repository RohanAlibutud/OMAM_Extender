# Alignment Extender

This is a script designed to add new target species to the MSAs from OrthoMaM by retrieving homologous sequences from NCBI using BLAST.

- Retrieves homologous protein sequences for specified species from NCBI
- In the event of multiple sequences being returned, a consensus sequence is calculated
- Adds high-confidence BLAST hits to the alignments
- Outputs alignments in 2-line FASTA format

## Requirements

- Python 3.x
- Biopython

Install Biopython using pip:

```bash
pip install biopython
```

## Input

- **OrthoMaM Alignments Folder:** Contains the original MSAs.
- **Species List:** Comma-separated list of species names to be added.

## Usage

Run the script using the following command:

```bash
python alignment_extender.py --spp "species1,species2,species3" --omam /path/to/orthomam_alignments --out /path/to/output_folder
```

### Arguments

- `--spp` (required): Comma-separated list of species to add, can be however long
- `--omam` (required): Path to the folder containing OrthoMaM alignments
- `--out` (required): Path to the output folder where extended alignments will be saved

### Example

```bash
python alignment_extender.py --spp "Homo_sapiens,Pan_troglodytes" --omam ./alignments --out ./extended_alignments
```

## Script Overview

1. **Argument Parsing:** Handles command-line arguments.
2. **BLAST Search:** Uses NCBI BLAST to find homologous protein sequences.
3. **Consensus Calculation:** Derives a consensus sequence from the input alignment.
4. **Species Addition:** Adds new sequences based on BLAST hits.
5. **Output Generation:** Saves the updated alignments in FASTA format.


## Notes

- The consensus sequence thing is a kludge and is probably highly questionable tbh
- The script filters BLAST hits based on E-value (< 0.001) and sequence coverage (> 70%).

## License

This project is licensed under the MIT License.

## Acknowledgments

Built using Biopython and NCBI BLAST APIs.

---

*Updated: 2025-02-03*

