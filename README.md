# MolecularMimicry
The algorithm of search for human T-cell epitopes in proteomes of microbiota

### Installation
```
git clone https://github.com/vasilisa-kovaliova/combinatorial_peptide_pooling.git
```

### Requirements

To use it you will need to install Python packages:
- Bio
- rapidfuzz
- pandas

You can install them with pip:
```
pip install -r requirements.txt
```

Or separately with pip:
```
pip install biopython
```
```
pip install rapidfuzz
```
```
pip install pandas
```
### Usage

Parameters:

-- microbiota_names -- path to the file containing names of the .fasta files with microbial proteomes
-- microbiota_path -- path to the directory with microbial proteomes
-- epitopes -- path to the file with epitopes
-- Hamming_distance -- Hamming distance for search
-- results -- path to the file for the resulting table

## Example
```
python bacteria_proteome_scanner.py -microbiota_names microbiota_names.str -microbiota_path microbiota_proteomes/ -epitopes epitopes.txt -Hamming_distance 4 -results scanning_results.tsv
```
