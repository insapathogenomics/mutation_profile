# Mutation profile

This repository comprises the script(s) developed during Monkeypox 2022 outbreak to explore the mutational profiles/signatures of this virus, but that can be of broad application to other species. Currently, it comprises the script(s):
- _get_mutation_profile.py_ that can be used to rapidly obtain the sequence context (size defined by the user) flanking SNPs of interest and determine their mutational profile according to the user's specifications (e.g. APOBEC3-mediated viral genome editing GA>AA and TC>TT replacements)

## Input/Output of _get_mutation_profile.py_
									
**OPTION1**   
_Inputs:_    
1. TSV file with the columns POS REF ALT (i.e. 1-indexed reference position, reference allele and alternative allele)
2. Fasta file including the reference genome

_Output:_    
1. TSV file with the mutation context and profile

**OPTION 2**     
_Inputs:_   
1. TSV file with the columns ID POS REF ALT (i.e. sample ID, 1-indexed reference position, reference allele and alternative allele)
2. Fasta file including the reference genome

_Outputs:_    
1. TSV file with the mutation context and profile for each sample present in the TSV input
2. TSV file with a summary report for each position of interest including the different patterns observed and their respective frequency

_NOTE: For options 1 and 2 the order of the columns in the input 1 is not important but their name is (ID, POS, REF, ALT)!!_

**OPTION 3**  
_Inputs:_    
1. Single-column file with a list of 1-indexed reference positions of interest   
2. Multiple Sequence Alignment (fasta) including the reference genome    
									
_Outputs:_     
1. TSV file with the mutation context and profile for each sample present in the alignment
2. TSV file with a summary report for each position of interest including the different patterns observed and their respective frequency

_TIP: If you do not know your positions of interest, you can run the script _alignment_processing.py_ of [ReporTree](https://github.com/insapathogenomics/ReporTree) and it will provide a list of positions of interest according to your specifications._


## Dependencies and installation
To run the _get_mutation_profile.py_ you will need:
- biopython
- pandas

If you want to create a conda environment with these dependencies:
```bash
conda create -n mutation_profile -c conda-forge -c anaconda biopython pandas
```

To install and run:
```bash
git clone https://github.com/insapathogenomics/mutation_profile.git
cd mutation_profile/
conda activate mutation_profile # if you created the conda environment
python get_mutation_profile.py -h
```

## Usage

```bash
  -h, --help            show this help message and exit

Mutation profile:
  Provide input/output specifications

  -f FASTA, --fasta FASTA
                        [MANDATORY] Input sequence file (fasta)
  -m MUTATION, --mutation_list MUTATION
                        [MANDATORY] Input mutation list that can be: 1)
                        single-column file with 1-based reference position
                        information (in this case the fasta file must be a
                        multiple sequence alignment of all the sequences of
                        interest); OR 2) tsv file with the columns POS, REF,
                        and ALT where POS = 1-based reference position. If you
                        want to include information for more than one sample
                        per position, add also the column 'ID' (note that the
                        order of the columns is not important but their name
                        is!)
  -r REF, --reference REF
                        [MANDATORY] Reference sequence name
  -b BEFORE, --before BEFORE
                        [OPTIONAL] Number of nucleotides to report BEFORE the
                        mutation (default = 5)
  -a AFTER, --after AFTER
                        [OPTIONAL] Number of nucleotides to report AFTER the
                        mutation (default = 5)
  -p PROFILES, --profiles PROFILES
                        [OPTIONAL] Comma-separated list of mutational profiles
                        of interest (upper-case!). Default = 'GA>AA,TC>TT'
  -o OUTPUT, --output OUTPUT
                        [OPTIONAL] Tag for output file name. Default =
                        Mutation_profile
```

## Examples using Monkeypox 2022 outbreak data available at [_examples/_](https://github.com/insapathogenomics/mutation_profile/tree/main/examples)

### Option 1 (this option reflects part of the analysis performed in the publication)
Providing a TSV file with the columns POS REF ALT (i.e. 1-indexed reference position, reference allele and alternative allele) and a fasta file including the reference genome (can be the same alignment or a normal fasta sequence).

```bash
python get_mutation_profile.py -f alignment_Figure1B.fasta -m positions_of_interest_POS_REF_ALT.txt -r 'MT903344.1_Monkeypox_virus_isolate_MPXVUK_P2_complete_genome' -b 10 -a 10 -o OPTION1
```

_Output:_    
1. TSV file with the mutation context and profile
<p align="center">
	<img width="636" alt="Captura de ecrã 2022-06-17, às 15 17 41" src="https://user-images.githubusercontent.com/19263468/174316512-056bb280-a89e-4d81-9ba3-6065f6a9a0f1.png">
</p>

### Option 2
Providing a TSV file with the columns ID POS REF ALT (i.e. samples id, 1-indexed reference position, reference allele and alternative allele) and a fasta file including the reference genome (can be the same alignment or a normal fasta sequence).

```bash
python get_mutation_profile.py -f alignment_Figure1B.fasta -m positions_of_interest_ID_POS_REF_ALT.txt -r 'MT903344.1_Monkeypox_virus_isolate_MPXVUK_P2_complete_genome' -b 10 -a 10 -o OPTION2
```

_Outputs:_    
1. TSV file with the mutation context and profile for each sample present in the TSV input
<p align="center">
	<img width="849" alt="Captura de ecrã 2022-06-17, às 15 21 20" src="https://user-images.githubusercontent.com/19263468/174317025-e090f7f8-17a4-4001-863f-cd965dfff9a6.png">
</p>

2. TSV file with a summary report for each position of interest including the different patterns observed and their respective frequency
<p align="center">
	<img width="1145" alt="Captura de ecrã 2022-06-17, às 15 23 07" src="https://user-images.githubusercontent.com/19263468/174317375-d17183e8-0f9c-46a6-ab83-9f1e212c4854.png">
</p>

### Option 3
Providing a single-column file with a list of 1-indexed reference positions of interest and a fasta Multiple Sequence Alignment including the reference genome.

```bash
python get_mutation_profile.py -f alignment_Figure1B.fasta -m Monkeypox_positions_of_interest.tsv -r 'MT903344.1_Monkeypox_virus_isolate_MPXVUK_P2_complete_genome' -b 10 -a 10 -o OPTION3
```

_Outputs:_     
1. TSV file with the mutation context and profile for each sample present in the alignment
<p align="center">
	<img width="1106" alt="Captura de ecrã 2022-06-17, às 15 31 56" src="https://user-images.githubusercontent.com/19263468/174319330-2702cf97-af7d-44e9-b7ca-2bf2f0f612ed.png">
</p>

2. TSV file with a summary report for each position of interest including the different patterns observed and their respective frequency
<p align="center">
<img width="1359" alt="Captura de ecrã 2022-06-17, às 15 33 18" src="https://user-images.githubusercontent.com/19263468/174319340-d95fd034-a740-48f1-8cd5-53ea77e785e8.png">	
</p>

_TIP: If you do not know your positions of interest, you can run the script _alignment_processing.py_ of [ReporTree](https://github.com/insapathogenomics/ReporTree) and it will provide a list of positions of interest according to your specifications. Example:_

```bash
python ReporTree/scripts/alignment_processing.py -align alignment_Figure1B.fasta -o Monkeypox --use-reference-coords -r 'MT903344.1_Monkeypox_virus_isolate_MPXVUK_P2_complete_genome' --keep-gaps --get-positions-interest
```

## Citation

If you use this script please cite the article where it was first described:

Isidro, J., Borges, V., Pinto, M. et al. Phylogenomic characterization and signs of microevolution in the 2022 multi-country outbreak of monkeypox virus.  
Nature Medicine (2022). https://doi.org/10.1038/s41591-022-01907-y
