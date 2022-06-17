# Mutation profile

This repository comprises the script(s) developed during Monkeypox 2022 outbreak to explore the mutational profiles/signatures of this virus, but that can be of broad application to other species. Currently, it comprises the script(s):
- _get_mutation_profile.py_ that can be used to rapidly obtain the sequence context (size defined by the user) flanking SNPs of interest and determine their mutational profile according to the user's specifications (e.g. APOBEC3-mediated viral genome editing GA>AA and TC>TT replacements)

## Input/Output of _get_mutation_profile.py_

**OPTION 1**  
_Inputs:_    
1. Single-column file with a list of 1-indexed reference positions of interest   
2. Multiple Sequence Alignment (fasta) including the reference genome    
									
_Outputs:_     
1. TSV file with the mutation context and profile for each sample present in the alignment
2. TSV file with a summary report for each position of interest including the different patterns observed and their respective frequency
									
**OPTION2**   
_Inputs:_    
1. TSV file with the columns POS REF ALT (i.e. 1-indexed reference position, reference allele and alternative allele)
2. Fasta file including the reference genome

_Output:_    
1. TSV file with the mutation context and profile

**OPTION 3**     
_Inputs:_   
1. TSV file with the columns ID POS REF ALT (i.e. sample ID, 1-indexed reference position, reference allele and alternative allele)
2. Fasta file including the reference genome

_Outputs:_    
1. TSV file with the mutation context and profile for each sample present in the TSV input
2. TSV file with a summary report for each position of interest including the different patterns observed and their respective frequency


_NOTE: For options 2 and 3 the order of the columns in the input 1 is not important but their name is (ID, POS, REF, ALT)!!_


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

## Examples using Monkeypox 2022 outbreak data available in [_examples/_](https://github.com/insapathogenomics/mutation_profile/tree/main/examples)

### Option 1
Providing a single-column file with a list of 1-indexed reference positions of interest and a fasta Multiple Sequence Alignment including the reference genome, 

_TIP:_ If you do not know your positions of interest, you can run the script _alignment_processing.py_ of [ReporTree](https://github.com/insapathogenomics/ReporTree) and it will provide you a list of positions of interest according to your specifications.

```bash
python get_mutation_profile.py -f alignment_Figure1B.fasta -m positions_of_interest.tsv -r 'MT903344.1_Monkeypox_virus_isolate_MPXVUK_P2_complete_genome' -b 10 -a 10 -o OPTION1
```

_Outputs:_     
1. TSV file with the mutation context and profile for each sample present in the alignment
<p align="center">
	<img width="926" alt="Captura de ecrã 2022-06-17, às 12 53 35" src="https://user-images.githubusercontent.com/19263468/174293274-6a24067d-7453-433d-a8cf-c30e72ba33e7.png">
</p>

2. TSV file with a summary report for each position of interest including the different patterns observed and their respective frequency
<p align="center">
	<img width="1316" alt="Captura de ecrã 2022-06-17, às 12 56 12" src="https://user-images.githubusercontent.com/19263468/174293631-3ef7a1db-1117-4908-9a98-bca98f692aef.png">
</p>

### Option 2
Providing a TSV file with the columns POS REF ALT (i.e. 1-indexed reference position, reference allele and alternative allele) and a fasta file including the reference genome (can be the same alignment or a normal fasta sequence).

```bash
python get_mutation_profile.py -f alignment_Figure1B.fasta -m positions_of_interest_with_REF_ALT_info.tsv -r 'MT903344.1_Monkeypox_virus_isolate_MPXVUK_P2_complete_genome' -b 10 -a 10 -o OPTION2
```

_Output:_    
1. TSV file with the mutation context and profile
<p align="center">
	<img width="697" alt="Captura de ecrã 2022-06-17, às 12 59 20" src="https://user-images.githubusercontent.com/19263468/174294057-354bca2d-8339-4546-8acd-f26625f519ae.png">
</p>

### Option 3
Providing a TSV file with the columns ID POS REF ALT (i.e. samples id, 1-indexed reference position, reference allele and alternative allele) and a fasta file including the reference genome (can be the same alignment or a normal fasta sequence).

```bash
python get_mutation_profile.py -f alignment_Figure1B.fasta -m positions_of_interest_with_ID_REF_ALT_info.tsv -r 'MT903344.1_Monkeypox_virus_isolate_MPXVUK_P2_complete_genome' -b 10 -a 10 -o OPTION3
```

_Outputs:_    
1. TSV file with the mutation context and profile for each sample present in the TSV input
<p align="center">
	<img width="733" alt="Captura de ecrã 2022-06-17, às 13 02 18" src="https://user-images.githubusercontent.com/19263468/174294457-10795e22-07b3-4b47-8572-57c381a0100f.png">
</p>

2. TSV file with a summary report for each position of interest including the different patterns observed and their respective frequency
<p align="center">
	<img width="1148" alt="Captura de ecrã 2022-06-17, às 13 03 43" src="https://user-images.githubusercontent.com/19263468/174294624-3dc29037-ec77-41e6-a644-135125699b76.png">
</p>

## Citation

If you use this script please cite this github page.
