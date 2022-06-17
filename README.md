# Mutation profile

This repository comprises the script(s) developed during Monkeypox 2022 outbreak to explore the mutational profiles/signatures of this virus, but that can be of broad application to other species. Currently, it comprises the script(s):
- _get_mutation_profile.py_ that can be used to rapidly obtain the sequence context (size defined by the user) flanking SNPs of interest and determine their mutational profile according to the user's specifications (e.g. APOBEC3-mediated viral genome editing GA>AA and TC>TT replacements)

## Input/Output of _get_mutation_profile.py_

**OPTION 1**
Inputs:
- Single-column file with a list of 1-indexed reference positions of interest
- Multiple Sequence Alignment (fasta) including the reference genome
									
Outputs:
- TSV file with the mutation context and profile for each sample present in the alignment
- TSV file with a summary report for each position of interest including the different patterns observed and their respective frequency
									
**OPTION2**
Inputs:
- TSV file with the columns POS REF ALT (i.e. 1-indexed reference position, reference allele and alternative allele)
- Fasta file including the reference genome

Output:
- TSV file with the mutation context and profile

**OPTION 3**
Inputs:
- TSV file with the columns ID POS REF ALT (i.e. sample ID, 1-indexed reference position, reference allele and alternative allele)
- Fasta file including the reference genome

Outputs:
- TSV file with the mutation context and profile for each sample present in the TSV input
- TSV file with a summary report for each position of interest including the different patterns observed and their respective frequency

_NOTE: For options 2 and 3 the order of the columns in the input 1 is not important but their name is (ID, POS, REF, ALT)!!_
