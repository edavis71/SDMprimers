# Automated site directed mutagenesis primers with validation

This program was created to automate the design of site directed mutagenesis (SDM) primers. This is particularly handy when designing these in bulk for multiple variants. For each variant, multiple primer options are returned along with their melting temperature and GC content. Additionally, a validation step is applied that cross references the HGVS change to the protein annotation and flags any mismatches. *The current functionality is for missense, insertion, and deletions made in cDNA and is applicable for any gene.*


### Prerequisites

This program was implemented for use with python version 2.7

A text file containing the annotation for the protein change and HGVS nomenclature for each variant you want to design primers for is required as input along with a text file containing the cDNA sequence for the respective gene. 


| Protein | HGVS |
| --- | --- |
| Y569D | c.1705T>G |
| L138ins | c.411_412insCTA |
| A1004_A1006del | c.3009_3017delAGCTATAGC |

### Running 

Run the program as shown below

```
./makeprimers.py /variants.txt /cdna.txt
```

### Output

Two files are returned, one containing any mismatches between protein annotation and HGVS change and one containing a list of forward and reverse primer pairs for each variant, along with their estimated melting temperature and GC content. Primers are formatted so that once you've chosen which pair to order based on your needs you can copy and paste the information directly into IDT for bulk ordering.

Y569D c.1705T>G

| length | primer sequence | Tm | GC | format for IDT |
| --- | --- | --- | --- | --- |
| 47 | 'GCAGTATACAAAGATGCTGATTTGGATTTATTAGACTCTCCTTTTGG' | 65.42 | 36.17 | Y569D SDM  F GCAGTATACAAAGATGCTGATTTGGATTTATTAGACTCTCCTTTTGG 25nm STD |
| 47 | 'CCAAAAGGAGAGTCTAATAAATCCAAATCAGCATCTTTGTATACTGC' | 65.42 | 36.17 | Y569D SDM  R CCAAAAGGAGAGTCTAATAAATCCAAATCAGCATCTTTGTATACTGC 25nm STD |


Example error: 

```
warning: found mismatch G971D c.2909G>A is actually G970D
```

## Feedback

Please send feedback to Emily Davis:
<edavis71@jhu.edu>


## Acknowledgments

This code was originally designed for simplifying and eliminating human error in designing primers used in creating variant CFTR cell lines. I thank everyone in the Cutting lab for their support and inspiration. 

## References

* [Functional Assays Are Essential for Interpretation of Missense Variants Associated with Variable Expressivity. Am J Hum Genet. 2018 Jun 7;102(6):1062-1077](https://www.sciencedirect.com/science/article/pii/S0002929718301356)



