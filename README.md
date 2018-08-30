# Automated site directed mutagenesis primers with validation

This program was created to automate the design of site directed mutagenesis (SDM) primers. This is particularly handy when designing these in bulk for multiple variants. For each variant, multiple primer options are returned along with their melting temperature and GC content. Additionally, a validation step is applied that cross references the HGVS change to the protein annotation and flags any mismatches. *The current functionality is for missense, insertion, and deletions made in cDNA and is applicable for any gene.*


### Prerequisites

This program was implemented for use with python version HERE

A text file containing the annotation for the protein change and HGVS nomenclature for each variant you want to design primers for is required as input along with a text file containing the cDNA sequence for the respective gene. 


| Protein | HGVS |
| --- | --- |
| Y569D | c.1705T>G |
| L138ins | c.411_412insCTA |
| A1004_A1006del | c.3009_3017delAGCTATAGC |

### Running 

Run the program as shown below

```
`./makeprimers.py ../test.txt ../cdna.txt`
```

### Output

Two files are returned, one containing any mismatches between protein annotation and HGVS change and one containing a list of forward and reverse primer pairs for each variant, along with their estimated melting temperature and GC content. Primers are formatted so that once you've chosen which to order based on your needs you can copy and paste the information directly into idt for bulk ordering.


| Protein | HGVS | 
| --- | --- |
| Y569D | c.1705T>G |
| L138ins | c.411_412insCTA |
| A1004_A1006del | c.3009_3017delAGCTATAGC |

## Feedback

Please send feedback to Emily Davis
<edavis71@jhu.edu>


## Acknowledgments

This code was originally designed for simplifying and eliminating human error in designing primers used in creating mutant CFTR cell lines. I thank everyone in the Cutting lab for their support and inspiration. 

## References

