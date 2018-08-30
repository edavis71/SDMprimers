# Automated site directed mutagenesis primers with validation

This program was created to automate the design of site directed mutagenesis (SDM) primers. This is particularly handy when designing these in bulk for multiple variants. For each variant, multiple primer options are returned along with their melting temperature and GC content. Additionally, a validation step is applied that cross references the HGVS change to the protein annotation and flags any mismatches. *The current functionality is for missense, insertion, and deletions made in cDNA and is applicable for any gene.*


### Prerequisites

This program was implemented for use with python version HERE

A text file containing the annotation for the protein change and HGVS nomenclature for each variant you want to design primers for is required as input along with a text file containing the cDNA sequence for the respective gene. 

```
Protein | HGVS 
--- | --- | ---
Y569D | c.1705T>G
L138ins | c.411_412insCTA
A1004_A1006del | c.3009_3017delAGCTATAGC

```

### Running 

Run the program as shown below

```
`./makeprimers.py ../test.txt ../cdna.txt`
```

### Output

Two files are returned, one containing any mismatches between protein annotation and HGVS change and one containing a list of forward and reverse primer pairs for each variant, along with their estimated melting temperature and GC content. Primers are formatted so that once you've chosen which to order based on your needs you can copy and paste the information directly into idt for bulk ordering.


Protein | HGVS 
--- | --- | ---
Y569D | c.1705T>G
L138ins | c.411_412insCTA
A1004_A1006del | c.3009_3017delAGCTATAGC


### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc

