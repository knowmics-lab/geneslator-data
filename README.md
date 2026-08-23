# geneslator-data

This repository contains all the scripts used to build annotation databases for geneslator R package [<https://github.com/knowmics-lab/geneslator>](<https://github.com/knowmics-lab/geneslator>). The scripts are used to download the most up-to-date raw annotation data from different sources, pre-process them and integrate them to build the final annotation tables for all species currently supported by **geneslator**. Final annotation tables will be returned as SQLite databases in `.sqlite` file format.

Sources used to download raw annotation data for all species are reported in the `sources.json` file. This file will be constantly updated, as long as new species will be supported by **geneslator**.

## Usage examples

To run the main script (`build_geneslator_data.R`), R (R>=4.6) is required. It can be installed from [<https://cran.r-project.org/>](<https://cran.r-project.org/>). The script requires several hours to complete.

```bash
#Run the script with default output folder ("Output") for final annotation tables.
Rscript build_geneslator_data.R

#Run the script with specified output folder for final annotation tables.
Rscript build_geneslator_data.R --folder AnnotFiles

#Run the script with specified output folder for final annotation tables
#Return the annotation tables as both SQLite databases and TSV text files.
Rscript build_geneslator_data.R --folder AnnotFiles --saveTxt
```

## How it works

**geneslator** is a comprehensive R package for gene identifier conversion and genome annotation across multiple model organisms. The package integrates data from several cross-organism databases and organism-specific resources within a single, coherent framework. Four different types of data about a gene are integrated: annotations from general databases (symbol, aliases, full name, genetype), annotations from species-specific databases, functional annotations (pathways and gene ontologies), and orthologs.

![Geneslator's workflow](https://github.com/user-attachments/assets/f6c741a6-0bed-4c3b-b04b-4620f2940b71)

Currently, annotation databases have been built for the following 19 model organisms: 

- *Homo sapiens* (Human)
- *Mus musculus* (Mouse)
- *Rattus norvegicus* (Rat)
- *Danio rerio* (Zebrafish)
- *Drosophila melanogaster* (Fly)
- *Caenorhabditis elegans* (Worm)
- *Saccharomyces cerevisiae* (Yeast)
- *Arabidopsis thaliana* (Arabidopsis)
- *Brassica oleracea* (Cabbage)
- *Brassica napus* (Rapeseed)
- *Solanum lycopersicum* (Tomato)
- *Vitis vinifera* (Grapevine)
- *Lupinus angustifolius* (Blue Lupin)
- *Phaseolus vulgaris* (Common Bean)
- *Macaca mulatta* (Macaque)
- *Apis mellifera* (Honey Bee)
- *Xenopus laevis* (African Clawed Frog)
- *Oryza sativa* (Rice)
- *Zea mays* (Maize)

More organisms will be included in future releases of **geneslator**.

### Data sources

General information about a gene (symbol, aliases, full name, and genetype) are extracted from NCBI Gene and Ensembl. Genetype represents the biotype classification of a gene (e.g., “protein-coding gene”, “non-coding RNA”, “pseudogene”, “lncRNA”). Databases for A.thaliana, C.elegans, D.melanogaster, and S.cerevisiae also include locus tag identifiers. 

Identifiers of a gene include Entrez GeneIDs (taken from NCBI), Ensembl GeneIDs (taken from NCBI and Ensembl), Uniprot IDs of its proteins (taken from Uniprot) and species-specific identifiers, coming from the most popular species-specific genome database, e.g HGNC for Human, MGI for Mouse, RGD for Rat, SGD for Yeast, WormBase for Worm, FlyBase for Fly, ZFIN for Zebrafish, TAIR for Arabidopsis, VGNC for Macaque and Xenbase for Frog. For Zebrafish, we also collect Ensembl GeneID and Gene symbols data from HCOP. 

**geneslator** annotation databases also integrates old discontinued and replaced gene identifiers from NCBI gene and Ensembl (starting from v.28 for Arabidopsis and from v.81 in the other organisms). These archived identifiers are stored in different columns with respect to current identifiers.

Genes’ orthologs are taken from NCBI, Ensembl and AllianceGenome. For Human, we also collect data from HCOP. Orthologs are represented by their gene symbols.

Pathway data include pathway ids and their names and are collected from Reactome and Wikipathways. Pathway data from KEGG can be retrieved on-the-fly using geneslator R package [<https://github.com/knowmics-lab/geneslator>](<https://github.com/knowmics-lab/geneslator>).

Gene ontology data are taken from GO and include GO IDs, full names, types (biological process, cellular component or molecular function) and evidence codes of gene annotations.

### Data integration

Integration of general information about genes and gene identifiers is done by prioritizing NCBI information over Ensembl data. For Zebrafish, integration of gene identifiers is done by giving the highest priority to NCBI, followed by HCOP and Ensembl.

Integration of orthologs data referring to the same gene has been done according to the following order: NCBI, HCOP (for Human), AllianceGenome and Ensembl. 

Annotation databases resulting from the integration of all gene are built as SQLite objects using the AnnotationForge R package.

### Database releases

**geneslator** annotation databases are stored as a Zenodo record and available at [<https://zenodo.org/records/20457977>](<https://zenodo.org/records/20457977>). Databases are updated on a monthly basis. At each update, annotation databases are stored in a new version of the Zenodo record.

## Citation

If you use geneslator in your work, please cite:

```r
citation("geneslator")
```

geneslator: an R package for comprehensive gene identifier conversion and annotation. Giulia Cavallaro, Giovanni Micale, Grete Francesca Privitera, Alfredo Pulvirenti, Stefano Forte, Salvatore Alaimo. bioRxiv 2026.03.30.714723; doi: https://doi.org/10.64898/2026.03.30.714723 

Micale G, Cavallaro G, Privitera GF (2026). geneslator: A Comprehensive Gene Identifier Conversion Tool. R package version 0.99.0. https://github.com/knowmics-lab/geneslator

## Authors

- **Giovanni Micale** - *Author and maintainer* - [ORCID](https://orcid.org/0000-0002-4953-026X)
- **Giulia Cavallaro** - *Author* - [ORCID](https://orcid.org/0009-0000-1212-8368)
- **Grete Francesca Privitera** - *Author* - [ORCID](https://orcid.org/0000-0003-1807-4780)

University of Catania

## Support

- **Issues**: https://github.com/knowmics-lab/geneslator/issues
- **Email**: giovanni.micale@unict.it

## References

- NCBI Gene: https://www.ncbi.nlm.nih.gov/gene
- Ensembl: https://www.ensembl.org
- UniProt: https://www.uniprot.org
- Gene Ontology: http://geneontology.org
- KEGG: https://www.kegg.jp
- Reactome: https://reactome.org
- WikiPathways: https://www.wikipathways.org
- Alliance of Genome Resources: https://www.alliancegenome.org
- AnnotationDbi: Pages H, Carlson M, Falcon S, Li N (2024). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor.
