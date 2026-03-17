# geneslator-data

Geneslator Database Repository for geneslator R package (https://github.com/knowmics-lab/geneslator).

## General overview

**geneslator** is a comprehensive R package for gene identifier conversion and genome annotation across multiple model organisms. The package integrates data from several cross-organism databases and organism-specific resources within a single, coherent framework. Four different types of data about a gene are integrated: annotations from general databases (symbol, aliases, full name, genetype), annotations from species-specific databases, functional annotations (pathways and gene ontologies), and orthologs.

![Geneslator's workflow](percorso/della/tua/immagine.png)

Currently, annotation databases have been built for the following 8 model organisms: Human (Homo sapiens), Mouse (Mus musculus), Rat (Rattus norvegicus), Fly (Drosophila melanogaster ), Zebrafish (Danio rerio), Yeast (Saccharomyces cerevisiae), Worm (Caenorhabditis elegans), and Arabidopsis (Arabidopsis thaliana). More organisms will be included in future releases of **geneslator**.

## Data sources

General information about a gene (symbol, aliases, full name, and genetype) are extracted from NCBI Gene and Ensembl. Genetype represents the biotype classification of a gene (e.g., “protein-coding gene”, “non-coding RNA”, “pseudogene”, “lncRNA”). Databases for A.thaliana, C.elegans, D.melanogaster, and S.cerevisiae also include locus tag identifiers. 

Identifiers of a gene include Entrez GeneIDs (taken from NCBI), Ensembl GeneIDs (taken from NCBI and Ensembl), Uniprot IDs of its proteins (taken from Uniprot) and species-specific identifiers, coming from the most popular species-specific genome database, namely HGNC for Human, MGI for Mouse, RGD for Rat, SGD for Yeast, WormBase for Worm, FlyBase for Fly, ZFIN for Zebrafish and TAIR for Arabidopsis. For Zebrafish, we also collect Ensembl GeneID and Gene symbols data from HCOP. 

**geneslator** annotation databases also integrates old discontinued and replaced gene identifiers from NCBI gene and Ensembl (starting from v.28 for Arabidopsis and from v.81 in the other organisms). These archived identifiers are stored in different columns with respect to current identifiers.

Genes’ orthologs are taken from NCBI, Ensembl and AllianceGenome. For Human, we also collect data from HCOP. Orthologs are represented by their gene symbols.

Pathway data include pathway ids and their names and are collected from KEGG Pathways, Reactome and Wikipathways.

Gene ontology data are taken from GO and include GO IDs, full names, types (biological process, cellular component or molecular function) and evidence codes of gene annotations.

## Data integration

Integration of general information about genes and gene identifiers is done by prioritizing NCBI information over Ensembl data. For Zebrafish, integration of gene identifiers is done by giving the highest priority to NCBI, followed by HCOP and Ensembl.

Integration of orthologs data referring to the same gene has been done according to the following order: NCBI, HCOP (for Human), AllianceGenome and Ensembl. 

Annotation databases resulting from the integration of all gene are built as SQLite objects using the AnnotationForge R package.

## Database releases

**geneslator** annotation databases are updated each month. At each update, databases are collected in a new Github release. All published releases are available at [<https://github.com/knowmics-lab/geneslator-data/releases>](<https://github.com/knowmics-lab/geneslator-data/releases>).
