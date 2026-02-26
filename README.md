# MMseqs-based Iterative Homology Search with Taxonomic Annotation

## Overview

This repository provides a reproducible pipeline to perform iterative distant homology searches using MMseqs2 as a high-performance substitute for PSI-BLAST, followed by automated taxonomic annotation of detected hits.

The script is designed for scalable protein homology searches against large reference databases (e.g., RefSeq-scale protein datasets) and integrates downstream taxonomy annotation using NCBI accession-to-TaxID mappings and taxonkit.

## Rationale
Why MMseqs instead of PSI-BLAST?

PSI-BLAST (Position-Specific Iterated BLAST) has historically been the standard method for detecting distant protein homologs through iterative profile refinement.

However, PSI-BLAST becomes computationally expensive and slow when applied to modern large-scale protein databases (e.g., RefSeq, nr).

MMseqs2 (Many-against-Many sequence searching) provides:

  - Orders-of-magnitude speed improvements

  - Efficient iterative profile-based searches

  - Scalable performance for very large databases

  - High sensitivity comparable to PSI-BLAST

For large-scale distant homology searches, MMseqs2 provides a computationally tractable alternative to PSI-BLAST while maintaining high sensitivity.

## References
PSI-BLAST

Altschul SF, Madden TL, Schäffer AA, et al.
Gapped BLAST and PSI-BLAST: a new generation of protein database search programs.
Nucleic Acids Research (1997) 25(17):3389–3402.
https://doi.org/10.1093/nar/25.17.3389

## MMseqs2

Steinegger M, Söding J.
MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets.
Nature Biotechnology (2017) 35:1026–1028.
https://doi.org/10.1038/nbt.3988

Steinegger M, Söding J.
MMseqs2: sensitive protein sequence searching for the analysis of massive data sets.
Bioinformatics (2018) 34(15): 2659–2661.

## Database Preparation (Required Before Running the Script)

Before running the pipeline, the target protein database must be formatted as an MMseqs database.

1. Prepare Protein FASTA Database

  The reference protein database must be in FASTA format.

  Example (RefSeq protein):

    refseq_protein.faa

2. Create MMseqs Database

Convert FASTA into an MMseqs database:

    mmseqs createdb refseq_protein.faa MM_refseq_protein --dbtype 1 --threads 108

--dbtype 1 → protein database

--threads → adjust to your available CPU cores

This step creates the internal MMseqs database files.

3. (Optional but Recommended) Create Prefilter Index

For large databases (RefSeq-scale) and repeated searches, building a prefilter index significantly accelerates subsequent searches.

    mkdir -p /thirdisk/mmseqs_tmp

    mmseqs createindex MM_refseq_protein /thirdisk/mmseqs_tmp \
      --threads 60 \
      --check-compatible 1 \
      -s 7.5
### Important Considerations

  Indexing a RefSeq-scale database can require very large disk space.

  In our tests, index creation consumed >1 TB of temporary disk space.

  Disk usage scales with database size and sensitivity parameter (-s).

  Attempting to control temporary files:

    mkdir -p mmseqs_tmp

    mmseqs createindex MM_refseq_protein mmseqs_tmp/ \
      --threads 60 \
      --remove-tmp-files 1 \
      --check-compatible 1 \
      -s 7.5

  However:

    --remove-tmp-files 1 does not reduce peak disk usage during index construction.

  It only removes temporary files after completion.

  Practical Implication

  -For small gene sets (e.g., a few query genes), building a full RefSeq index may be:

  - Disk-expensive (>1 TB)

  - Potentially unnecessary

In such cases, direct searches without index precomputation may be more practical.

### Taxonomic Annotation Requirements

The script requires:

  - NCBI protein accession → TaxID mapping file
  (e.g., prot.accession2taxid or .gz)

  - NCBI taxonomy dump directory for taxonkit
  (containing nodes.dmp, names.dmp, etc.)

These are used to:

  Map MMseqs hits to TaxIDs
  Retrieve full lineage using taxonkit
  Summary of Pipeline
  Create MMseqs query database
  Perform iterative MMseqs search (PSI-BLAST equivalent strategy)
  Export alignment results  
  Map protein accessions to TaxIDs
  Retrieve taxonomic lineage
  Append taxonomy columns to final result table

### Final output:

    PREFIX_Res.with_taxonomy.tsv

Containing:

  - Original MMseqs alignment metrics
  - TaxID
  - Full taxonomic lineage

## Running the Pipeline

Basic Execution Example

After preparing your MMseqs target database (see Database Preparation section), you can run the pipeline as follows:

    chmod +x mmseqs_sting_taxonomy.sh

    ./mmseqs_sting_taxonomy.sh \
      -p sting_MM \
      -q sting.faa \
      -T MM_refseq_protein \
      -t 60 \
      -n 3 \
      -e 1e-5 \
      -m TaxonomyDB/prot.accession2taxid.gz \
      -d TaxonomyDB/ \
      -o results_sting
## Parameter Explanation

Parameter	Description

    -p	Prefix used to standardize output file names
    -q	Protein FASTA file containing one or multiple query sequences
    -T	Target MMseqs protein database (must be preformatted)
    -t	Number of CPU threads
    -n	Number of MMseqs iterations (PSI-BLAST equivalent behavior)
    -e	E-value cutoff
    -m	Path to prot.accession2taxid mapping file (plain or gzipped)
    -d	Path to NCBI taxonomy dump directory for taxonkit
    -o	Output directory

Optional flags:

Flag	Description

    --keep-tmp	Keep MMseqs temporary directory
    --keep-tax	Keep intermediate taxonomy files
Output Files

The pipeline produces:

    results_sting/
    ├── sting_MM_Res.tsv
    └── sting_MM_Res.with_taxonomy.tsv
    sting_MM_Res.tsv

Raw MMseqs alignment output including:

    Query accession
    Target accession
    Alignment coordinates
    Coverage
    Percent identity
    E-value
    Bit score
    sting_MM_Res.with_taxonomy.tsv
    Same table as above, with two additional columns appended:
    taxid
    lineage

The lineage column contains the full NCBI taxonomy path, e.g.:

    cellular organisms;Eukaryota;Metazoa;Arthropoda;Insecta;Diptera;Drosophilidae;Drosophila;Drosophila melanogaster

### Minimal Example (Small Test Run)

For testing with a small query file and reduced resources:

    ./mmseqs_sting_taxonomy.sh \
      -p test_run \
      -q small_queries.faa \
      -T MM_refseq_protein \
      -t 8 \
      -n 2 \
      -e 1e-3 \
      -m prot.accession2taxid.gz \
      -d taxonomy_db \
      -o test_output

This is recommended to validate:

  - Database formatting
  - Mapping file compatibility
  - Taxonomy directory correctness
  = Available disk space

### Recommended Strategy

For exploratory distant homology searches:

Start with:

    -n 3
    -s 7.5
    -e 1e-5

Increase sensitivity (-s) or iterations (-n) for deeper searches.

For very large databases:

  - Ensure sufficient disk space.
  - Consider whether building a full MMseqs index is justified.
  - Monitor temporary directory size.
