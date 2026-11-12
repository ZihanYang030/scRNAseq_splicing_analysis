# scRNAseq Alternative Splicing Analysis

single-cell RNA-seq workflow to perform alternative splicing analysis using the public [GSE185862 dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185862) 

What it does:

- retrieves SRA reads
- creates FASTQ files
- trims reads
- aligns with STAR
- performs splicing analysis for [*Pakap*](https://www.informatics.jax.org/marker/MGI:5141924) gene and terminal exon-6 evidence in L2/3 IT CTX cells and oligodendrocytes

note: exact target coordinates and working data / file paths are intentionally **de-identified**

All scripts here use one single project root, `PROJECT_DIRECTORY`. Replace it as needed

## Requirements

- [SRA Toolkit](https://github.com/ncbi/sra-tools) (`prefetch`, `fasterq-dump`)
- [Trimmomatic](https://github.com/usadellab/Trimmomatic)
- [STAR](https://github.com/alexdobin/STAR) and a compatible mouse reference-genome index
- `samtools`, Python 3, and R for the targeted analysis
- Optional: `pigz` for faster compression; AWS CLI (not available for this dataset, provided for reference)

Main analysis pipeline: SRA Toolkit 3.2.0, Trimmomatic 0.39, and STAR 2.7.11b. \
Optional web scraper needs: `pandas`, `playwright`, and a Playwright browser installation

## Inputs

If data is stored separately from the scripts, set `PROJECT_DIRECTORY` to a local data directory outside the repository. Also, only one accession should be placed in each line in `SRR_Acc_List.txt`. 

The target-specific scripts read `TARGET_REGIONS_FILE`, see [`config/target_regions.example.tsv`](config/target_regions.example.tsv) for example. Replace the placeholders with your own annotation values (please use one-based inclusive coordinates) and store the filled file outside this repository. Make sure to keep genome build and chromosome names consistent across the annotation, STAR index, and BAM files.

```bash
export PROJECT_DIRECTORY=/path/to/local_analysis_data
export TARGET_REGIONS_FILE=/path/to/target_regions.tsv
python3 scripts/regions.py validate
```

The commands below assume they are run from the repository root. \
**Please make sure to check available storage before downloading or aligning a large dataset.**

## Workflow

1. **Prepare accession lists** \
   Download the public GSE185862 metadata (which is large, check your storage before download), confirm the run-to-cell mapping, and write the three accession lists. \
   Each SRR run is treated as a cell in this analysis (this might be different for different datasets). Please review the source metadata before reusing this assumption. 

2. **Download and convert reads** \
   Run `bash scripts/prefetch_missing_sra.sh "$PROJECT_DIRECTORY/SRR_Acc_List.txt"`, then `bash scripts/prefetch_fasterq_generic.sh "$PROJECT_DIRECTORY/SRR_Acc_List.txt"`. \
   Use `scripts/check_missing_sra.sh` and `scripts/list_missing_fastq.sh` to find incomplete inputs. 

3. **Prepare each population and align** \
   Use `scripts/link_fastq_by_run.sh` to create run-specific FASTQ links (optional). \
   Build a compatible STAR index in `$PROJECT_DIRECTORY/star_index` using the reference genome. \
   Then set `RUN_NAME` and `LIST` and run `scripts/run_trim_star.sh` for each population. \
   The trimming settings used for this project: `HEADCROP:10 TRAILING:15 AVGQUAL:20 MINLEN:30`. \
   Use `scripts/summarize_star_mapping.sh` summarizes the mapping quality 

4. **Define target exons and collect junctions** \
   Run `python3 scripts/write_manual_Pakap_exons_current.py` to create annotation tables from the coordinates. \
   Run `python3 scripts/find_Pakap_junctions_current.py` to collect unique-read junction evidence from STAR output. 

5. **Calculate EIP and targeted junction measures** \
   Run `Rscript scripts/calc_Pakap_manual_EIP_current.R`. calculates the regular internal-exon EIP and exon-6 junction measures for each cell and population (as exon-6 appears as the terminal exon, we then switched to classifying evidences for entry, internal continuation, and skipping)

6. **Measure terminal-region coverage** \
   Run `bash scripts/combine_accession_lists.sh`, `bash scripts/run_star_region_bam_L23ITCTX_Oligo_all.sh`, and `bash scripts/count_Pakap_exon6_region_coverage_all.sh` \
   The BAM script temporarily creates a full sorted BAM, retains the configured `bam_window`, and removes the temporary full BAM. \
   This step expects paired-end trimmed reads.

7. **Classify evidence and chekc for overlap** \
   Run `Rscript scripts/calc_Pakap_terminal_exon6_evidence_all_cells.R` and `Rscript scripts/check_category_overlap.R`. \
   Then check `$PROJECT_DIRECTORY/results` for results.

Example commands for the two population-specific alignment runs:

```bash
bash scripts/link_fastq_by_run.sh L23ITCTX_matched "$PROJECT_DIRECTORY/L23ITCTX_Acc_List.txt"
RUN_NAME=L23ITCTX_matched LIST="$PROJECT_DIRECTORY/L23ITCTX_Acc_List.txt" bash scripts/run_trim_star.sh

bash scripts/link_fastq_by_run.sh Oligo_matched "$PROJECT_DIRECTORY/Oligo_Acc_List.txt"
RUN_NAME=Oligo_matched LIST="$PROJECT_DIRECTORY/Oligo_Acc_List.txt" bash scripts/run_trim_star.sh
```

## Exon Inclusiong Proportion

For an internal exon, exon inclusion proportion (EIP) uses junction reads that support entry and exit from the exon (`I`) and reads supporting exon skipping (`E`):

$$EIP = \frac{I}{I + 2E}$$

For *Pakap* exon 6, `E5 → E6a` supports entry from exon 5 to exon 6 part a, `E6a → E7` supports internal continuation from exon 6 to exon 7, and `E5 → E7` supports skipping exon 6. 

**As entry alone cannot establish transcript termination,** we then used the targeted analysis to combine junction evidence with read overlap in exon 6b and the 3′ UTR. We define strong, moderate, or weak terminal candidates, internal-continuation evidence, and skip evidence at **cell level** (a cell can support more than one category through different reads).

## Scripts

| Steps                                | Files |
|--------------------------------------|---|
| Input checks and retrieval           | [`check_missing_sra.sh`](scripts/check_missing_sra.sh), [`prefetch_missing_sra.sh`](scripts/prefetch_missing_sra.sh), [`prefetch_fasterq_generic.sh`](scripts/prefetch_fasterq_generic.sh), [`list_missing_fastq.sh`](scripts/list_missing_fastq.sh) |
| (Optional) submitted-name lookup     | [`list_original_fastq_s3_for_srrs.sh`](scripts/list_original_fastq_s3_for_srrs.sh), [`scrape_original_fastq_names.py`](scripts/scrape_original_fastq_names.py) |
| FASTQ alignment                      | [`link_fastq_by_run.sh`](scripts/link_fastq_by_run.sh), [`run_trim_star.sh`](scripts/run_trim_star.sh), [`summarize_star_mapping.sh`](scripts/summarize_star_mapping.sh) |
| Junction and EIP analysis            | [`write_manual_Pakap_exons_current.py`](scripts/write_manual_Pakap_exons_current.py), [`find_Pakap_junctions_current.py`](scripts/find_Pakap_junctions_current.py), [`calc_Pakap_manual_EIP_current.R`](scripts/calc_Pakap_manual_EIP_current.R) |
| Regional coverage and classification | [`combine_accession_lists.sh`](scripts/combine_accession_lists.sh), [`run_star_region_bam_L23ITCTX_Oligo_all.sh`](scripts/run_star_region_bam_L23ITCTX_Oligo_all.sh), [`count_Pakap_exon6_region_coverage_all.sh`](scripts/count_Pakap_exon6_region_coverage_all.sh), [`calc_Pakap_terminal_exon6_evidence_all_cells.R`](scripts/calc_Pakap_terminal_exon6_evidence_all_cells.R), [`check_category_overlap.R`](scripts/check_category_overlap.R) |
| reader                               | [`regions.py`](scripts/regions.py) |

The optional AWS and browser lookup scripts are kept for reference, although it does not work for this specific dataset. When downloading data, be aware of its structure, and do not treat an empty lookup result as proof that a run lacks data.

Suggestions and feedbacks are greatly appreciated! Please contact: [zy56@rice.edu](mailto:zy56@rice.edu). Thank you!

## References

- [GSE185862, NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185862): source sequencing dataset and metadata
- [*Pakap* Gene, Mouse Genome Informatics](https://www.informatics.jax.org/marker/MGI:5141924): gene of interest
- [King et al., *Journal of Neuroscience* (2026)](https://www.jneurosci.org/content/46/16/e0019262026): alternative-splicing analysis approach
- [STAR](https://github.com/alexdobin/STAR), [Trimmomatic](https://github.com/usadellab/Trimmomatic), and [SRA Toolkit](https://github.com/ncbi/sra-tools): softwares used


## Citation

Xiaoyun Ding*, Yuki Ogawa, Victoria L. Palfini, Wei Zhang, Andrew P. Anderson, Lucia Gancedo-López, Zian Liao, Zihan Yang, Seth G. Haddix, Yueheng Xing, Duc V. M. Nguyen, Jackson R. Curtis, Davis J. Church, Yu Wu, Antrix Jain, Alexander B. Saltzman, Juan A. Oses-Prieto, Matthew Torre, Trinh T. Tat, Daniel C. Kraushaar, Alma L. Burlingame, Anna Malovannaya, Lulu Shang, Pengfei Liu, Jiaxing Li, Yudong Gao, Matthew N. Rasband*. ***Paralemmin 2 is a shared, alternatively-spliced regulator of axon initial segments and paranodal junctions in oligodendrocytes.*** (2026). [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.08.18.745582v1.abstract).
