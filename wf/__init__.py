from typing import List, Optional

from latch.resources.launch_plan import LaunchPlan
from latch.resources.workflow import workflow
from latch.types import metadata
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile

from wf.entrypoint import (
    CPrimerPosition,
    LibraryGenerationMethod,
    LineageTreeBuilder,
    ProcessingMode,
    SampleSheet,
    UMIPosition,
    initialize,
    nextflow_runtime,
)


@workflow(metadata._nextflow_metadata)
def nf_nf_core_airrflow(
    input: LatchFile,
    run_name: str,
    email: Optional[str],
    library_generation_method: LibraryGenerationMethod,
    race_linker: Optional[LatchFile],
    vprimers: Optional[LatchFile],
    cprimers: Optional[LatchFile],
    primer_revpr: bool,
    index_file: bool,
    adapter_fasta: Optional[LatchFile],
    trim_nextseq: bool,
    save_trimmed: bool,
    maskprimers_align: bool,
    assemblepairs_sequential: bool,
    align_cregion: bool,
    internal_cregion_sequences: Optional[str],
    save_databases: bool,
    reference_fasta: Optional[LatchFile],
    reference_igblast: Optional[LatchFile],
    fetch_imgt: bool,
    detect_contamination: bool,
    remove_chimeric: bool,
    lineage_trees: bool,
    skip_all_clones_report: bool,
    skip_report_threshold: bool,
    reference_10x: Optional[str],
    skip_report: bool,
    skip_multiqc: bool,
    multiqc_title: Optional[str],
    multiqc_methods_description: Optional[str],
    mode: ProcessingMode = ProcessingMode.FASTQ,
    miairr: Optional[
        str
    ] = "${projectDir}/assets/reveal/mapping_MiAIRR_BioSample_v1.3.1.tsv",
    vprimer_start: Optional[int] = 0,
    cprimer_start: Optional[int] = 0,
    cprimer_position: CPrimerPosition = CPrimerPosition.R1,
    umi_position: UMIPosition = UMIPosition.R1,
    umi_length: Optional[int] = -1,
    umi_start: Optional[int] = 0,
    trim_fastq: bool = True,
    clip_r1: Optional[int] = 0,
    clip_r2: Optional[int] = 0,
    three_prime_clip_r1: Optional[int] = 0,
    three_prime_clip_r2: Optional[int] = 0,
    filterseq_q: Optional[int] = 20,
    primer_consensus: Optional[float] = 0.6,
    primer_mask_mode: Optional[str] = "cut",
    buildconsensus_maxerror: Optional[float] = 0.1,
    buildconsensus_maxgap: Optional[float] = 0.5,
    cluster_sets: bool = True,
    primer_r1_maxerror: Optional[float] = 0.2,
    primer_r2_maxerror: Optional[float] = 0.2,
    primer_extract_len: Optional[int] = 0,
    primer_maxlen: Optional[int] = 50,
    cregion_maxlen: Optional[int] = 100,
    cregion_maxerror: Optional[float] = 0.3,
    cregion_mask_mode: Optional[str] = "tag",
    reassign: bool = True,
    productive_only: bool = True,
    isotype_column: Optional[str] = "c_call",
    collapseby: Optional[str] = "sample_id",
    clonal_threshold: Optional[str] = "auto",
    cloneby: Optional[str] = "subject_id",
    crossby: Optional[str] = "subject_id",
    lineage_tree_builder: LineageTreeBuilder = LineageTreeBuilder.RAXML,
    lineage_tree_exec: Optional[str] = "/usr/local/bin/raxml-ng",
    singlecell: Optional[str] = "single_cell",
    report_rmd: Optional[str] = "${projectDir}/assets/repertoire_comparison.Rmd",
    report_css: Optional[str] = "${projectDir}/assets/nf-core_style.css",
    report_logo: Optional[str] = "${projectDir}/assets/nf-core-airrflow_logo_light.png",
    report_logo_img: Optional[
        str
    ] = "${projectDir}/assets/nf-core-airrflow_logo_reports.png",
    outdir: LatchOutputDir = LatchOutputDir("latch:///Airrflow"),
) -> None:
    """
    nf-core/airrflow is a bioinformatics best-practice pipeline to analyze B-cell or T-cell repertoire sequencing data.

    <html>
    <p align="center">
    <img src="https://user-images.githubusercontent.com/31255434/182289305-4cc620e3-86ae-480f-9b61-6ca83283caa5.jpg" alt="Latch Verified" width="100">
    </p>

    <p align="center">
    <strong>
    Latch Verified
    </strong>
    </p>

    <p align="center">

    [![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.2642009-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.2642009)
    [![AIRR compliant](https://img.shields.io/static/v1?label=AIRR-C%20sw-tools%20v1&message=compliant&color=008AFF&labelColor=000000&style=plastic)](https://docs.airr-community.org/en/stable/swtools/airr_swtools_standard.html)

    ## Introduction

    **nf-core/airrflow** is a bioinformatics best-practice pipeline to analyze B-cell or T-cell repertoire sequencing data. It makes use of the [Immcantation](https://immcantation.readthedocs.io) toolset. The input data can be targeted amplicon bulk sequencing data of the V, D, J and C regions of the B/T-cell receptor with multiplex PCR or 5' RACE protocol, single-cell VDJ sequencing using the 10xGenomics libraries, or assembled reads (bulk or single-cell).

    <html>
    <p align="center">
    <img src="https://github.com/nf-core/airrflow/blob/4.1.0/docs/images/airrflow_workflow_overview.png" alt="airflow" width="100">
    </p>

    The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

    This workflow is hosted on Latch Workflows, using a native Nextflow integration, with a graphical interface for accessible analysis by scientists. There is also an integration with Latch Registry so that batched workflows can be launched from “graphical sample sheets” or tables associating raw sequencing files with metadata.

    On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/airrflow/results).

    ## Pipeline summary

    nf-core/airrflow allows the end-to-end processing of BCR and TCR bulk and single cell targeted sequencing data. Several protocols are supported, please see the [usage documentation](https://nf-co.re/airrflow/usage) for more details on the supported protocols. The pipeline has been certified as [AIRR compliant](https://docs.airr-community.org/en/stable/swtools/airr_swtools_compliant.html) by the AIRR community, which means that it is compatible with downstream analysis tools also supporting this format.

    <html>
    <p align="center">
    <img src="https://github.com/nf-core/airrflow/blob/4.1.0/docs/images/metro-map-airrflow.png" alt="airrflow map" width="100">
    </p>

    1. QC and sequence assembly

    - Bulk
    - Raw read quality control, adapter trimming and clipping (`Fastp`).
    - Filter sequences by base quality (`pRESTO FilterSeq`).
    - Mask amplicon primers (`pRESTO MaskPrimers`).
    - Pair read mates (`pRESTO PairSeq`).
    - For UMI-based sequencing:
        - Cluster sequences according to similarity (optional for insufficient UMI diversity) (`pRESTO ClusterSets`).
        - Build consensus of sequences with the same UMI barcode (`pRESTO BuildConsensus`).
    - Assemble R1 and R2 read mates (`pRESTO AssemblePairs`).
    - Remove and annotate read duplicates (`pRESTO CollapseSeq`).
    - Filter out sequences that do not have at least 2 duplicates (`pRESTO SplitSeq`).
    - single cell
    - cellranger vdj
        - Assemble contigs
        - Annotate contigs
        - Call cells
        - Generate clonotypes

    2. V(D)J annotation and filtering (bulk and single-cell)

    - Assign gene segments with `IgBlast` using a germline reference (`Change-O AssignGenes`).
    - Annotate alignments in AIRR format (`Change-O MakeDB`)
    - Filter by alignment quality (locus matching v_call chain, min 200 informative positions, max 10% N nucleotides)
    - Filter productive sequences (`Change-O ParseDB split`)
    - Filter junction length multiple of 3
    - Annotate metadata (`EnchantR`)

    3. QC filtering (bulk and single-cell)

    - Bulk sequencing filtering:
    - Remove chimeric sequences (optional) (`SHazaM`, `EnchantR`)
    - Detect cross-contamination (optional) (`EnchantR`)
    - Collapse duplicates (`Alakazam`, `EnchantR`)
    - Single-cell QC filtering (`EnchantR`)
    - Remove cells without heavy chains.
    - Remove cells with multiple heavy chains.
    - Remove sequences in different samples that share the same `cell_id` and nucleotide sequence.
    - Modify `cell_id`s to ensure they are unique in the project.

    4. Clonal analysis (bulk and single-cell)

    - Find threshold for clone definition (`SHazaM`, `EnchantR`).
    - Create germlines and define clones, repertoire analysis (`SCOPer`, `EnchantR`).
    - Build lineage trees (`Dowser`, `IgphyML`, `RAxML`, `EnchantR`).

    5. Repertoire analysis and reporting

    - Custom repertoire analysis pipeline report (`Alakazam`).
    - Aggregate QC reports (`MultiQC`).

    ## Usage

    To run nf-core/airrflow with your data, prepare a tab-separated samplesheet with your input data. Depending on the input data type (bulk or single-cell, raw reads or assembled reads) the input samplesheet will vary. Please follow the [documentation on samplesheets](https://nf-co.re/airrflow/usage#input-samplesheet) for more details. An example samplesheet for running the pipeline on bulk BCR / TCR sequencing data in fastq format looks as follows:

    | sample_id | filename_R1                     | filename_R2                     | filename_I1                     | subject_id | species | pcr_target_locus | tissue | sex    | age | biomaterial_provider | single_cell | intervention   | collection_time_point_relative | cell_subset  |
    | --------- | ------------------------------- | ------------------------------- | ------------------------------- | ---------- | ------- | ---------------- | ------ | ------ | --- | -------------------- | ----------- | -------------- | ------------------------------ | ------------ |
    | sample01  | sample1_S8_L001_R1_001.fastq.gz | sample1_S8_L001_R2_001.fastq.gz | sample1_S8_L001_I1_001.fastq.gz | Subject02  | human   | IG               | blood  | NA     | 53  | sequencing_facility  | FALSE       | Drug_treatment | Baseline                       | plasmablasts |
    | sample02  | sample2_S8_L001_R1_001.fastq.gz | sample2_S8_L001_R2_001.fastq.gz | sample2_S8_L001_I1_001.fastq.gz | Subject02  | human   | TR               | blood  | female | 78  | sequencing_facility  | FALSE       | Drug_treatment | Baseline                       | plasmablasts |

    Each row represents a sample with fastq files (paired-end).

    For common **bulk sequencing protocols** we provide pre-set profiles that specify primers, UMI length, etc for common commercially available sequencing protocols. Please check the [Supported protocol profiles](#supported-protocol-profiles) for a full list of available profiles. An example command running the NEBNext UMI protocol profile with docker containers is:

    See the [usage documentation](https://nf-co.re/airrflow/usage) and the [parameter documentation](https://nf-co.re/airrflow/parameters) for more details on how to use the pipeline and all the available parameters.

    For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/airrflow/usage) and the [parameter documentation](https://nf-co.re/airrflow/parameters).

    ## Pipeline output

    To see the the results of a test run with a full size dataset refer to the [results](https://nf-co.re/airrflow/results) tab on the nf-core website pipeline page.
    For more details about the output files and reports, please refer to the
    [output documentation](https://nf-co.re/airrflow/output).

    ## Credits

    nf-core/airrflow was originally written by:

    - [Gisela Gabernet](https://github.com/ggabernet)
    - [Susanna Marquez](https://github.com/ssnn-airr)
    - [Alexander Peltzer](@apeltzer)
    - [Simon Heumos](@subwaystation)

    We thank the following people for their extensive assistance in the development of the pipeline:

    - [David Ladd](https://github.com/dladd)
    - [Friederike Hanssen](https://github.com/ggabernet/friederikehanssen)

    ## Contributions and Support

    If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

    For further information or help, don't hesitate to get in touch on the [Slack `#airrflow` channel](https://nfcore.slack.com/channels/airrflow) (you can join with [this invite](https://nf-co.re/join/slack)).

    ## Citations

    If you use nf-core/airrflow for your analysis, please cite the preprint as follows:

    > **nf-core/airrflow: an adaptive immune receptor repertoire analysis workflow employing the Immcantation framework**
    >
    > Gisela Gabernet, Susanna Marquez, Robert Bjornson, Alexander Peltzer, Hailong Meng, Edel Aron, Noah Y. Lee, Cole Jensen, David Ladd, Friederike Hanssen, Simon Heumos, nf-core community, Gur Yaari, Markus C. Kowarik, Sven Nahnsen, Steven H. Kleinstein.
    >
    > BioRxiv. 2024. doi: [10.1101/2024.01.18.576147](https://doi.org/10.1101/2024.01.18.576147).

    The specific pipeline version using the following DOI: [10.5281/zenodo.2642009](https://doi.org/10.5281/zenodo.2642009)

    Please also cite all the tools that are being used by the pipeline. An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

    You can cite the `nf-core` publication as follows:

    > **The nf-core framework for community-curated bioinformatics pipelines.**
    >
    > Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
    >
    > _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

    """

    pvc_name: str = initialize(run_name=run_name)

    nextflow_runtime(
        run_name=run_name,
        pvc_name=pvc_name,
        input=input,
        mode=mode,
        outdir=outdir,
        email=email,
        miairr=miairr,
        library_generation_method=library_generation_method,
        race_linker=race_linker,
        vprimers=vprimers,
        cprimers=cprimers,
        vprimer_start=vprimer_start,
        cprimer_start=cprimer_start,
        cprimer_position=cprimer_position,
        primer_revpr=primer_revpr,
        umi_position=umi_position,
        umi_length=umi_length,
        umi_start=umi_start,
        index_file=index_file,
        trim_fastq=trim_fastq,
        adapter_fasta=adapter_fasta,
        clip_r1=clip_r1,
        clip_r2=clip_r2,
        three_prime_clip_r1=three_prime_clip_r1,
        three_prime_clip_r2=three_prime_clip_r2,
        trim_nextseq=trim_nextseq,
        save_trimmed=save_trimmed,
        filterseq_q=filterseq_q,
        primer_consensus=primer_consensus,
        primer_mask_mode=primer_mask_mode,
        buildconsensus_maxerror=buildconsensus_maxerror,
        buildconsensus_maxgap=buildconsensus_maxgap,
        cluster_sets=cluster_sets,
        primer_r1_maxerror=primer_r1_maxerror,
        primer_r2_maxerror=primer_r2_maxerror,
        maskprimers_align=maskprimers_align,
        primer_extract_len=primer_extract_len,
        primer_maxlen=primer_maxlen,
        assemblepairs_sequential=assemblepairs_sequential,
        align_cregion=align_cregion,
        internal_cregion_sequences=internal_cregion_sequences,
        cregion_maxlen=cregion_maxlen,
        cregion_maxerror=cregion_maxerror,
        cregion_mask_mode=cregion_mask_mode,
        reassign=reassign,
        productive_only=productive_only,
        save_databases=save_databases,
        reference_fasta=reference_fasta,
        reference_igblast=reference_igblast,
        fetch_imgt=fetch_imgt,
        isotype_column=isotype_column,
        collapseby=collapseby,
        detect_contamination=detect_contamination,
        remove_chimeric=remove_chimeric,
        clonal_threshold=clonal_threshold,
        lineage_trees=lineage_trees,
        cloneby=cloneby,
        crossby=crossby,
        lineage_tree_builder=lineage_tree_builder,
        lineage_tree_exec=lineage_tree_exec,
        singlecell=singlecell,
        skip_all_clones_report=skip_all_clones_report,
        skip_report_threshold=skip_report_threshold,
        reference_10x=reference_10x,
        report_rmd=report_rmd,
        report_css=report_css,
        report_logo=report_logo,
        report_logo_img=report_logo_img,
        skip_report=skip_report,
        skip_multiqc=skip_multiqc,
        multiqc_title=multiqc_title,
        multiqc_methods_description=multiqc_methods_description,
    )


LaunchPlan(
    nf_nf_core_airrflow,
    "Small Test",
    {
        "input": [
            SampleSheet(
                sampleid="sampleID_1a",
                forwardreads=LatchFile(
                    "s3://latch-public/nf-core/ampliseq/test_data/1a_S103_L001_R1_001.fastq.gz"
                ),
                reversereads=LatchFile(
                    "s3://latch-public/nf-core/ampliseq/test_data/1a_S103_L001_R2_001.fastq.gz"
                ),
                run=None,
            ),
            SampleSheet(
                sampleid="sampleID_1",
                forwardreads=LatchFile(
                    "s3://latch-public/nf-core/ampliseq/test_data/1_S103_L001_R1_001.fastq.gz"
                ),
                reversereads=LatchFile(
                    "s3://latch-public/nf-core/ampliseq/test_data/1_S103_L001_R2_001.fastq.gz"
                ),
                run=None,
            ),
            SampleSheet(
                sampleid="sampleID_2a",
                forwardreads=LatchFile(
                    "s3://latch-public/nf-core/ampliseq/test_data/2a_S115_L001_R1_001.fastq.gz"
                ),
                reversereads=LatchFile(
                    "s3://latch-public/nf-core/ampliseq/test_data/2a_S115_L001_R2_001.fastq.gz"
                ),
                run=None,
            ),
            SampleSheet(
                sampleid="sampleID_2",
                forwardreads=LatchFile(
                    "s3://latch-public/nf-core/ampliseq/test_data/2_S115_L001_R1_001.fastq.gz"
                ),
                reversereads=LatchFile(
                    "s3://latch-public/nf-core/ampliseq/test_data/2_S115_L001_R2_001.fastq.gz"
                ),
                run=None,
            ),
        ],
    },
)
