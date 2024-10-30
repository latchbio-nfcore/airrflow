import typing

from latch.types.directory import LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import (
    LatchRule,
    NextflowParameter,
    Params,
    Section,
    Spoiler,
    Text,
)

flow = [
    Section(
        "Input/Output Options",
        Text("Define where the pipeline should find input data and save output data."),
        Params(
            "input",
            "mode",
        ),
    ),
    Section(
        "Protocol",
        Text("Experimental protocol used to generate the data"),
        Params("library_generation_method", "race_linker"),
    ),
    Section(
        "Output Directory",
        Params("run_name"),
        Text("Parent directory for outputs"),
        Params("outdir"),
    ),
    Spoiler(
        "Primer Handling",
        Text(
            "Define the primer region start and how to deal with the primer alignment."
        ),
        Params(
            "vprimers",
            "cprimers",
            "vprimer_start",
            "cprimer_start",
            "cprimer_position",
            "primer_revpr",
        ),
    ),
    Spoiler(
        "UMI Barcode Handling",
        Text("Define how UMI barcodes should be treated."),
        Params("umi_position", "umi_length", "umi_start", "index_file"),
    ),
    Spoiler(
        "Adapter Trimming",
        Text("Options for adapter trimming and read clipping"),
        Params(
            "trim_fastq",
            "adapter_fasta",
            "clip_r1",
            "clip_r2",
            "three_prime_clip_r1",
            "three_prime_clip_r2",
            "trim_nextseq",
            "save_trimmed",
        ),
    ),
    Spoiler(
        "Sequence Assembly Options",
        Text("Options for the pRESTO sequence assembly processes"),
        Params(
            "filterseq_q",
            "primer_consensus",
            "primer_mask_mode",
            "buildconsensus_maxerror",
            "buildconsensus_maxgap",
            "cluster_sets",
            "primer_r1_maxerror",
            "primer_r2_maxerror",
            "maskprimers_align",
            "primer_extract_len",
            "primer_maxlen",
            "assemblepairs_sequential",
            "align_cregion",
            "internal_cregion_sequences",
            "cregion_maxlen",
            "cregion_maxerror",
            "cregion_mask_mode",
        ),
    ),
    Spoiler(
        "VDJ Annotation Options",
        Text("Options for the VDJ annotation processes."),
        Params(
            "reassign",
            "productive_only",
            "save_databases",
            "reference_fasta",
            "reference_igblast",
            "fetch_imgt",
            "isotype_column",
        ),
    ),
    Spoiler(
        "Bulk Filtering Options",
        Text("Options for bulk sequence filtering after VDJ assignment."),
        Params("collapseby", "detect_contamination", "remove_chimeric"),
    ),
    Spoiler(
        "Clonal Analysis Options",
        Text("Define how the B-cell clonal trees should be calculated."),
        Params(
            "clonal_threshold",
            "lineage_trees",
            "cloneby",
            "crossby",
            "lineage_tree_builder",
            "lineage_tree_exec",
            "singlecell",
            "skip_all_clones_report",
            "skip_report_threshold",
        ),
    ),
    Spoiler(
        "Single Cell Analysis Options",
        Text("Options specific for raw single cell input."),
        Params("reference_10x"),
    ),
    Spoiler(
        "Report Options",
        Params(
            "report_rmd",
            "report_css",
            "report_logo",
            "report_logo_img",
            "skip_report",
            "skip_multiqc",
        ),
    ),
    Spoiler(
        "Generic Options",
        Params("multiqc_title", "multiqc_methods_description", "email", "miairr"),
    ),
]

generated_parameters = {
    "run_name": NextflowParameter(
        type=str,
        display_name="Run Name",
        description="Name of run",
        batch_table_column=True,
        rules=[
            LatchRule(
                regex=r"^[a-zA-Z0-9_-]+$",
                message="Run name must contain only letters, digits, underscores, and dashes. No spaces are allowed.",
            )
        ],
    ),
    "input": NextflowParameter(
        type=LatchFile,
        default=None,
        section_title="Input/output options",
        description="Path to comma-separated file containing information about the samples in the experiment.",
    ),
    "mode": NextflowParameter(
        type=typing.Optional[str],
        default="fastq",
        section_title=None,
        description='Specify the processing mode for the pipeline. Available options are "fastq" and "assembled".',
    ),
    "outdir": NextflowParameter(
        type=LatchOutputDir,
        default=None,
        section_title=None,
        description="The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.",
    ),
    "email": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description="Email address for completion summary.",
    ),
    "miairr": NextflowParameter(
        type=typing.Optional[str],
        default="${projectDir}/assets/reveal/mapping_MiAIRR_BioSample_v1.3.1.tsv",
        section_title=None,
        description="Path to MiAIRR-BioSample mapping",
    ),
    "library_generation_method": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title="Protocol",
        description="Protocol used for the V(D)J amplicon sequencing library generation.",
    ),
    "race_linker": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description="Path to fasta file containing the linker sequence, if no V-region primers were used but a linker sequence is present (e.g. 5' RACE SMARTer TAKARA protocol).",
    ),
    "vprimers": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title="Primer handling",
        description="Path to a fasta file containinc the V-region primer sequences.",
    ),
    "cprimers": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description="Path to a fasta file containing the C-region primer sequences.",
    ),
    "vprimer_start": NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        description="Start position of V region primers (without counting the UMI barcode).",
    ),
    "cprimer_start": NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        description="Start position of C region primers (without counting the UMI barcode).",
    ),
    "cprimer_position": NextflowParameter(
        type=typing.Optional[str],
        default="R1",
        section_title=None,
        description="Indicate if C region primers are in the R1 or R2 reads.",
    ),
    "primer_revpr": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Specify to match the tail-end of the sequence against the reverse complement of the primers. This also reverses the behavior of the --start argument, such that start position is relative to the tail-end of the sequence. (default: False)Maximum scoring error for the Presto MaxPrimer process for the C and/or V region primers identification.",
    ),
    "umi_position": NextflowParameter(
        type=typing.Optional[str],
        default="R1",
        section_title="UMI barcode handling",
        description="Indicate if UMI indices are recorded in the R1 (default) or R1 fastq file.",
    ),
    "umi_length": NextflowParameter(
        type=typing.Optional[int],
        default=-1,
        section_title=None,
        description="UMI barcode length in nucleotides. Set to 0 if no UMIs present.",
    ),
    "umi_start": NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        description="UMI barcode start position in the index read.",
    ),
    "index_file": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Indicate if UMI indices are recorded in a separate index file.",
    ),
    "trim_fastq": NextflowParameter(
        type=typing.Optional[bool],
        default=True,
        section_title="Adapter trimming",
        description="Whether to trim adapters in fastq reads with fastp.",
    ),
    "adapter_fasta": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description="Fasta file with adapter sequences to be trimmed.",
    ),
    "clip_r1": NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        description="Number of bases to clip 5' in R1 reads.",
    ),
    "clip_r2": NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        description="Number of bases to clip 5' in R2 reads.",
    ),
    "three_prime_clip_r1": NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        description="Number of bases to clip 3' in R1 reads.",
    ),
    "three_prime_clip_r2": NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        description="Number of bases to clip 3' in R2 reads.",
    ),
    "trim_nextseq": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Trim adapters specific for Nextseq sequencing",
    ),
    "save_trimmed": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Option to save trimmed reads.",
    ),
    "filterseq_q": NextflowParameter(
        type=typing.Optional[int],
        default=20,
        section_title="sequence assembly options",
        description="Quality threshold for pRESTO FilterSeq sequence filtering.",
    ),
    "primer_consensus": NextflowParameter(
        type=typing.Optional[float],
        default=0.6,
        section_title=None,
        description="Maximum error for building the primer consensus in the pRESTO Buildconsensus step.",
    ),
    "primer_mask_mode": NextflowParameter(
        type=typing.Optional[str],
        default="cut",
        section_title=None,
        description="Masking mode for the pRESTO MaskPrimer step. Available: cut, mask, trim, tag.",
    ),
    "buildconsensus_maxerror": NextflowParameter(
        type=typing.Optional[float],
        default=0.1,
        section_title=None,
        description="Maximum error for building the sequence consensus in the pRESTO BuildConsensus step.",
    ),
    "buildconsensus_maxgap": NextflowParameter(
        type=typing.Optional[float],
        default=0.5,
        section_title=None,
        description="Maximum gap for building the sequence consensus in the pRESTO BuildConsensus step.",
    ),
    "cluster_sets": NextflowParameter(
        type=typing.Optional[bool],
        default=True,
        section_title=None,
        description="Cluster sequences by similarity regardless of any annotation with pRESTO ClusterSets and annotate the cluster ID additionally to the UMI barcode.",
    ),
    "primer_r1_maxerror": NextflowParameter(
        type=typing.Optional[float],
        default=0.2,
        section_title=None,
        description="Maximum allowed error for R1 primer alignment.",
    ),
    "primer_r2_maxerror": NextflowParameter(
        type=typing.Optional[float],
        default=0.2,
        section_title=None,
        description="Maximum allowed error for R2 primer alignment.",
    ),
    "maskprimers_align": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Align primers instead of scoring them. Used for protocols without primer fixed positions.",
    ),
    "primer_extract_len": NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        description="Length of the extracted primers with MaskPrimer extract.",
    ),
    "primer_maxlen": NextflowParameter(
        type=typing.Optional[int],
        default=50,
        section_title=None,
        description="Maximum allowed primer length when aligning the primers.",
    ),
    "assemblepairs_sequential": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Use AssemblePairs sequential instead of AssemblePairs align when assembling read pairs.",
    ),
    "align_cregion": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Align internal C-region for a more precise isotype characterization.",
    ),
    "internal_cregion_sequences": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description="Provide internal C-region sequences for a more precise C-region characterization. Then also set the `align_cregion` flag.",
    ),
    "cregion_maxlen": NextflowParameter(
        type=typing.Optional[int],
        default=100,
        section_title=None,
        description="Maximum allowed length when aligning the internal C-region.",
    ),
    "cregion_maxerror": NextflowParameter(
        type=typing.Optional[float],
        default=0.3,
        section_title=None,
        description="Maximum allowed error when aligning the internal C-region.",
    ),
    "cregion_mask_mode": NextflowParameter(
        type=typing.Optional[str],
        default="tag",
        section_title=None,
        description="Mask mode for C-region alignment.",
    ),
    "reassign": NextflowParameter(
        type=typing.Optional[bool],
        default=True,
        section_title="VDJ annotation options",
        description="Whether to reassign genes if the input file is an AIRR formatted tabulated file.",
    ),
    "productive_only": NextflowParameter(
        type=typing.Optional[bool],
        default=True,
        section_title=None,
        description="Subset to productive  sequences.",
    ),
    "save_databases": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Save databases so you can use the cache in future runs.",
    ),
    "reference_fasta": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description="Path to the germline reference fasta.",
    ),
    "reference_igblast": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description="Path to the cached igblast database.",
    ),
    "fetch_imgt": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Set this flag to fetch the IMGT reference data at runtime.",
    ),
    "isotype_column": NextflowParameter(
        type=typing.Optional[str],
        default="c_call",
        section_title=None,
        description="Set the column in the AIRR rearrangement file that isotype information should be gathered from.",
    ),
    "collapseby": NextflowParameter(
        type=typing.Optional[str],
        default="sample_id",
        section_title="Bulk filtering options",
        description="Name of the field used to collapse duplicated sequences.",
    ),
    "detect_contamination": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Whether to run the process to detect contamination.",
    ),
    "remove_chimeric": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Whether to apply the chimera removal filter.",
    ),
    "clonal_threshold": NextflowParameter(
        type=typing.Optional[str],
        default="auto",
        section_title="Clonal analysis options",
        description="Set the clustering threshold Hamming distance value. Default: 'auto'",
    ),
    "lineage_trees": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Perform clonal lineage tree analysis.",
    ),
    "cloneby": NextflowParameter(
        type=typing.Optional[str],
        default="subject_id",
        section_title=None,
        description="Name of the field used to group data files to identify clones.",
    ),
    "crossby": NextflowParameter(
        type=typing.Optional[str],
        default="subject_id",
        section_title=None,
        description="Name of the field used to identify external groups used to identify a clonal threshold.",
    ),
    "lineage_tree_builder": NextflowParameter(
        type=typing.Optional[str],
        default="raxml",
        section_title=None,
        description="Lineage tree software to use to build trees within Dowser. If you change the default, also set the `lineage_tree_exec` parameter.",
    ),
    "lineage_tree_exec": NextflowParameter(
        type=typing.Optional[str],
        default="/usr/local/bin/raxml-ng",
        section_title=None,
        description="Path to lineage tree building executable.",
    ),
    "singlecell": NextflowParameter(
        type=typing.Optional[str],
        default="single_cell",
        section_title=None,
        description="Name of the field used to determine if a sample is single cell sequencing or not.",
    ),
    "skip_all_clones_report": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Skip report of EnchantR DefineClones for all samples together.",
    ),
    "skip_report_threshold": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Skip report of EnchantR FindThreshold for all samples together.",
    ),
    "reference_10x": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title="Single cell analysis options",
        description="Path to the reference directory required by cellranger. Can either be directory or tar.gz.",
    ),
    "report_rmd": NextflowParameter(
        type=typing.Optional[str],
        default="${projectDir}/assets/repertoire_comparison.Rmd",
        section_title="Report options",
        description="Custom report Rmarkdown file.",
    ),
    "report_css": NextflowParameter(
        type=typing.Optional[str],
        default="${projectDir}/assets/nf-core_style.css",
        section_title=None,
        description="Custom report style file in css format.",
    ),
    "report_logo": NextflowParameter(
        type=typing.Optional[str],
        default="${projectDir}/assets/nf-core-airrflow_logo_light.png",
        section_title=None,
        description="Custom logo for the report.",
    ),
    "report_logo_img": NextflowParameter(
        type=typing.Optional[str],
        default="${projectDir}/assets/nf-core-airrflow_logo_reports.png",
        section_title=None,
        description="Custom logo for the EnchantR reports.",
    ),
    "skip_report": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Skip repertoire analysis and report generation.",
    ),
    "skip_multiqc": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Skip multiqc report.",
    ),
    "multiqc_title": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title="Generic options",
        description="MultiQC report title. Printed as page header, used for filename if not otherwise specified.",
    ),
    "multiqc_methods_description": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description="Custom MultiQC yaml file containing HTML including a methods description.",
    ),
}
