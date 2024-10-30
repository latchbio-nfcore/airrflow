import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional

import requests
from latch.executions import rename_current_execution, report_nextflow_used_storage
from latch.ldata.path import LPath
from latch.resources.tasks import custom_task, nextflow_runtime_task
from latch.types.directory import LatchOutputDir
from latch.types.file import LatchFile
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.nextflow.workflow import get_flag
from latch_cli.services.register.utils import import_module_by_path
from latch_cli.utils import urljoins

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)
import latch_metadata


@dataclass
class SampleSheet:
    sample: str
    fastq_1: LatchFile
    fastq_2: Optional[LatchFile]
    antibody: Optional[str]
    control: Optional[str]


class LibraryGenerationMethod(Enum):
    SPECIFIC_PCR_UMI = "specific_pcr_umi"
    SPECIFIC_PCR = "specific_pcr"
    DT_5P_RACE = "dt_5p_race"
    DT_5P_RACE_UMI = "dt_5p_race_umi"
    SC_10X_GENOMICS = "sc_10x_genomics"


class ProcessingMode(Enum):
    FASTQ = "fastq"
    ASSEMBLED = "assembled"


class CPrimerPosition(Enum):
    R1 = "R1"
    R2 = "R2"


class UMIPosition(Enum):
    R1 = "R1"
    R2 = "R2"


class LineageTreeBuilder(Enum):
    RAXML = "raxml"
    IGPHYML = "igphyml"


@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize(run_name: str) -> str:
    rename_current_execution(str(run_name))

    token = os.environ.get("FLYTE_INTERNAL_EXECUTION_ID")
    if token is None:
        raise RuntimeError("failed to get execution token")

    headers = {"Authorization": f"Latch-Execution-Token {token}"}

    print("Provisioning shared storage volume... ", end="")
    resp = requests.post(
        "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage-ofs",
        headers=headers,
        json={
            "storage_expiration_hours": 0,
            "version": 2,
        },
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]


@nextflow_runtime_task(cpu=4, memory=8, storage_gib=100)
def nextflow_runtime(
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
    shared_dir = Path("/nf-workdir")

    exec_name = _get_execution_name()
    if exec_name is None:
        print("Failed to get execution name.")
        exec_name = "unknown"

    latch_log_dir = urljoins("latch:///your_log_dir/nf_nf_core_airrflow", exec_name)
    print(f"Log directory: {latch_log_dir}")

    ignore_list = [
        "latch",
        ".latch",
        ".git",
        "nextflow",
        ".nextflow",
        "work",
        "results",
        "miniconda",
        "anaconda3",
        "mambaforge",
    ]

    shutil.copytree(
        Path("/root"),
        shared_dir,
        ignore=lambda src, names: ignore_list,
        ignore_dangling_symlinks=True,
        dirs_exist_ok=True,
    )

    cmd = [
        "/root/nextflow",
        "run",
        str(shared_dir / "main.nf"),
        "-work-dir",
        str(shared_dir),
        "-profile",
        "docker",
        "-c",
        "latch.config",
        "-resume",
        *get_flag("input", input),
        *get_flag("mode", mode),
        *get_flag("outdir", LatchOutputDir(f"{outdir.remote_path}/{run_name}")),
        *get_flag("email", email),
        *get_flag("miairr", miairr),
        *get_flag("library_generation_method", library_generation_method),
        *get_flag("race_linker", race_linker),
        *get_flag("vprimers", vprimers),
        *get_flag("cprimers", cprimers),
        *get_flag("vprimer_start", vprimer_start),
        *get_flag("cprimer_start", cprimer_start),
        *get_flag("cprimer_position", cprimer_position),
        *get_flag("primer_revpr", primer_revpr),
        *get_flag("umi_position", umi_position),
        *get_flag("umi_length", umi_length),
        *get_flag("umi_start", umi_start),
        *get_flag("index_file", index_file),
        *get_flag("trim_fastq", trim_fastq),
        *get_flag("adapter_fasta", adapter_fasta),
        *get_flag("clip_r1", clip_r1),
        *get_flag("clip_r2", clip_r2),
        *get_flag("three_prime_clip_r1", three_prime_clip_r1),
        *get_flag("three_prime_clip_r2", three_prime_clip_r2),
        *get_flag("trim_nextseq", trim_nextseq),
        *get_flag("save_trimmed", save_trimmed),
        *get_flag("filterseq_q", filterseq_q),
        *get_flag("primer_consensus", primer_consensus),
        *get_flag("primer_mask_mode", primer_mask_mode),
        *get_flag("buildconsensus_maxerror", buildconsensus_maxerror),
        *get_flag("buildconsensus_maxgap", buildconsensus_maxgap),
        *get_flag("cluster_sets", cluster_sets),
        *get_flag("primer_r1_maxerror", primer_r1_maxerror),
        *get_flag("primer_r2_maxerror", primer_r2_maxerror),
        *get_flag("maskprimers_align", maskprimers_align),
        *get_flag("primer_extract_len", primer_extract_len),
        *get_flag("primer_maxlen", primer_maxlen),
        *get_flag("assemblepairs_sequential", assemblepairs_sequential),
        *get_flag("align_cregion", align_cregion),
        *get_flag("internal_cregion_sequences", internal_cregion_sequences),
        *get_flag("cregion_maxlen", cregion_maxlen),
        *get_flag("cregion_maxerror", cregion_maxerror),
        *get_flag("cregion_mask_mode", cregion_mask_mode),
        *get_flag("reassign", reassign),
        *get_flag("productive_only", productive_only),
        *get_flag("save_databases", save_databases),
        *get_flag("reference_fasta", reference_fasta),
        *get_flag("reference_igblast", reference_igblast),
        *get_flag("fetch_imgt", fetch_imgt),
        *get_flag("isotype_column", isotype_column),
        *get_flag("collapseby", collapseby),
        *get_flag("detect_contamination", detect_contamination),
        *get_flag("remove_chimeric", remove_chimeric),
        *get_flag("clonal_threshold", clonal_threshold),
        *get_flag("lineage_trees", lineage_trees),
        *get_flag("cloneby", cloneby),
        *get_flag("crossby", crossby),
        *get_flag("lineage_tree_builder", lineage_tree_builder),
        *get_flag("lineage_tree_exec", lineage_tree_exec),
        *get_flag("singlecell", singlecell),
        *get_flag("skip_all_clones_report", skip_all_clones_report),
        *get_flag("skip_report_threshold", skip_report_threshold),
        *get_flag("reference_10x", reference_10x),
        *get_flag("report_rmd", report_rmd),
        *get_flag("report_css", report_css),
        *get_flag("report_logo", report_logo),
        *get_flag("report_logo_img", report_logo_img),
        *get_flag("skip_report", skip_report),
        *get_flag("skip_multiqc", skip_multiqc),
        *get_flag("multiqc_title", multiqc_title),
        *get_flag("multiqc_methods_description", multiqc_methods_description),
    ]

    print("Launching Nextflow Runtime")
    print(" ".join(cmd))
    print(flush=True)

    failed = False
    try:
        env = {
            **os.environ,
            "NXF_ANSI_LOG": "false",
            "NXF_HOME": "/root/.nextflow",
            "NXF_OPTS": "-Xms1536M -Xmx6144M -XX:ActiveProcessorCount=4",
            "NXF_DISABLE_CHECK_LATEST": "true",
            "NXF_ENABLE_VIRTUAL_THREADS": "false",
            "NXF_ENABLE_FS_SYNC": "true",
        }

        if False:
            env["LATCH_LOG_DIR"] = latch_log_dir

        subprocess.run(
            cmd,
            env=env,
            check=True,
            cwd=str(shared_dir),
        )
    except subprocess.CalledProcessError:
        failed = True
    finally:
        print()

        nextflow_log = shared_dir / ".nextflow.log"
        if nextflow_log.exists():
            remote = LPath(urljoins(latch_log_dir, "nextflow.log"))
            print(f"Uploading .nextflow.log to {remote.path}")
            remote.upload_from(nextflow_log)

        print("Computing size of workdir... ", end="")
        try:
            result = subprocess.run(
                ["du", "-sb", str(shared_dir)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=5 * 60,
            )

            size = int(result.stdout.split()[0])
            report_nextflow_used_storage(size)
            print(f"Done. Workdir size: {size / 1024 / 1024 / 1024: .2f} GiB")
        except subprocess.TimeoutExpired:
            print(
                "Failed to compute storage size: Operation timed out after 5 minutes."
            )
        except subprocess.CalledProcessError as e:
            print(f"Failed to compute storage size: {e.stderr}")
        except Exception as e:
            print(f"Failed to compute storage size: {e}")

    if failed:
        sys.exit(1)
