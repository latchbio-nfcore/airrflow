from latch.types.directory import LatchDir
from latch.types.metadata import LatchAuthor, NextflowMetadata, NextflowRuntimeResources

from .parameters import generated_parameters

NextflowMetadata(
    display_name="nf-core/airrflow",
    author=LatchAuthor(
        name="nf-core",
    ),
    repository="https://github.com/latchbio-nfcore/airrflow",
    parameters=generated_parameters,
    runtime_resources=NextflowRuntimeResources(
        cpus=4,
        memory=8,
        storage_gib=100,
    ),
    log_dir=LatchDir("latch:///your_log_dir"),
)
