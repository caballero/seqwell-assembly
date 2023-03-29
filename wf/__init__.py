"""
SeqWell assembly workflow
"""

import os
import subprocess
from pathlib import Path

from latch import small_task, workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import LatchAuthor, LatchFile, LatchMetadata, LatchParameter, LatchDir

@small_task
def main_task(
    readsdir: LatchDir = LatchDir(
        "s3://seqWell_plasmid_fastqs/"
    ),
    plates: str = 'SeqWell-07SEP21_SO10828_FASTQ',
    run: str = '20210917_Example',
    downsample: int = 10000,
    analysis: str = "assembly"
) -> LatchDir:

    # directory of output
    local_dir = "/root/work/" 
    local_prefix = os.path.join(local_dir, run, analysis)

    main_cmd = [
        "/root/bin/micromamba",
        "run",
        "-n",
        "seqWell-nf",
        "/bin/bash",
        "-c",
        f"""
        nextflow run \
            /root/wf-nf/fastq_assembly.nf \
            --plates {plates} \
            --downsample {downsample} \
            --run {run} \
            --analysis {analysis} \
            --readsdir {readsdir.local_path} \
            -work-dir {local_dir}
        """,
    ]

    subprocess.run(main_cmd, check=True)
    return LatchDir(local_prefix, f"latch:///seqwell/{run}/{analysis}")


"""The metadata included here will be injected into your interface."""
metadata = LatchMetadata(
    display_name="SeqWell assembly workflow",
    documentation="",
    author=LatchAuthor(
        name="Juan Caballero",
        email="juan.caballero.perez@gmail.com",
        github="github.com/caballero",
    ),
    repository="https://github.com/your-repo",
    license="",
    parameters={
        "run": LatchParameter(
            display_name="Run name",
            description="Run name",
            placeholder="run_001",
            batch_table_column=True,  # Show this parameter in batched mode.
        ),
        "plates": LatchParameter(
            display_name="Plate name",
            description="Plate name",
            placeholder="plate_001",
            batch_table_column=True,  # Show this parameter in batched mode.
        ),
        "downsample": LatchParameter(
            display_name="Down sample Fastq",
            description="Limit the total reads to be processed, 0 will process all reads",
            batch_table_column=True,  # Show this parameter in batched mode.
        ),
        "readsdir": LatchParameter(
            display_name="Select the Fastq directory location",
            description="Specify the location of the fastq files.",
            batch_table_column=True,  # Show this parameter in batched mode.
        ),
    },
    tags=[],
)


@workflow(metadata)
def run_wf(
    readsdir: LatchDir,
    plates: str,
    run: str,
    downsample: int,
) -> LatchDir:
    """SeqWell assembly workflow

    Performs assembly of raw short-reads and identify circularized sequences (plasmids)
    ----

    ## Inputs

    - Paired-ends Fastq reads

    ## Outputs

    - Summary CSV with assembly statistics
    - Sequences in Fasta format
    - Assembly graphs (GFA) 
    - Coverage plots (PNG)
    - Intermediate files for read aligns

    ## Workflow

    1. Down sample reads if required
    2. Preprocess raw reads with BBmap (adapter trimming, merge and clean)
    3. Aseembly with Uniclycler
    4. Circularize sequences
    5. Align reads to assembled sequences
    6. Generate statistics

    """
    return main_task(
                readsdir=readsdir, 
                plates=plates, 
                run=run, 
                downsample=downsample, 
                analysis='assembly')


"""
Add test data with a LaunchPlan. Provide default values in a dictionary with
the parameter names as the keys. These default values will be available under
the 'Test Data' dropdown at console.latch.bio.
"""
LaunchPlan(
    run_wf,
    "Test Data",
    {
        "readsdir": LatchDir(
            "s3://seqWell_plasmid_fastqs/"
        ),
        "plates": 'SeqWell-07SEP21_SO10828_FASTQ',
        "run": '20210917_Example',
        "downsample": 10000,
    },
)
