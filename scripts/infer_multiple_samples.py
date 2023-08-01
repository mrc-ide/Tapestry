import os
import click
import subprocess
import pandas as pd
from dataclasses import dataclass


# ================================================================================
# Utilities
#
# ================================================================================


@dataclass(frozen=False)
class SimulatedSequenceData:
    """
    On top of the simulated VCF, simulated sequence data
    includes a summary CSV file and a CSV holding all IBD
    segmenents.

    These are assumed to have the same `prefix` as the VCF,
    with fixed suffixes.

    """

    vcf: str
    source_dir: str = ""
    vcf_prefix: str = ""
    summary_csv: str = ""
    segment_csv: str = ""
    results_dir: str = ""

    def __post_init__(self):
        self.source_dir = os.path.dirname(self.vcf)
        self.vcf_prefix = os.path.basename(self.vcf).replace(".gz", "").replace(".vcf", "")
        self.summary_csv = f"{self.source_dir}/{self.vcf_prefix}.summary.csv"
        self.segment_csv = f"{self.source_dir}/{self.vcf_prefix}.ibd_segments.csv"
        self.results_dir = produce_dir(
            self.source_dir.replace("example_data", "results")
        )


def produce_dir(*args):
    """
    Produce a new directory by concatenating `args`,
    if it does not already exist
    params
        *args: str1, str2, str3 ...
            Comma-separated strings which will
            be combined to produce the directory,
            e.g. str1/str2/str3
    returns
        dir_name: str
            Directory name created from *args.
    """

    # Define directory path
    dir_name = os.path.join(*args)

    # Create if doesn't exist
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    return dir_name


# ================================================================================
# Running Tapestry from python
#
# ================================================================================


def run_tapestry_infer(input_vcf: str, 
                       sample_name: str, 
                       build="release", 
                       dry_run=False, 
                       **kwargs):
    """
    Run Tapestry
    
    """
    
    # Build cmd
    cmd = f"{build}/tapestry infer"
    cmd += f" -i {input_vcf}"
    cmd += f" -s {sample_name}"
    cmd += f" {' '.join([f'{k} {v}' for k,v in kwargs.items()])}"
    
    if dry_run:
        print(cmd)
        return
    
    # Run
    subprocess.run(cmd, shell=True, check=True)


# ================================================================================
# Entry point
#
# ================================================================================


@click.command(short_help="Run Tapestry for multiple sample from a VCF")
@click.option(
    "-i",
    "--input_vcf",
    type=click.Path(exists=True),
    required=True,
    help="Path to input VCF.",
)
@click.option(
    "-n",
    "--n_samples",
    type=int,
    default=None,
    required=False,
    help="Number of samples.",
)
@click.option(
    "-w",
    "--w_proposal",
    type=float,
    default=0.05,
    required=False,
    help="Proportion-titre proposal SD ~N(0, w)."
)
@click.option(
    "-v",
    "--var_wsaf",
    type=float,
    default=100,
    required=False,
    help="Variance in WSAF."
)
def main(input_vcf, n_samples, w_proposal, var_wsaf):
    """
    Run `Tapestry` for multiple samples

    """

    # Collect sequencing data   
    data = SimulatedSequenceData(input_vcf)

    # Get samples to run
    sample_df = pd.read_csv(data.summary_csv)
    print("Loaded {summary_df.shape[0]} samples from {data.summary_csv}.")
    if n_samples is not None:
        print(f"  Downsamling to {n_samples} randomly.")
        sample_df.sample(n_samples, inplace=True)
    sample_names = sample_df["sample_id"].tolist()

    # Prepare arguments accessible through this API
    args_dt = {
        "--w_proposal": w_proposal,
        "--var_wsaf": var_wsaf
    }

    # Run over multiple samples
    for sample_name in sample_names:
        args_dt["--output_dir"] = f"{data.results_dir}/{sample_name}"
        run_tapestry_infer(
            input_vcf=input_vcf,
            sample_name=sample_name,
            **args_dt
        )


if __name__ == "__main__":
    main()

    