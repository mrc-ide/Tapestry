
import os
import click
import subprocess
import pandas as pd
from typing import List



# ================================================================================
# Parameters
#
# ================================================================================


OUTPUT_DIR = "results"
INPUT_VCF = "example_data/simulated_infections.DRCongo.K02.vcf"
SUMMARY_CSV = "example_data/simulated_infections.DRCongo.K02.summary.csv"
N_SAMPLES = 5


# ================================================================================
# Fuctions
#
# ================================================================================


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


def run_tapestry_infer_multiple_samples(input_vcf: str, 
                                        sample_names: List[str], 
                                        build="release", 
                                        **kwargs):
    """
    Run Tapestry across multiple samples

    """

    for sample_name in sample_names:
        run_tapestry_infer(
            input_vcf=input_vcf,
            sample_name=sample_name,
            build=build,
            **kwargs
        )


# ================================================================================
# Entry point
#
# ================================================================================


@click.command(short_help="Run Tapestry for multiple sample from a VCF")
@click.option(
    "-i",
    "--input_vcf",
    type=click.Path(exists=True),
    default=INPUT_VCF,
    required=False,
    help="Path to input VCF.",
)
@click.option(
    "-c",
    "--summary_csv",
    type=click.Path(exists=True),
    default=SUMMARY_CSV,
    required=False,
    help="Path to summary CSV.",
)
@click.option(
    "-n",
    "--n_samples",
    type=int,
    default=5,
    required=False,
    help="Number of samples.",
)
@click.option(
    "-K",
    "--coi",
    type=int,
    default=2,
    required=False,
    help="COI to run."
)
@click.option(
    "-w",
    "--w_proposal",
    type=float,
    default=0.5,
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
def main(input_vcf, summary_csv, n_samples, coi, w_proposal, var_wsaf):
    """
    Run `Tapestry` over a selection of samples

    """

    # Get sample names
    summary_df = (
        pd.read_csv(summary_csv)
        .head(n_samples)
    )
    sample_names = summary_df["sample_id"].tolist()

    # Produce output directory
    results_dir = produce_dir(OUTPUT_DIR)

    # Preparer arguments
    args_dt = {
        "--output_dir": "",
        "--COI": coi,
        "--w_proposal": w_proposal,
        "--var_wsaf": var_wsaf
    }

    # Run over multiple samples
    for sample_name in sample_names:
        args_dt["--output_dir"] = f"{results_dir}/{sample_name}"
        run_tapestry_infer(
            input_vcf=input_vcf,
            sample_name=sample_name,
            **args_dt
        )


if __name__ == "__main__":
    main()

