import pandas as pd
from unravel.lib.dirs import TapestrySampleOutputDir
from unravel.lib.chromosomes import ChromosomeInformation
from unravel.lib.plot_mcmc import DiagnosticsTracePlotter, ProportionsTracePlotter
from unravel.lib.plot_ibd import WSAFPlotter, IBDSegmentPlotter, CombinedPlotter


BED_COLS = ["chrom", "start", "end", "name"]


def sample(sample_dir: str) -> None:
    """
    Make informative plots for a individual tapestry sample

    """

    # Define tapestry directories
    tapestry_dirs = TapestrySampleOutputDir(sample_dir)

    # Iterate over COIs and plot
    for coi, coi_dirs in tapestry_dirs.coi.items():
        png_temp = f"{coi_dirs.coi_dir}/plot.{'{name}'}.png"
        print(f"Plotting for COI={coi}...")

        # MCMC Plots ---------------------------------------------
        prop_plotter = ProportionsTracePlotter(
            mcmc_path=coi_dirs.mcmc_diagnostics_csv, 
            prop_path=coi_dirs.mcmc_parameters_csv
        )
        prop_plotter.plot(
            output_path=png_temp.format(name="mcmc.proportions")
        )
        diag_plotter = DiagnosticsTracePlotter(
            mcmc_path=coi_dirs.mcmc_diagnostics_csv
        )
        diag_plotter.plot(
            output_path=png_temp.format(name="mcmc.diagnostics")
        )

        # IBD Plots ---------------------------------------------
        # Note that these start from the dataframes, rather than
        # from paths. Might want to change that?
        # Load
        ibd_pairwise_df = pd.read_csv(coi_dirs.ibd_pairwise_csv)
        try:
            ibd_segments_df = pd.read_csv(
                coi_dirs.ibd_segs_bed,
                sep="\t",
                header=None
            )
            ibd_segments_df.columns = BED_COLS
        except pd.errors.EmptyDataError:  # no segments
            ibd_segments_df = pd.DataFrame(columns=BED_COLS)

        # Plot
        chrom_info = ChromosomeInformation(ibd_pairwise_df)
        wsaf_plotter = WSAFPlotter(wsaf_df=ibd_pairwise_df, chrom_info=chrom_info)
        ibd_plotter = IBDSegmentPlotter(segs_df=ibd_segments_df, chrom_info=chrom_info)
        combined_plotter = CombinedPlotter(
            chrom_info=chrom_info, 
            wsaf_plotter=wsaf_plotter, 
            ibd_plotter=ibd_plotter
        )
        combined_plotter.plot(
            output_path=png_temp.format(name="fit.ibd_segments")
        )

