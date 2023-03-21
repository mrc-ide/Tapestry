import os
import click
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec


# ================================================================================
# Parameters
#
# ================================================================================


OUTPUT_DIR = "results"
SUMMARY_CSV = "example_data/simulated_infections.summary.csv"


# ================================================================================
# Fuctions
#
# ================================================================================


class ProportionsTracePlotter:
    
    # PARAMETERS
    COL_PALETTE = "Set1"
    BURN_COL = "darkgrey"
    
    def __init__(self, mcmc_path, prop_path):
        """
        Class to plot proportion trace from an MCMC
        
        """
        
        # Dataframes
        self.mcmc_df = pd.read_csv(mcmc_path)
        self.prop_df = pd.read_csv(prop_path)
        self.results_df = self._merge_results()
        self.burn_df = self.results_df.query("phase == 'burn'")
        self.sample_df = self.results_df.query("phase == 'sample'")
        self.n_burn = self.burn_df.shape[0]
        self.n_sample = self.sample_df.shape[0]
        assert self.n_burn > 0, "No burn phase found."
        assert self.n_sample > 0, "No sample phase found."
        
        # Derived statistics
        self.K = None
        self.prop_columns = None
        self.col_dt = None
        
        # Derived statistics set here
        self._set_proportions()
        
    @classmethod
    def from_output_dir(cls, output_dir):
        mcmc_path = f"{output_dir}/mcmc.diagnostics.csv"
        prop_path = f"{output_dir}/mcmc.parameters.csv"
        return cls(mcmc_path, prop_path)
    
    def _merge_results(self):
        """
        Merge the MCMC diagnostic and proportion data
        
        """
        
        return pd.merge(self.mcmc_df, self.prop_df, on="iter")
    
    def _set_proportions(self):
        """
        Set proportion-related information
        
        """
        
        # Proportions
        self.prop_columns = [c for c in self.results_df.columns
                             if c.startswith("prop")]
        assert len(self.prop_columns) > 1, "Must be a column starting `prop`."
        
        # COI
        self.K = len(self.prop_columns)
        
        # Colors
        self.col_dt = dict(zip(
            self.prop_columns,
            sns.color_palette(self.COL_PALETTE, self.K)
        ))
        
    def _plot_trace(self, ax, title=None):
        """
        Build the proportion trace plot
        
        """
        
        for prop in self.prop_columns:
            # BURN
            # Initial value
            ax.scatter(
                x=self.burn_df["iter"].iloc[0],
                y=self.burn_df[prop].iloc[0],
                color=self.BURN_COL,
                marker='o',
                s=50,
                clip_on=False,
                zorder=10,
            )
            # Trace
            ax.plot(
                self.burn_df["iter"],
                self.burn_df[prop],
                color=self.BURN_COL,
                lw=1.5,
                clip_on=False,
                zorder=10
            )
            
            # SAMPLE
            ax.plot(
                self.sample_df["iter"],
                self.sample_df[prop],
                color=self.col_dt[prop],
                lw=1.5,
                clip_on=False,
                zorder=10
            )
        
        # Axis and Limits
        ax.set_axisbelow(True)
        ax.set_ylim((0, 1))
        ax.set_xlim((0, self.n_burn + self.n_sample))

        # Ticks
        ax.xaxis.set_major_locator(plt.MultipleLocator(200))
        ax.xaxis.set_minor_locator(plt.MultipleLocator(100))
        ax.yaxis.set_major_locator(plt.MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.025))

        # Grid
        ax.grid(which='minor', axis='y', ls='dotted', alpha=0.5)
        ax.grid(which='major', axis='y', ls='solid', alpha=0.5)

        # Burn
        ax.add_patch(Rectangle(xy=(0, 0),
                               height=1.0, 
                               width=self.n_burn,
                               color=self.BURN_COL,
                               alpha=0.5,
                               zorder=-10
                              )
                    )
        ax.annotate(xy=(0, 0.95), 
                    xycoords="data",
                    text=" $\it{Burn}$", 
                    color="white")
        ax.annotate(xy=(self.n_burn, 0.95), 
                    xycoords="data", 
                    text=" $\it{Sample}$", 
                    color="black")

        # Labels
        ax.set_xlabel("Iterations")
        ax.set_ylabel("Proportions")
        if title is not None:
            ax.set_title(title, loc="left")
            
    
    def _plot_histogram(self, ax, show_mean=True):
        """
        Build the mariginal histogram
        
        """
        
        data = [self.burn_df[c] for c in self.prop_columns]
        data += [self.sample_df[c] for c in self.prop_columns]

        ax.hist(data,
                bins=np.linspace(0, 1, 50),
                ec='black',
                stacked=True,
                orientation='horizontal',
                color=[self.BURN_COL]*self.K + list(self.col_dt.values())
               )


        # Posterior mean
        if show_mean:
            for c in self.prop_columns:
                mu = self.sample_df[c].mean()
                ax.axhline(y=mu, 
                           lw=2,
                           color=self.col_dt[c],
                           zorder=-3
                          )
                ax.annotate(xy=(ax.get_xlim()[1], mu),
                            va="center",
                            xycoords="data",
                            text=f" $\mu=${mu:.04f}",
                            color=self.col_dt[c]
                           )

        # Ticks
        ax.yaxis.set_major_locator(plt.MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.025))

        # Grid
        ax.grid(which='minor', axis='y', ls='dotted', alpha=0.5)
        ax.grid(which='major', axis='y', ls='solid', alpha=0.5)

        # Hide ticks
        ax.tick_params(axis='both', labelleft=False, labelbottom=False, left=False, bottom=False)
        
    
    def plot(self, title=None, output_path=None):
        """
        Plot trace and histogram in combined figure
        
        """
        
        # Prepare canvas
        fig = plt.figure(figsize=(10, 4))
        fig.subplots_adjust(wspace=0.075)

        gs = GridSpec(nrows=1, ncols=6)
        ax_trace = plt.subplot(gs[:5])
        ax_hist = plt.subplot(gs[5], sharey=ax_trace)
        
        # Plot
        self._plot_trace(ax=ax_trace, title=title)
        self._plot_histogram(ax=ax_hist)
        
        # Save
        if output_path is not None:
            fig.savefig(
                output_path,
                bbox_inches="tight",
                pad_inches=0.5
            )
            plt.close(fig)



class DiagnosticsTracePlotter:
    
    # PARAMETERS
    COL_PALETTE = "Set1"
    BURN_COL = "darkgrey"
    SAMPLE_AR_COL = "forestgreen"
    SAMPLE_LLK_COL = "steelblue"
    
    def __init__(self, mcmc_path):
        """
        Class to plot proportion trace from an MCMC
        
        """
        
        # Dataframes
        self.mcmc_df = pd.read_csv(mcmc_path)
        self.burn_df = self.mcmc_df.query("phase == 'burn'")
        self.sample_df = self.mcmc_df.query("phase == 'sample'")
        self.n_burn = self.burn_df.shape[0]
        self.n_sample = self.sample_df.shape[0]
        self.n_iter = self.n_burn + self.n_sample
        assert self.n_burn > 0, "No burn phase found."
        assert self.n_sample > 0, "No sample phase found."
        
    @classmethod
    def from_output_dir(cls, output_dir):
        mcmc_path = f"{output_dir}/mcmc.diagnostics.csv"
        return cls(mcmc_path)
    
    
    def _acceptance_plotter(self, ax, title=None):
        """
        Plot trace of running acceptance rate
        
        """
        # BURN
        ax.plot(
            self.burn_df["iter"],
            self.burn_df["acceptance_rate"],
            color=self.BURN_COL
        )

        # SAMPLE
        ax.plot(
            self.sample_df["iter"],
            self.sample_df["acceptance_rate"],
            color=self.SAMPLE_AR_COL
        )

        # Sample posterior
        show_mean = True
        if show_mean:
            mu = self.sample_df.iloc[-1]["acceptance_rate"]
            ax.annotate(xy=(self.n_iter, mu),
                        va="center",
                        xycoords="data",
                        text=f" $\mu=${mu:.04f}",
                        color=self.SAMPLE_AR_COL
                       )

        # Axis and Limits
        ax.set_axisbelow(True)
        ax.set_ylim((0, 1))
        ax.set_xlim((0, self.n_iter))

        # Ticks
        ax.xaxis.set_major_locator(plt.MultipleLocator(200))
        ax.xaxis.set_minor_locator(plt.MultipleLocator(100))
        ax.yaxis.set_major_locator(plt.MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.025))

        # Grid
        ax.grid(which='minor', axis='y', ls='dotted', alpha=0.5)
        ax.grid(which='major', axis='y', ls='solid', alpha=0.5)

        # Burn
        ax.add_patch(Rectangle(xy=(0, 0),
                               height=1.0, 
                               width=self.n_burn,
                               color=self.BURN_COL,
                               alpha=0.5,
                               zorder=-10
                              )
                    )
        ax.annotate(xy=(0, 0.95), 
                    xycoords="data",
                    text=" $\it{Burn}$", 
                    color="white")
        ax.annotate(xy=(self.n_burn, 0.95), 
                    xycoords="data", 
                    text=" $\it{Sample}$", 
                    color="black")

        # Labels
        ax.set_ylabel("Acceptance Rate")
        
        if title is not None:
            ax.set_title(title, loc="left")
        
        # Hide ticks
        ax.tick_params(axis='both', labelbottom=False)
        
    def _logposterior_plotter(self, ax):
        """
        Plot trace of log-posterior
        
        """

        # BURN
        ax.plot(
            self.burn_df["iter"],
            self.burn_df["logposterior"],
            color=self.BURN_COL
        )

        # SAMPLE
        ax.plot(
            self.sample_df["iter"],
            self.sample_df["logposterior"],
            color=self.SAMPLE_LLK_COL
        )

        # Sample posterior
        show_mean = True
        if show_mean:
            mu = self.sample_df.iloc[-1]["logposterior"]
            ax.annotate(xy=(self.n_iter, mu),
                        va="center",
                        xycoords="data",
                        text=f" $\mu=${mu:.04f}",
                        color=self.SAMPLE_LLK_COL
                       )

        # Axis and Limits
        ax.set_axisbelow(True)
        ax.set_xlim((0, self.n_iter))

        # Ticks
        ax.xaxis.set_major_locator(plt.MultipleLocator(200))
        ax.xaxis.set_minor_locator(plt.MultipleLocator(100))

        # Grid
        ax.grid(which='minor', axis='y', ls='dotted', alpha=0.5)
        ax.grid(which='major', axis='y', ls='solid', alpha=0.5)

        # Burn
        ax.add_patch(Rectangle(xy=(ax.get_xlim()[0], ax.get_ylim()[0]),
                               height=ax.get_ylim()[1] - ax.get_ylim()[0], 
                               width=self.n_burn,
                               color=self.BURN_COL,
                               alpha=0.5,
                               zorder=-10
                              )
                    )

        # Labels
        ax.set_xlabel("Iterations")
        ax.set_ylabel("Log-posterior")
        
        
    def plot(self, title=None, output_path=None):
        """
        Plot trace and histogram in combined figure
        
        """
        
        # Prepare canvas
        fig = plt.figure(figsize=(10, 7))
        fig.subplots_adjust(hspace=0.05)

        gs = GridSpec(nrows=2, ncols=1)
        ax_accept = plt.subplot(gs[0])
        ax_llk = plt.subplot(gs[1], sharex=ax_accept)
        
        # Plot
        self._acceptance_plotter(ax=ax_accept, title=title)
        self._logposterior_plotter(ax=ax_llk)
        
        # Save
        if output_path is not None:
            fig.savefig(
                output_path,
                bbox_inches="tight",
                pad_inches=0.5
            )
            plt.close(fig)


# ================================================================================
# Entry point
#
# ================================================================================

@click.command(short_help="Make plots of Tapestry outputs.")
@click.option(
    "-c",
    "--summary_csv",
    type=click.Path(exists=True),
    default=SUMMARY_CSV,
    required=False,
    help="Path to summary CSV.",
)
def main(summary_csv):

    # Prepare
    output_dir = OUTPUT_DIR
    summary_df = pd.read_csv(summary_csv)

    # Iterate over samples
    sample_names = [s for s in os.listdir(output_dir) if s.startswith("SMI")]
    for sample_name in sample_names:
        
        # Create sample directory
        print(f"Plotting sample: {sample_name}...")
        sample_dir = f"{output_dir}/{sample_name}"
        
        # Extract parameters
        sample_info = summary_df.query("sample_id == @sample_name").squeeze()
        
        # Prepare title
        title = f"{sample_name} | "
        title += ", ".join([f"{sample_info[p]:.03f}" for p in sample_info.keys() if p.startswith("prop")])
        title += " | $f_{IBD}=$" + f"{sample_info['f_ibd']:.03}"
        
        # Plot proportions
        prop_plotter = ProportionsTracePlotter.from_output_dir(sample_dir)
        prop_plotter.plot(
            title=title, 
            output_path=f"{sample_dir}/plot.parameters.{sample_name}.pdf"
        )

        # Plot diagnostics
        diag_plotter = DiagnosticsTracePlotter.from_output_dir(sample_dir)
        diag_plotter.plot(
            title=sample_name,
            output_path=f"{sample_dir}/plot.diagnostics.{sample_name}.pdf"
        )

if __name__ == "__main__":
    main()

