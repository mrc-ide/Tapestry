import os
import sys
import click
import pandas as pd
pd.options.mode.chained_assignment = None 

from collections import namedtuple
from dataclasses import dataclass
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
import seaborn as sns

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.infer_multiple_samples import SimulatedSequenceData, produce_dir



# ================================================================================
# Functions
#
# ================================================================================


@dataclass
class BEDRecord:
    chrom: str
    start: int
    end: int
    name: str = ""


def get_ibd_segments(chroms: List[str],
                     pos: List[int],
                     ibd_states: List[bool],
                     name: Optional[str] = ""
                    ) -> pd.DataFrame:
    """
    Get all IBD segments implied from a list indicating the IBD
    state at each position in a genome
    
    """
        
    # Init
    t = 0
    chrom = chroms[t]
    ibd_state = ibd_states[t]
    if ibd_state:
        start = pos[t]
    
    # Iterate
    bed_records = []
    for t in range(1, len(pos)):
        
        # Handle end of chromosomes
        if chrom != chroms[t]:
            if ibd_state:
                end = pos[t-1]
                bed_records.append(BEDRecord(chrom, start, end, name))
            if ibd_states[t]:
                start = pos[t]
            chrom = chroms[t]
            ibd_state = ibd_states[t]
            continue
            
        # Handle chromosome internal
        if ibd_state:
            if not ibd_states[t]: # segement end
                end = (pos[t-1] + pos[t]) / 2
                bed_records.append(BEDRecord(chrom, start, end, name))
        elif ibd_states[t]: # segment start
            start = (pos[t-1] + pos[t]) / 2
        
        # Update memory
        chrom = chroms[t]
        ibd_state = ibd_states[t]
        
    # Terminate
    if ibd_state:
        end = pos[t-1]  # as with end of chromosome
        bed_records.append(BEDRecord(chrom, start, end, name))
        
    return pd.DataFrame(bed_records)




# ================================================================================
# Plotting
#
# ================================================================================




class ChromosomeInformation:
    
    CHROM_GAP = 10**4
    EDGE_GAP = 10**5
    CMAP = "Paired"  # colors for chromosomes
    
    def __init__(self, site_df: pd.DataFrame):
        """
        Hold information about individual chromosomes;
        including how much to shift sites to produce a genome-wide
        x-axis
        
        """
        
        # Check inputs
        if not "pos" in site_df.columns:
            raise ValueError("Column 'pos' must be in dataframe.")        
        if not "chrom" in site_df.columns:
            raise ValueError("Column 'chrom' must be in dataframe.")
        
        # Store
        self.site_df = site_df
        
        # Compute length information
        self.lengths = self._get_chromosome_lengths()
        self.genome_length = sum(self.lengths.values())
        
        # Chromosome names & colors
        self.names = list(self.lengths.keys())
        self.cols = self._get_chromosome_colors()
        
        # Compute shifts
        self._shift_lookup = self._calc_chromosome_shifts(self.lengths)
        
        
    def _get_chromosome_lengths(self) -> Dict[str, int]:
        """
        Create a dictionary of chromosome lengths based on the passed sites
        
        """
        
        dt = {}
        for chrom_name, chrom_df in self.site_df.groupby("chrom"):
            dt[chrom_name] = chrom_df["pos"].max()
    
        return dt
    
    
    def _get_chromosome_colors(self) -> Dict[set, tuple]:
        """
        Create a dictionary mapping each chrosome to a color
        
        """
        
        # Prepare colors
        assert len(self.names) == 14, "Only support 14 chromosomes."
        dt = dict(zip(
            self.names,
            sns.color_palette(self.CMAP, 10) + sns.color_palette(self.CMAP, 4)
        ))
        
        return dt
    
    
    def _calc_chromosome_shifts(self, chrom_lengths) -> Dict[str, int]:
        """
        Compute shifted site position to make a genome-wide x-axis plot
        from chromosomal positions
        
        """
        
        lookup = {}
        running_shift = self.EDGE_GAP
        for chrom_name, chrom_length in chrom_lengths.items():
            lookup[chrom_name] = running_shift
            running_shift += chrom_length + self.CHROM_GAP
        
        return lookup
    
    
    def get_shifted_positions(self, chroms: List[str], pos: List[int]) -> List[int]:
        """
        Get shifted positions given an array of chromosomes and positions
        
        """
        
        return [
            p + self._shift_lookup[c]
            for (c, p) in zip(chroms, pos)
        ]
    
    
    def set_genome_axis(self, ax):
        """
        Set the limits, ticks, and labels for a genome-wide axis

        """
        # Limits
        xmin = 0 - self.EDGE_GAP
        xmax = self._shift_lookup[self.names[-1]] 
        xmax += self.lengths[self.names[-1]] 
        xmax += self.EDGE_GAP
        ax.set_xlim(xmin, xmax)

        # Ticks and labels
        ax.xaxis.set_major_locator(plt.MultipleLocator(10**6))
        ax.xaxis.set_minor_locator(plt.MultipleLocator(2.5*10**5))
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda k,v: f"{int(k/10**6):00d}"))


class WSAFPlotter:
    def __init__(self, wsaf_df, chrom_info):
        """
        Create a plot of WSAF along the genome
        
        """
        
        # Storage
        self.wsaf_df = wsaf_df
        self.chrom_info = chrom_info
        
        # Compute adjusted positions
        self.wsaf_df["ppos"] = chrom_info.get_shifted_positions(
            chroms=wsaf_df["chrom"],
            pos=wsaf_df["pos"]
        )
        
        
    def plot(self, ax, title=None, add_grid=True):
        """
        Plot WSAF along the genome
        
        """
        
        # Plot
        ax.scatter(
            x=self.wsaf_df["ppos"],
            y=self.wsaf_df["wsafs"],
            c=[self.chrom_info.cols[c] for c in self.wsaf_df["chrom"]],
            s=3
        )
        
        # Ticks
        ax.yaxis.set_major_locator(plt.MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.05))
                
        # Labels
        ax.set_ylabel("WSAF")
        
        # Limits
        ax.set_ylim((-0.02, 1.02))
        
        # Grid
        if add_grid:
            ax.grid(which='major',
                    ls='solid', 
                    zorder=-10,
                    alpha=0.25)
            ax.grid(which='minor', 
                    ls='dotted', 
                    zorder=-10, 
                    alpha=0.25)
            
        # Genome-axis
        self.chrom_info.set_genome_axis(ax)


# TODO: Some clear code duplication here
class IBDSegmentPlotter:

    IBD_PAIR_HEIGHT = 0.5    
    IBD_BAR_HEIGHT = 0.4
    
    def __init__(self, segs_df, chrom_info):
        """
        Plot IBD segments across the genome
        - I probably want to refactor to show all groups, regardless of
        whether they have any IBD segments
        
        """
        
        # Storage
        self.segs_df = segs_df
        self.chrom_info = chrom_info
        
        # Compute adjusted positions
        for c in ["start", "end"]:
            self.segs_df[f"p{c}"] = chrom_info.get_shifted_positions(
                chroms=segs_df["chrom"], pos=segs_df[c]
            )
        
        # Group by name
        self.pairwise_grps = self.segs_df.groupby("name")
        self.n_grps = len(self.pairwise_grps)
        
        # Optional true IBD segmennts
        self.true_segs = False
        
        
    @classmethod
    def from_ibd_states(cls, ibd_df, chrom_info):
        """
        Instantiate the plotter from a dataframe containing IBD state
        information
        
        This will definitely involve assuption
            - mainly about the way thee columns are named
        
        """

        pair_names = [c for c in ibd_df.columns if c.startswith("strain")]
        segs_df = pd.concat(
                    [get_ibd_segments(
                        chroms=ibd_df["chrom"], 
                        pos=ibd_df["pos"], 
                        ibd_states=ibd_df[c],
                        name=c)
                     for c in pair_names]

                )
    
        return cls(segs_df, chrom_info)
    
    
    def load_true_segments(self, true_segs_df):
        """
        If `true` IBD segments are known, load and process
        for inclusion in plot
        
        """
        
        self.true_segs = True
        self.true_segs_df = true_segs_df
        
        # Compute adjusted positions
        for c in ["start", "end"]:
            self.true_segs_df[f"p{c}"] = self.chrom_info.get_shifted_positions(
                chroms=true_segs_df["chrom"], pos=true_segs_df[c]
            )
        
        # Group by name
        self.true_pairwise_grps = self.true_segs_df.groupby("name")
        
    
    def _build_y_axis(self, ax):
        """
        Create y-axis ticks and labels for
        pairwise compairsons
        
        """
        
        # Limits
        min_val = -1
        max_val = len(self.pairwise_grps)
        ax.set_ylim(min_val, max_val)
        
        # Prepare labels
        labels = [""] + list(self.pairwise_grps.groups.keys()) + [""]
        
        # Create ticks
        ax.set_yticks(range(min_val, max_val+1))
        ax.set_yticklabels(labels, fontsize=9)
        
        
    def plot(self, ax, title=None, add_grid=True):
        """
        Make a plot of all IBD segments for each pairwise
        comparison in the segments dataframe
        
        """
        
        # Plot all IBD segments
        ax.plot()
        for j, (name, pair_df) in enumerate(self.pairwise_grps):
            ax.axhline(j, lw=0.5, color='black', zorder=-10)
            for _, seg in pair_df.iterrows():
                rect = Rectangle(
                    xy=(seg["pstart"], j - self.IBD_BAR_HEIGHT/2),
                    height=self.IBD_BAR_HEIGHT,
                    width=seg["pend"] - seg["pstart"],
                    facecolor=self.chrom_info.cols[seg["chrom"]],
                    edgecolor="black",
                    lw=0.5
                )
                ax.add_patch(rect)
            
        # Create y-axis
        self._build_y_axis(ax)
        
        # Genome-axis
        self.chrom_info.set_genome_axis(ax)
        
        # Grid
        if add_grid:
            ax.grid(which='major',
                    ls='solid', 
                    zorder=-10,
                    alpha=0.25)
            ax.grid(which='minor',
                    ls='dotted', 
                    zorder=-10, 
                    alpha=0.25)
            
        # Add overlay if true
        if self.true_segs:
            for j, (name, pair_df) in enumerate(self.true_pairwise_grps):
                ax.axhline(j, lw=0.5, color='black', zorder=-10)
                for _, seg in pair_df.iterrows():
                    rect = Rectangle(
                        xy=(seg["pstart"], j - self.IBD_BAR_HEIGHT/2),
                        height=self.IBD_BAR_HEIGHT,
                        width=seg["pend"] - seg["pstart"],
                        facecolor="none",
                        edgecolor="black",
                        lw=1.5
                    )
                    ax.add_patch(rect)


class CombinedPlotter:
    
    # Constants
    SCALING = 0.08
    ROWS_PER_IBD_SEG = 4
    
    def __init__(self, chrom_info, wsaf_plotter, ibd_plotter):
        self.chrom_info = chrom_info
        self.wsaf_plotter = wsaf_plotter
        self.ibd_plotter = ibd_plotter
        
    def plot(self, title=None, output_path=None):
        """
        Create a combined plot
        
        """
        
        SCALING = 0.15
        ROWS_PER_IBD_SEG = 2
        n_ibd_pairs = self.ibd_plotter.n_grps
        
        # Axes ordering and sizing
        AxesFrame = namedtuple("AxesFrame", ["name", "rows", "plot_func"])
        axes_order = [
            AxesFrame("WSAF", 24, self.wsaf_plotter.plot),
            AxesFrame("IBDSegs", n_ibd_pairs*ROWS_PER_IBD_SEG, self.ibd_plotter.plot),
        ]
        
        # Figue sizes
        total_rows = sum([a.rows for a in axes_order]) + len(axes_order)
        height = total_rows * SCALING
        width = 10
        
        # Create
        fig = plt.figure(figsize=(width, height))
        fig.subplots_adjust(hspace=0.2)
        
        # Grid
        gs = GridSpec(nrows=total_rows, ncols=1)
        
        # Iterate and plot
        l = 0
        for i, axis_item in enumerate(axes_order):
            
            # Plot
            print(f"  {axis_item.name}...")
            ax = plt.subplot(gs[l:(l+axis_item.rows+1)])
            axis_item.plot_func(ax)
            
            # Some bespoke extra stuff
            if axis_item.name == "WSAF":
                ax.set_title(title, loc="left")
                ax.tick_params(labelbottom=False)
                
            if axis_item.name == "IBDSegs":
                ax.set_xlabel("Genome Position (Mbp)")
                
            # Define boundary of next plot
            l = l + axis_item.rows + 1

        # Optionally write
        if output_path is not None:
            fig.savefig(output_path, bbox_inches="tight", pad_inches=0.5)
            plt.close(fig)




# ================================================================================
# Entry
#
# ================================================================================




@click.command(short_help="Make plots of Tapestry outputs.")
@click.option(
    "-i",
    "--input_vcf",
    type=click.Path(exists=True),
    required=True,
    help="Path to input VCF.",
)
def main(input_vcf):

    # Data prepare
    data = SimulatedSequenceData(input_vcf)
    summary_df = pd.read_csv(data.summary_csv)
    true_ibd_segments_df = pd.read_csv(data.segment_csv)

    # Iterate over samples
    sample_names = [s for s in os.listdir(data.results_dir) if s.startswith("SMI")]
    for sample_name in sample_names:
        
        # Create sample directory
        print(f"Plotting sample: {sample_name}...")
        sample_dir = f"{data.results_dir}/{sample_name}"
        
        # Extract parameters
        sample_info = summary_df.query("sample_id == @sample_name").squeeze()
        if (sample_info.shape[0] == 0):
            print("  Couldn't find truth. Skipping.")
            continue

        # Prepare title
        props = ", ".join([f"{float(p):.03f}" for p in sample_info["props"].split(";")])

        title = f"{sample_name} | "
        title += props
        #title += ", ".join([f"{sample_info[p]:.03f}" for p in sample_info.keys() if p.startswith("prop")])
        title += " | $f_{IBD}=$" + f"{sample_info['f_ibd']:.03}"

        # LOAD IBD
        # Pairwise matrix
        ibd_pairwise_df = pd.read_csv(f"{sample_dir}/fit.ibd.pairwise.csv")
        # Segments
        try:
            ibd_segments_df = pd.read_csv(f"{sample_dir}/fit.ibd.segments.bed",
                                        sep="\t",
                                        header=None)
            ibd_segments_df.columns = ["chrom", "start", "end", "name"]

            true_sample_segments_df = true_ibd_segments_df.query("sample_id == @sample_name")
            true_sample_segments_df.insert(0, "name", 
                                           [f"strains{row['strain1']}-{row['strain2']}"
                                            for _, row in true_sample_segments_df.iterrows()])
        except pd.errors.EmptyDataError:  # no segments
            ibd_segments_df = pd.DataFrame(columns=["chrom", "start", "end", "name"])
        

        
        # Plot
        chrom_info = ChromosomeInformation(ibd_pairwise_df)
        wsaf_plotter = WSAFPlotter(ibd_pairwise_df, chrom_info)
        ibd_plotter = IBDSegmentPlotter(ibd_segments_df, chrom_info)
        ibd_plotter.load_true_segments(true_sample_segments_df)
        comb_plotter = CombinedPlotter(chrom_info, wsaf_plotter, ibd_plotter)
        comb_plotter.plot(title=title,
                          output_path=f"{sample_dir}/plot.ibd.pairwise.{sample_name}.pdf")


if __name__ == "__main__":
    main()

