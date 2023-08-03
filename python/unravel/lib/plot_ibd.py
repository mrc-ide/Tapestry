import pandas as pd

from collections import namedtuple

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

from unravel.lib.ibd import get_ibd_segments

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
