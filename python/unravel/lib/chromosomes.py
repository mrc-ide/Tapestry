import pandas as pd
from typing import Dict, List

import seaborn as sns
import matplotlib.pyplot as plt


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
