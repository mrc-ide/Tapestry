import pandas as pd
from dataclasses import dataclass
from typing import List, Optional


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
