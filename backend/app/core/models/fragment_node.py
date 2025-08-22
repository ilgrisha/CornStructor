# File: backend/app/core/models/fragment_node.py
# Version: v0.2.0

"""
Shared dataclass representing a fragment in the assembly tree.
"""

from dataclasses import dataclass, field
from typing import List


@dataclass
class FragmentNode:
    """
    Represents one node in the assembly tree.
    """
    fragment_id:   str                # unique ID for locating this fragment (uniform format)
    level:         int
    start:         int                # absolute start in the root sequence
    end:           int                # absolute end in the root sequence
    seq:           str                # fragment sequence (body+overlaps)
    strand:        str                # '+' or '-' orientation
    overlap_prev:  str                # upstream overlap in node orientation
    overlap_next:  str                # downstream overlap
    is_oligo:      bool               # True if this is a leaf (no further subdivision)
    ga_log:        List[float]        = field(default_factory=list)
    ga_detail:     List[dict]         = field(default_factory=list)  # per-generation stats: {gen,best,mean,std}
    ga_cluster_id: str                = ""  # stable id to group siblings' GA run (unique per parent cluster)
    children:      List["FragmentNode"] = field(default_factory=list)
