DEFAULT_TOTAL_COLS = [f"group{i}" for i in range(1, 7)]
DEFAULT_PHOSPHO_COLS = [f"p_group{i}" for i in range(1, 7)]
DEFAULT_CORRECTED_COLS = [f"phospho_corrected_{i}" for i in range(1, 7)]

DEFAULT_COMPARISONS = [
    ("group1", "group4", "2277.IgM.5min", "2277.IgM.60min"),
    ("group2", "group5", "2277.IgM.cAMP.5min", "2277.IgM.cAMP.60min"),
    ("group3", "group6", "2282.IgM.5min", "2282.IgM.60min"),
    ("group1", "group2", "2277.IgM.5min", "2277.IgM.cAMP.5min"),
    ("group1", "group3", "2277.IgM.5min", "2282.IgM.5min"),
    ("group2", "group3", "2277.IgM.cAMP.5min", "2282.IgM.5min"),
    ("group4", "group5", "2277.IgM.60min", "2277.IgM.cAMP.60min"),
    ("group4", "group6", "2277.IgM.60min", "2282.IgM.60min"),
    ("group5", "group6", "2277.IgM.cAMP.60min", "2282.IgM.60min"),
]
