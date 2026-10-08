// Expected: accept. The author explicitly bounds both table indices.
i = ((+(1) : %(200)) ~ _) : -(50);
bounded = max(0, min(i, 99));
process = rwtable(100, 0.0, bounded, _, bounded);
