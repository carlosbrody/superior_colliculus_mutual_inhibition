include("svd_cluster.jl")
include("rule_encoding.jl")

data, filenames, encoding, error_types = SVD_interactive("C17"; farmdir="MiniOptimized", backend_mode=true);

# Control, Pro/Anti Index, SVD Dim 1
scatter_svd_by_index(data,encoding, error_types, 1, 1, 1)

PLOT_A_ZILLION_INDEXES= false;

# Just here for convience 
if PLOT_A_ZILLION_INDEXES
# Control, Pro/Anti Index, SVD Dim 2
scatter_svd_by_index(data,encoding, error_types, 1, 1, 2)

# Control, Right/Left Index, SVD Dim 1
scatter_svd_by_index(data,encoding, error_types, 1, 2, 1)

# Control, Right/Left Index, SVD Dim 2
scatter_svd_by_index(data,encoding, error_types, 1, 2, 2)

# Control, Diagonal Index, SVD Dim 1
scatter_svd_by_index(data,encoding, error_types, 1, 3, 1)

# Control, Diagonal Index, SVD Dim 2
scatter_svd_by_index(data,encoding, error_types, 1, 3, 2)


# Delay, Pro/Anti Index, SVD Dim 1
scatter_svd_by_index(data,encoding, error_types, 2, 1, 1)

# Delay, Pro/Anti Index, SVD Dim 2
scatter_svd_by_index(data,encoding, error_types, 2, 1, 2)

# Delay, Right/Left Index, SVD Dim 1
scatter_svd_by_index(data,encoding, error_types, 2, 2, 1)

# Delay, Right/Left Index, SVD Dim 2
scatter_svd_by_index(data,encoding, error_types, 2, 2, 2)

# Delay, Diagonal Index, SVD Dim 1
scatter_svd_by_index(data,encoding, error_types, 2, 3, 1)

# Delay, Diagonal Index, SVD Dim 2
scatter_svd_by_index(data,encoding, error_types, 2, 3, 2)







end
