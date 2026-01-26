# base functions

# -------------------- Loading Stuff -------------------- #

using CodecZlib
using MatrixMarket
using SparseArrays
using CSV
using DataFrames

# rows = genes
# cols = cells


function gunzip_file(gz_path::AbstractString)

    @assert endswith(gz_path, ".gz")
    new_path = replace(gz_path, ".gz" => "")

    open(gz_path) do fin
        open(new_path, write=true) do fout
            write(fout, GzipDecompressorStream(fin))
        end
    end

    return new_path

end

function read_mtx(path::AbstractString)

    @assert endswith(path, ".mtx")

    return sparse(MatrixMarket.mmread(path))

end

function read_tsv(path::AbstractString)

    @assert endswith(path, ".tsv")

    return CSV.read(path, DataFrame; delim='\t', header=false)

end


function load_10x()

    # load features, counts and metadata
    # initialize matrices for raw counts, normalize/log1p counts, and scaled counts
    # also for PCA and then a graph object?

end

# -------------------- QC and Filtering -------------------- #

mito_prefix = na
ribo_prefix = na

function qc()
    # add n_counts, n_features/n_genes, percent_mito, percent_ribo
end

function filter_cells()

    # filter min count per cell
    # min number of features/genes
    # max number of features (if mentioned)
    # max percent mito

end

# remove genes from count matrix that have less than the number min_count reads in cell_cutoff number of cells
# aka if cell_cutoff = 3 and min_count is 1, only genes with a count of at least 1 in at least 3 cells will be kept.
function filter_genes(gene_counts::AbstractMatrix, cell_cutoff::Int, expression_cutoff::Int)
    n_genes = size(gene_counts, 1)
    keep = falses(n_genes)

    for g in 1:n_genes
        low_cells = 0

        for c in 1:size(gene_counts, 2)
            if gene_counts[g, c] < expression_cutoff
                low_cells += 1
            end
        end

        keep[g] = low_cells < cell_cutoff
    end

    return gene_counts[keep, :]
end

# standard normalization
function normalize(gene_counts::AbstractMatrix, scaling_factor::Int64 = 10000)

    cell_totals = sum(gene_counts, dims = 1) # sum columns = cells

    cell_totals[cell_totals .== 0] .= 1

    normalized_counts = (gene_counts ./ cell_totals) .* (scaling_factor)

    return normalized_counts

end

# normal log1p transform
function log1p_transform(gene_counts::AbstractMatrix)

    logged_counts = log1p.(gene_counts)

    return logged_counts

end

# log (x + E) transform
function log_e_transform(gene_counts::AbstractMatrix, epsilon::Int64)

    logged_counts = log.(epsilon .+ gene_counts) # for custom log transform

    return logged_counts

end

safe_sqrt(x) = sqrt(x < 0 ? 0.0 : x)

# finds the top no_variable_genes number of variable genes
function find_variable_genes(gene_counts::AbstractMatrix, no_variable_genes::Int)

    n_genes = size(gene_counts, 1)

    gene_variances = Vector{Float64}(undef, no_variable_genes)

    for g in 1:n_genes

        gene_variance[g] = var(gene_counts[g, :])

    end

    no_keep = min(n_hvg, n_genes)
    sorted_indices = sortperm(gene_variances; rev = true)
    variable_genes = sorted_indices[1:no_keep]

    return variable_genes, gene_variances

end

function scale_data()

end

# deviance residual method
function deviance_residuals(counts::AbstractMatrix; epsilon::Float64 = 1e-12)

    cell_library_sizes = vec(sum(counts, dims = 1))
    mean_library_size  = mean(cell_library_sizes)

    size_factors = cell_library_sizes ./ mean_library_size

    for i in eachindex(size_factors)
        if size_factors[i] == 0.0
            size_factors[i] = 1.0
        end
    end

    depth_corrected_counts = counts ./ reshape(size_factors, 1, :)

    gene_baseline_means = mean(depth_corrected_counts, dims = 2)

    for i in eachindex(gene_baseline_means)
        if gene_baseline_means[i] == 0.0
            gene_baseline_means[i] = epsilon
        end
    end

    expected_counts = gene_baseline_means .* reshape(size_factors, 1, :)

    residuals = similar(expected_counts, Float64)

    for g in axes(counts, 1), c in axes(counts, 2)
        observed = counts[g, c]
        expected = expected_counts[g, c]

        if observed == 0.0
            deviance = 2.0 * expected
        else
            deviance = 2.0 * (
                observed * log(observed / expected) -
                (observed - expected)
            )
        end

        # numerical guard
        if deviance < 0.0
            deviance = 0.0
        end

        residuals[g, c] = sign(observed - expected) * sqrt(deviance)
    end

    return residuals

end

function pca()
end

# maybe import harmony

function find_neighbors()
end

function make_graph()
end

function leiden_cluster()
end

function louvain_cluster()
end

# plotting


d = read_10x_mtx("filtered_feature_bc_matrix/")
qc!(d)
filter_cells!(d; min_counts=500, min_genes=200, max_pct_mito=20)
filter_genes!(d; min_cells=3)

normalize_total!(d)
log1p!(d)

select_hvgs!(d; n=2000)
scale!(d)

pca!(d; n_pcs=50)
neighbors!(d; k=15, n_pcs=30)

leiden!(d; resolution=1.0)
