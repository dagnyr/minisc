# base functions

# -------------------- Loading Stuff -------------------- #

using CodecZlib
using MatrixMarket
using SparseArrays
using CSV
using DataFrames

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
function filter_genes(gene_counts::AbstractMatrix, cell_cutoff::Int, min_count::Int)
    filtered = vec(sum(gene_counts .< cell_cutoff, dims = 1) .< min_count)
    return X[:, filtered]
end

function normalize()
end

function log1p()
end

function find_variable_genes()
end

function scale_data()
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
