# Network and cache helpers for Gviz plotting (sourced by locus_summary_plots_gviz.R).
# locusGviz() calls BiomartGeneRegionTrack() and IdeogramTrack() unqualified, so the
# cached versions defined here replace the Gviz constructors for that script only.
# Set options(targpred.gviz_cache_dir=<dir>) before plotting.

gviz_cache_dir <- function() {
    dir <- getOption("targpred.gviz_cache_dir")
    if (is.null(dir)) stop("options(targpred.gviz_cache_dir) is not set")
    dir.create(dir, recursive=TRUE, showWarnings=FALSE)
    dir
}

# Writes via a temporary file and rename, so concurrent workers never read a
# partially written cache file.
save_cache_rds <- function(object, path) {
    temp <- paste0(path, ".", Sys.getpid())
    saveRDS(object, temp)
    stopifnot(file.rename(temp, path))
}

# biomaRt's default query cache is ~/.cache/biomaRt. On a network file system,
# concurrent processes serialise on its sqlite lock. Each process gets its own
# local cache; gene tracks are cached persistently in gviz_cache_dir().
use_local_biomart_cache <- function() {
    Sys.setenv(BIOMART_CACHE=file.path(tempdir(), paste0("biomaRt_", Sys.getpid())))
}

# biomaRt 2.54 hard-codes a 10-second timeout in bmRequest. Allow the same
# timeout as gene queries; a no-op for versions without the hard-coded value.
extend_biomart_timeout <- function(seconds=180) {
    options(timeout=seconds)
    bm_request <- getFromNamespace("bmRequest", "biomaRt")
    body(bm_request) <- parse(text=gsub("timeout(10)", paste0("timeout(", seconds, ")"),
                                      paste(deparse(body(bm_request)), collapse="\n"),
                                      fixed=TRUE))[[1]]
    assignInNamespace("bmRequest", bm_request, ns="biomaRt")
}

# Ensembl intermittently answers gene queries with server-side failures
# (mostly 504 gateway timeouts). Gviz re-queries Ensembl at plot time whenever
# the plot window extends past the cached gene track, so callers retry the whole
# plot. Other errors are not retried.
transient_error <- paste("server error", "HTTP error code", "Timeout was reached",
                         "timed out", "Could not resolve host", "Connection reset",
                         sep="|")
with_retry <- function(f, attempts=5, wait=60) {
    for (attempt in seq_len(attempts)) {
        result <- tryCatch(f(), error=function(e) e)
        if (!inherits(result, "error")) return(result)
        if (attempt == attempts || !grepl(transient_error, conditionMessage(result))) stop(result)
        delay <- wait * 2^(attempt - 1)
        message("Transient Ensembl error (attempt ", attempt, "/", attempts, "), retrying in ",
                delay, " s: ", conditionMessage(result))
        Sys.sleep(delay)
    }
}

# The GRCh37 Mart is connected lazily, so fully cached runs need no Ensembl access.
gene_mart <- NULL
get_gene_mart <- function() {
    if (!is.null(gene_mart)) return(gene_mart)
    path <- file.path(gviz_cache_dir(), "biomart_grch37.rds")
    if (file.exists(path)) {
        cached <- readRDS(path)
        stopifnot(grepl("GRCh37", cached$assembly))
        gene_mart <<- cached$mart
        return(gene_mart)
    }
    for (attempt in seq_len(5)) {
        gene_mart <<- tryCatch(useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl",
                                         host="https://grch37.ensembl.org"),
                               error=function(e) { message("Attempt ", attempt, " failed: ", conditionMessage(e)); NULL })
        if (!is.null(gene_mart)) break
    }
    if (is.null(gene_mart)) stop("Failed to connect to the GRCh37 Ensembl gene service after 5 attempts")
    datasets <- listDatasets(gene_mart)
    assembly <- datasets$version[datasets$dataset == "hsapiens_gene_ensembl"]
    if (length(assembly) != 1L || !grepl("GRCh37", assembly)) {
        stop("BioMart gene dataset is not GRCh37: ", paste(assembly, collapse=", "))
    }
    save_cache_rds(list(mart=gene_mart, assembly=assembly, retrieved=Sys.time(),
                        host="https://grch37.ensembl.org"), path)
    message("Cached BioMart metadata for ", assembly)
    gene_mart
}

# One gene track per GRB region, reused across traits and reruns.
BiomartGeneRegionTrack <- function(...) {
    dots <- list(...)
    key <- paste(dots$chromosome, dots$start, dots$end, sep="_")
    path <- file.path(gviz_cache_dir(), paste0("genes_", key, ".rds"))
    if (file.exists(path)) return(readRDS(path))
    dots$biomart <- get_gene_mart()
    result <- do.call(Gviz::BiomartGeneRegionTrack, dots)
    save_cache_rds(result, path)
    result
}

# Ideograms from local hg19 cytobands instead of a live UCSC request per plot.
IdeogramTrack <- function(...) {
    dots <- list(...)
    path <- file.path(gviz_cache_dir(), paste0("ideogram_", dots$chromosome, ".rds"))
    if (file.exists(path)) return(readRDS(path))
    bands <- file.path(gviz_cache_dir(), "hg19_cytoBandIdeo.txt.gz")
    if (!file.exists(bands)) {
        temp <- paste0(bands, ".", Sys.getpid())
        download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBandIdeo.txt.gz",
                      temp, mode="wb", quiet=TRUE)
        stopifnot(file.rename(temp, bands))
    }
    dots$bands <- read.delim(bands, header=FALSE,
                             col.names=c("chrom", "chromStart", "chromEnd", "name", "gieStain"))
    result <- do.call(Gviz::IdeogramTrack, dots)
    save_cache_rds(result, path)
    result
}
