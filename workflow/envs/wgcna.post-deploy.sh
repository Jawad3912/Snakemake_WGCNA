#!/usr/bin/env bash
set -euo pipefail

fonts="${CONDA_PREFIX:?}/fonts/wgcna"
mkdir -p "$fonts"
curl -fsSL --retry 3 https://github.com/ArtifexSoftware/urw-base35-fonts/archive/refs/tags/20200910.tar.gz |
  tar xz -C "$fonts" --strip-components=2 --wildcards '*/fonts/NimbusSans-*.otf'
curl -fsSL --retry 3 https://github.com/dejavu-fonts/dejavu-fonts/releases/download/version_2_37/dejavu-fonts-ttf-2.37.tar.bz2 |
  tar xj -C "$fonts" --strip-components=2 --wildcards '*/ttf/DejaVuSans.ttf' '*/ttf/DejaVuSans-Bold*.ttf' '*/ttf/DejaVuSans-Oblique.ttf'
(cd "$fonts" && sha256sum -c --quiet) <<'EOF'
7c25be4d78155523080ab85b10277150657ff7dabbcad7037bdd536c9b6d0d08  NimbusSans-Regular.otf
7f33328e6b4d4cd21b45fa625791928c9407dc702db6780e56b09ca9a3ecaa67  NimbusSans-Bold.otf
7b0bef5686aa58c0fd0f0d01beeae56664208490e30ca6a25431281c9a0c6402  NimbusSans-Italic.otf
3f47fb34fcb7de09f8cbc9f305191340ddebf7a068419f4bb5f49287dea59b87  NimbusSans-BoldItalic.otf
7da195a74c55bef988d0d48f9508bd5d849425c1770dba5d7bfc6ce9ed848954  DejaVuSans.ttf
e6476c1b80502924294eed40894c5b18e06c181444ca953e5334262df9c27724  DejaVuSans-Bold.ttf
4af75fa16ee6d3ad43e1ecec41862c24954af26a55c6bb1ebb27bd486a50f5f4  DejaVuSans-Oblique.ttf
eb436dca0c2594b73d8b603b892e374fdfd8d885d25ffb4f18df4c4c0b49e50f  DejaVuSans-BoldOblique.ttf
EOF

alias_sans() {
  printf '  <alias binding="strong">\n    <family>%s</family>\n' "$1"
  printf '    <prefer><family>Nimbus Sans</family><family>DejaVu Sans</family></prefer>\n  </alias>\n'
}
{
  printf '<?xml version="1.0"?>\n<!DOCTYPE fontconfig SYSTEM "urn:fontconfig:fonts.dtd">\n<fontconfig>\n'
  printf '  <reset-dirs/>\n  <dir>%s</dir>\n' "$fonts"
  alias_sans sans-serif
  alias_sans Helvetica
  alias_sans Arial
  printf '</fontconfig>\n'
} >"$CONDA_PREFIX/etc/fonts/conf.d/01-wgcna-fonts.conf"
fc-cache -f "$fonts" >/dev/null

R --vanilla --no-echo <<'EOF'
cran_snapshot <- "https://packagemanager.posit.co/cran/2026-11-16"
bioc <- "https://bioconductor.org/packages/3.23"

repos <- c(
  CRAN = cran_snapshot,
  BioCsoft = file.path(bioc, "bioc"),
  BioCann = file.path(bioc, "data/annotation")
)

pinned <- data.frame(
  package = c("DESeq2", "WGCNA", "metafor", "gprofiler2"),
  version = c("1.52.0", "1.74", "5.2-1", "0.2.4"),
  repo = c(repos[["BioCsoft"]], cran_snapshot, cran_snapshot, cran_snapshot)
)

packages <- c(
  pinned$package, "BiocParallel", "tidyverse", "patchwork", "pheatmap",
  "ggrepel", "png", "openxlsx"
)

options(Ncpus = max(1L, parallel::detectCores() - 1L), timeout = 600)
lib <- .libPaths()[1]

is_pinned_version <- function(pkg, version) {
  installed <- installed_version(pkg)
  !is.na(installed) && package_version(installed) == package_version(version)
}

installed_version <- function(pkg) {
  if (requireNamespace(pkg, lib.loc = lib, quietly = TRUE)) {
    as.character(packageVersion(pkg, lib.loc = lib))
  } else {
    NA_character_
  }
}

install_from_archive <- function(package, version, repo) {
  tarball <- sprintf("%s_%s.tar.gz", package, version)
  urls <- file.path(repo, "src/contrib", c(tarball, file.path("Archive", package, tarball)))
  for (url in urls) {
    try(install.packages(url, lib = lib, repos = NULL, type = "source"))
    if (is_pinned_version(package, version)) {
      return(invisible(TRUE))
    }
  }
}

for (attempt in 1:3) {
  todo <- packages[is.na(vapply(packages, installed_version, character(1)))]
  if (!length(todo)) break
  install.packages(todo, lib = lib, repos = repos, type = "source")
}

pinned_ok <- function() mapply(is_pinned_version, pinned$package, pinned$version)

for (i in which(!pinned_ok())) {
  install_from_archive(pinned$package[i], pinned$version[i], pinned$repo[i])
}

versions <- vapply(packages, installed_version, character(1))
missing <- packages[is.na(versions)]
wrong <- pinned$package[!pinned_ok()]
if (length(missing)) stop("Failed to install: ", paste(missing, collapse = ", "))
if (length(wrong)) stop("Wrong version installed for: ", paste(wrong, collapse = ", "))

print(versions)
EOF
