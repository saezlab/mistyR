FROM rocker/r-base:latest

WORKDIR /home/docker
RUN mkdir mistyR
COPY . mistyR
WORKDIR mistyR

RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y curl
RUN Rscript -e 'install.packages("remotes")'
RUN $(Rscript -e 'cat("apt-get install -y", paste(gsub("apt-get install -y ", "", remotes::system_requirements("ubuntu", "20.04"))))')
RUN apt-get install -y ghostscript libmagick++-dev libhdf5-dev libgeos-c1v5
RUN Rscript -e 'remotes::install_cran("BiocManager", ask = FALSE); BiocManager::install(version = "3.15", ask = FALSE)'
RUN Rscript -e 'remotes::install_deps(dependencies = TRUE, repos = BiocManager::repositories(), upgrade = "always")'
RUN Rscript -e 'remotes::install_cran(c("Seurat", "rhdf5", "tinytex"), ask = FALSE)'
RUN Rscript -e 'tinytex::install_tinytex(); tinytex::tlmgr_update(); tinytex::tlmgr_install("pdfcrop")'

RUN Rscript -e 'remotes::install_local(build_vignettes = TRUE)'

WORKDIR /home/docker
RUN rm -rf mistyR
