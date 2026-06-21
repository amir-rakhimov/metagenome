# Kraken2:
knitr::spin("code/r-scripts/kraken2-pipeline/001-kraken2-classification-stats.R",
              knit = FALSE)
file.rename("code/r-scripts/kraken2-pipeline/001-kraken2-classification-stats.Rmd",
            "markdown/001-kraken2-classification-stats.Rmd")

knitr::spin("code/r-scripts/kraken2-pipeline/002-phyloseq-kraken2.R", knit = FALSE)
file.rename("code/r-scripts/kraken2-pipeline/002-phyloseq-kraken2.Rmd",
            "markdown/002-phyloseq-kraken2.Rmd")

knitr::spin("code/r-scripts/kraken2-pipeline/003-summary-stats-kraken2.R",
            knit = FALSE)
file.rename("code/r-scripts/kraken2-pipeline/003-summary-stats-kraken2.Rmd",
            "markdown/003-summary-stats-kraken2.Rmd")

knitr::spin("code/r-scripts/kraken2-pipeline/004-barplots-kraken2.R", knit = FALSE)
file.rename("code/r-scripts/kraken2-pipeline/004-barplots-kraken2.Rmd",
            "markdown/004-barplots-kraken2.Rmd")

knitr::spin("code/r-scripts/kraken2-pipeline/005-kraken2-maaslin2.R", knit = FALSE)
file.rename("code/r-scripts/kraken2-pipeline/005-kraken2-maaslin2.Rmd",
            "markdown/005-kraken2-maaslin2.Rmd")

rmarkdown::render('./markdown/001-kraken2-classification-stats.Rmd',
                  'html_document',
                  here::here())

rmarkdown::render('./markdown/002-phyloseq-kraken2.Rmd',
                  'html_document',
                  here::here())

rmarkdown::render('./markdown/003-summary-stats-kraken2.Rmd',
                  'html_document',
                  here::here())

rmarkdown::render('./markdown/004-barplots-kraken2.Rmd',
                  'html_document',
                  here::here())

rmarkdown::render('./markdown/005-kraken2-maaslin2.Rmd',
                  'html_document',
                  here::here())



# MAG:
  
knitr::spin("code/r-scripts/mag-assembly/001-mag-combine-all-data.R",
              knit = FALSE)
file.rename("code/r-scripts/mag-assembly/001-mag-combine-all-data.Rmd",
            "markdown/006-mag-combine-all-data.Rmd")

knitr::spin("code/r-scripts/mag-assembly/002-summary-stats-mag.R",
            knit = FALSE)
file.rename("code/r-scripts/mag-assembly/002-summary-stats-mag.Rmd",
            "markdown/007-summary-stats-mag.Rmd")

knitr::spin("code/r-scripts/mag-assembly/003-mag-taxonomy-plots.R",
            knit = FALSE)
file.rename("code/r-scripts/mag-assembly/003-mag-taxonomy-plots.Rmd",
            "markdown/008-mag-taxonomy-plots.Rmd")

knitr::spin("code/r-scripts/mag-assembly/004-mag-gene-annotation-combine-data.R",
            knit = FALSE)
file.rename("code/r-scripts/mag-assembly/004-mag-gene-annotation-combine-data.Rmd",
            "markdown/009-mag-gene-annotation-combine-data.Rmd")

knitr::spin("code/r-scripts/mag-assembly/005-mag-gene-annotation-plots.R",
            knit = FALSE)
file.rename("code/r-scripts/mag-assembly/005-mag-gene-annotation-plots.Rmd",
            "markdown/010-mag-gene-annotation-plots.Rmd")


rmarkdown::render('./markdown/006-mag-combine-all-data.Rmd',
                  'html_document',
                  here::here())

rmarkdown::render('./markdown/007-summary-stats-mag.Rmd',
                  'html_document',
                  here::here())

rmarkdown::render('./markdown/008-mag-taxonomy-plots.Rmd',
                  'html_document',
                  here::here())

rmarkdown::render('./markdown/009-mag-gene-annotation-combine-data.Rmd',
                  'html_document',
                  here::here())

rmarkdown::render('./markdown/010-mag-gene-annotation-plots.Rmd',
                  'html_document',
                  here::here())

render_book(input = "markdown", 
            output_format = "bookdown::html_document2",
            clean = TRUE,
            envir = parent.frame(),
            output_dir = "Metagenome_bind", 
            new_session = TRUE, preview = FALSE,
            config_file = "_bookdown.yml")
