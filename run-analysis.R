markdown.dir <- file.path("markdown")
if(!dir.exists(markdown.dir)){
  dir.create(markdown.dir)
}

# Kraken2:
knitr::spin("code/r-scripts/kraken2-pipeline/001-kraken2-classification-stats.R",
              knit = FALSE)
file.rename("code/r-scripts/kraken2-pipeline/001-kraken2-classification-stats.Rmd",
            file.path(markdown.dir, "001-kraken2-classification-stats.Rmd"))
rmarkdown::render(file.path(markdown.dir, "001-kraken2-classification-stats.Rmd"),
                  'html_document',
                  knit_root_dir = here::here())

knitr::spin("code/r-scripts/kraken2-pipeline/002-phyloseq-kraken2.R", knit = FALSE)
file.rename("code/r-scripts/kraken2-pipeline/002-phyloseq-kraken2.Rmd",
            file.path(markdown.dir, "002-phyloseq-kraken2.Rmd"))
rmarkdown::render(file.path(markdown.dir, "002-phyloseq-kraken2.Rmd"),
                  'html_document',
                  knit_root_dir = here::here())

knitr::spin("code/r-scripts/kraken2-pipeline/003-summary-stats-kraken2.R",
            knit = FALSE)
file.rename("code/r-scripts/kraken2-pipeline/003-summary-stats-kraken2.Rmd",
            file.path(markdown.dir, "003-summary-stats-kraken2.Rmd"))
rmarkdown::render(file.path(markdown.dir, "003-summary-stats-kraken2.Rmd"),
                  'html_document',
                  knit_root_dir = here::here())

knitr::spin("code/r-scripts/kraken2-pipeline/004-barplots-kraken2.R", knit = FALSE)
file.rename("code/r-scripts/kraken2-pipeline/004-barplots-kraken2.Rmd",
            file.path(markdown.dir, "004-barplots-kraken2.Rmd"))
rmarkdown::render(file.path(markdown.dir, "004-barplots-kraken2.Rmd"),
                  'html_document',
                  knit_root_dir = here::here())

knitr::spin("code/r-scripts/kraken2-pipeline/005-kraken2-maaslin2.R", knit = FALSE)
file.rename("code/r-scripts/kraken2-pipeline/005-kraken2-maaslin2.Rmd",
            file.path(markdown.dir, "005-kraken2-maaslin2.Rmd"))
rmarkdown::render(file.path(markdown.dir, "005-kraken2-maaslin2.Rmd"),
                  'html_document',
                  knit_root_dir = here::here())

# MAG:
knitr::spin("code/r-scripts/mag-assembly/001-mag-combine-all-data.R",
              knit = FALSE)
file.rename("code/r-scripts/mag-assembly/001-mag-combine-all-data.Rmd",
            file.path(markdown.dir, "006-mag-combine-all-data.Rmd"))
rmarkdown::render(file.path(markdown.dir, "006-mag-combine-all-data.Rmd"),
                  'html_document',
                  knit_root_dir = here::here())

knitr::spin("code/r-scripts/mag-assembly/002-summary-stats-mag.R",
            knit = FALSE)
file.rename("code/r-scripts/mag-assembly/002-summary-stats-mag.Rmd",
            file.path(markdown.dir, "007-summary-stats-mag.Rmd"))
rmarkdown::render(file.path(markdown.dir, "007-summary-stats-mag.Rmd"),
                  'html_document',
                  knit_root_dir = here::here())

knitr::spin("code/r-scripts/mag-assembly/003-mag-taxonomy-plots.R",
            knit = FALSE)
file.rename("code/r-scripts/mag-assembly/003-mag-taxonomy-plots.Rmd",
            file.path(markdown.dir, "008-mag-taxonomy-plots.Rmd"))
rmarkdown::render(file.path(markdown.dir, "008-mag-taxonomy-plots.Rmd"),
                  'html_document',
                  knit_root_dir = here::here())

knitr::spin("code/r-scripts/mag-assembly/004-mag-gene-annotation-combine-data.R",
            knit = FALSE)
file.rename("code/r-scripts/mag-assembly/004-mag-gene-annotation-combine-data.Rmd",
            file.path(markdown.dir, "009-mag-gene-annotation-combine-data.Rmd"))
rmarkdown::render(file.path(markdown.dir, "009-mag-gene-annotation-combine-data.Rmd"),
                  'html_document',
                  knit_root_dir = here::here())

knitr::spin("code/r-scripts/mag-assembly/005-mag-gene-annotation-plots.R",
            knit = FALSE)
file.rename("code/r-scripts/mag-assembly/005-mag-gene-annotation-plots.Rmd",
            file.path(markdown.dir, "010-mag-gene-annotation-plots.Rmd"))
rmarkdown::render(file.path(markdown.dir, "010-mag-gene-annotation-plots.Rmd"),
                  'html_document',
                  knit_root_dir = here::here())

library(bookdown)
render_book(input = file.path(markdown.dir),  
            output_format = "bookdown::html_document2",
            clean = TRUE,
            envir = parent.frame(),
            output_dir = file.path( "Metagenome_bind") , 
            new_session = TRUE, preview = FALSE,
            config_file = file.path(markdown.dir,"_bookdown.yml"))
