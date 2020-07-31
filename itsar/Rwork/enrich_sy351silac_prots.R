library(SY351SILAC)
library(enrichR)
library(tidyverse)


dbs <- enrichR::listEnrichrDbs()
enriched <- enrichr(genes =  psitestats.protein.totals %>% select(Gene) %>% pull(), dbs$libraryName)
numsig <- unlist(lapply(enriched, function(x) sum(x$Adjusted.P.value < 0.01))) %>% sort()

endf <- dplyr::bind_rows(enriched, .id = "library_name")
endf <- dplyr::as_tibble(endf)

save(endf, file = 'endf.RData')

endf <- endf %>% dplyr::arrange(Adjusted.P.value)
