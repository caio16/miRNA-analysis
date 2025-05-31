# Análise de miRNAs
O presente projeto visa realizar algumas análises envolvendo miRNAs: significância (ANOVA, Limma, t-test, wilcox-test), target genes e enriquecimento. O script pode ser aplicado para qualquer dataset que contenha valores de expressão de miRNAs em 2 ou mais grupos.  
Todas as análises foram realizadas em R [\[1\]](#r) através do RStudio [\[2\]](#rstudio).

# Instalando Pacotes
Antes de tudo, lembre-se de instalar os pacotes necessários através do `install.packages("pacote")`: readxl [\[3\]](#readxl), dplyr [\[4\]](#dplyr), tidyr [\[5\]](#tidyr), broom [\[6\]](#broom), limma [\[7\]](#limma), pheatmap [\[8\]](#pheatmap), purrr [\[9\]](#purrr), tibble [\[10\]](#tibble), knitr [\[11\]](#knitr), ggplot2 [\[12\]](#ggplot2) e enrichR [\[13\]](#enrichr).  

# Comentários
Todas as etapas necessárias estão descritas ao longo do código. Algumas seções estão inteiramente como comentário por um dos dois motivos: o trecho é opcional ou o teste realizado não produziu resultados significativos.

# Referências
<a id="r"></a>
[1] R Core Team (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.  

<a id="rstudio"></a>
[2] Posit team (2025). RStudio: Integrated Development Environment for R. Posit Software, PBC, Boston, MA. URL http://www.posit.co/.  

<a id="readxl"></a>
[4] Wickham H, Bryan J (2025). readxl: Read Excel Files. R package version 1.4.5, https://github.com/tidyverse/readxl, https://readxl.tidyverse.org.  

<a id="tidyr"></a>
[5] Wickham H, Vaughan D, Girlich M (2024). tidyr: Tidy Messy Data. R package version 1.3.1, https://github.com/tidyverse/tidyr, https://tidyr.tidyverse.org.  

<a id="broom"></a>
[6] Robinson D, Hayes A, Couch S (2025). broom: Convert Statistical Objects into Tidy Tibbles. R package version 1.0.8, https://github.com/tidymodels/broom, https://broom.tidymodels.org/.

<a id="limma"></a>
[7] Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7), e47. doi:10.1093/nar/gkv007.

<a id="pheatmap"></a>
[8] Kolde R (2018). pheatmap: Pretty Heatmaps. R package version 1.0.12, https://github.com/raivokolde/pheatmap.

<a id="purrr"></a>
[9] Wickham H, Henry L (2025). purrr: Functional Programming Tools. R package version 1.0.4, https://github.com/tidyverse/purrr, https://purrr.tidyverse.org/.

<a id="tibble"></a>
[10] Müller K, Wickham H (2024). tibble: Simple Data Frames. R package version 3.2.1, https://github.com/tidyverse/tibble, https://tibble.tidyverse.org/.

<a id="knitr"></a>
[11] Xie Y (2025). knitr: A General-Purpose Package for Dynamic Report Generation in R. R package version 1.50, https://yihui.org/knitr/.

<a id="ggplot2"></a>
[12] Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org.

<a id="enrichr"></a>
[13] Jawaid W (2025). enrichR: Provides an R Interface to 'Enrichr'. R package version 3.4, https://CRAN.R-project.org/package=enrichR.
