require(scales)
library(tidyverse)
library(cowplot)

build_summary <- function(betas,pip, pattern){
  betas %>% 
    inner_join(EWAS_Catalog_03_07_2019) %>%
    inner_join(pip) %>%
    dplyr::select(Author, PMID, CpG,Trait, Beta, Tissue, effect, PIP) %>%
    filter(str_detect(tolower(Trait), pattern = pattern))  %>%
    mutate(congruency = ifelse(sign(as.numeric(Beta)) == sign(effect),"Same","Different")) %>%
    mutate_at('Beta',as.numeric) %>%
    arrange(desc(PIP)) %>%
   mutate(Tissue = ifelse((str_detect(tolower(Tissue),"blood")), "Blood", Tissue))
    
}

plot_effects_cat <- function( trait_summary){
 trait_summary$CpG <- factor(trait_summary$CpG,levels = unique(trait_summary$CpG))
 
 trait.effects <- trait_summary %>% 
   mutate(Beta= abs(Beta))  %>% 
   mutate(effect=abs(effect)) %>% 
   filter(Beta<1) %>% 
   filter(effect > 0) %>%
   mutate(Tissue=ifelse((str_detect(tolower(Tissue),"blood")), "Blood", "Other")) %>%
   ggplot() + 
   geom_point(aes(x=CpG, y = effect), shape =4) + 
   theme_light() + 
   theme(axis.title.x = element_blank(),legend.position="top", axis.text.x=element_blank()) + 
   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
   geom_point(aes(x = CpG, y = Beta,colour = congruency, shape=Tissue), alpha=0.3, stroke=0.4) + scale_x_discrete(breaks = NULL) +
   annotation_logticks(sides = "l")  +
   labs(y= "Absolute effect-size (log10 scale)", colour = "Sign congruency")
 
 trait.pip <- trait_summary %>% 
   mutate(Beta= abs(Beta))  %>% 
   mutate(effect=abs(effect)) %>% 
   filter(Beta<1) %>% 
   filter(effect > 0) %>% 
   ggplot() + 
   geom_line(aes(x=CpG, y = PIP, group =1)) + 
   scale_x_discrete(breaks = NULL) + 
   theme_light() + 
   theme(axis.text.x = element_blank()) +
   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
   annotation_logticks(sides="l")  +
   labs(y= "PIP (log10 scale)")
 
 
  list(effects = trait.effects, pip = trait.pip)
}

