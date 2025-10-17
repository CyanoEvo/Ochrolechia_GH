library(tidyr)
library(stringr)
library(forcats)
library(taxonomizr)
library(stringr)
library(dplyr)
library(ggplot2)
library(cowplot)

theme_nc <- function(){
  theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 10), 
          axis.title.y = element_text(size = 10))
}
setwd("/home/nathan/Projects/CryptFunc_general/Ochrolechia_metagenomes/R")
#setwd("D:/Projects/CryptFunc_general/Ochrolechia_metagenomes/R")


###start_from here_with _new_file
#get a table just with unique pfams (no tax)
pfam_och_CPM.raw_no_tax <- pfam_och_CPM.raw
pfam_och_CPM.raw_no_tax <- pfam_och_CPM.raw_no_tax [!grepl("g__", pfam_och_CPM.raw_no_tax$function.),]
pfam_och_CPM.raw_no_tax <- pfam_och_CPM.raw_no_tax [!grepl("unclassified", pfam_och_CPM.raw_no_tax$function.),]

O_totals <- pfam_och_CPM.raw

row.names(O_totals) <- O_totals$function.
O_totals <- O_totals %>% rename(GH_PFAM = function.)
O_totals <-O_totals %>% filter(rowSums(across(where(is.numeric)))!=0)


GH_vector<-c(
  "PF00232", #GH1
  "PF00703", #GH2
  "PF00933", #GH3
  # "PF02506", #GH4
  "PF00150", #GH5
  #low"PF01341", #GH6
  "PF00840", #GH7
  #low"PF01270", #GH8, lichenase
  #low"PF00759", #GH9
  #low"PF00331", #GH10
  #low"PF00457", #GH11
  #low"PF01670", #GH12
  #"PF03701", #GH13
  #low"PF01373", #GH14
  "PF09137", #Glucodextranase (GH15)
  "PF00722", #GH16
  #low"PF00332", #GH17
  "PF00704", #GH18
  #  "PF00182", #Prokaryote chitinase (GH19)
  #low"PF00728", #GH20
  #  "PF00062", #GH22
  #"PF00959", #GH24
  #low"PF01183", #GH25
  #low"PF02156", #GH26
  "PF16499", #GH27
  #  "PF02065", #GH27 NOT SURE THIS IS THE RIGHT ONE TO USE
  "PF00295", #GH28
  #low"PF01120", #GH29
  #low"PF02055", #GH30
  "PF01055", #GH31
  "PF00251", #GH32
  #  "PF02973", #GH33
  # "PF00064", #GH34
  #low"PF01301", #GH35
  #low"PF17167", #GH36
  "PF01204", #Trehalase (GH37)
  #"PF18438", #GH38
  #low"PF01229", #GH39
  "PF04616", #GH43
  #low"PF02015", #GH45
  #"PF01374", #GH46
  "PF01532", #GH47
  # "PF03718", #GH49
  "PF06964", #GH51
  #low"PF03664", #GH62
  "PF03200", #GH63
  "PF16483", #GH64
  "PF03632", #GH65
  #  "PF13199", #GH66
  #low"PF03648", #GH67
  #"PF02435", #GH68
  #low"PF03659", #GH71
  "PF03198", #GH72
  #low"PF07335", #GH75
  "PF03663", #GH76
  #low"PF02446", #GH77
  # "PF21104", #GH78
  #"PF03662", #GH79
  #"PF13647", #GH80
  #low"PF03639", #GH81
  #"PF21809", #GH84
  #"PF03644", #GH85
  #low"PF07470", #GH88
  "PF07971", #GH92
  #"PF22124", #GH95
  #"PF10566", #GH97
  #"PF08306", #GH98
  #"PF16317", #GH99
  #"PF12899", #GH100
  #"PF12905", #GH101
  #low"PF17132", #GH106
  #"PF05838", #GH108
  #"PF21252", #GH109
  #low"PF03537", #GH114
  #low"PF15979", #GH115#
  # "PF21258", #GH120
  #"PF22680", #GH123
  #low"PF06824", #GH125
  "PF11790" #GH128
  #"PF11308", #GH129
  #low"PF04041", #GH130
  #"PF21087" #GH134
  # "PF01232" #mtdh rossmasn domain
)


#############
O_pfams_totals<-O_totals[GH_vector,]
#O_pfams_totals<-na.omit(O_pfams_totals)
O_pfams_totals_long <- gather(O_pfams_totals, sample, val, OFS_1:OFC_3, factor_key=TRUE)

O_pfams_totals_long$rep = O_pfams_totals_long$sample 
O_pfams_totals_long <- O_pfams_totals_long %>% mutate(across('rep', str_replace, 'OFS_1', 'R1'))
O_pfams_totals_long <- O_pfams_totals_long %>% mutate(across('rep', str_replace, 'OFS_2', 'R2'))
O_pfams_totals_long <- O_pfams_totals_long %>% mutate(across('rep', str_replace, 'OFS_3', 'R3'))
O_pfams_totals_long <- O_pfams_totals_long %>% mutate(across('rep', str_replace, 'OFC_1', 'R1'))
O_pfams_totals_long <- O_pfams_totals_long %>% mutate(across('rep', str_replace, 'OFC_2', 'R2'))
O_pfams_totals_long <- O_pfams_totals_long %>% mutate(across('rep', str_replace, 'OFC_3', 'R3'))

O_pfams_totals_long$growth_form = O_pfams_totals_long$sample 
O_pfams_totals_long <- O_pfams_totals_long %>% mutate(across('growth_form', str_replace, 'OFS_1', '1_Spines'))
O_pfams_totals_long <- O_pfams_totals_long %>% mutate(across('growth_form', str_replace, 'OFS_2', '1_Spines'))
O_pfams_totals_long <- O_pfams_totals_long %>% mutate(across('growth_form', str_replace, 'OFS_3', '1_Spines'))
O_pfams_totals_long <- O_pfams_totals_long %>% mutate(across('growth_form', str_replace, 'OFC_1', '2_Crust'))
O_pfams_totals_long <- O_pfams_totals_long %>% mutate(across('growth_form', str_replace, 'OFC_2', '2_Crust'))
O_pfams_totals_long <- O_pfams_totals_long %>% mutate(across('growth_form', str_replace, 'OFC_3', '2_Crust'))

O_pfams_totals_long$val_log<-log(O_pfams_totals_long$val+1)
O_pfams_totals_long.summary <- O_pfams_totals_long %>%
  group_by(growth_form,GH_PFAM) %>%
  summarise(
    sd = sd(val, na.rm = FALSE),
    val = mean(val)
  )
O_pfams_totals_long.summary

O_pfams_totals_long.summary_abundant <- subset(O_pfams_totals_long.summary, val >= 1)

kruskal.test(val ~ growth_form,
             data = O_pfams_totals_long.summary_abundant)

cazymes_bar <- ggplot(O_pfams_totals_long.summary, aes(x = GH_PFAM, y = val, ymin = val-sd, ymax = val+sd))
cazymes_plot<-cazymes_bar + 
  geom_errorbar(
    aes(ymin = val, ymax = val+sd, group = growth_form),
    width = 0.2, position = position_dodge(0.8))+
  geom_col(aes(fill = growth_form), position = position_dodge(0.8), width = 0.8)+
  scale_x_discrete(labels=c("PF00232"="GH1",
                            "PF00703"="GH2",
                            "PF00933"="GH3",
                            "PF02506"="GH4",
                            "PF00150"="GH5",
                            #low        "PF01341"="GH6",
                            "PF00840"="GH7",
                            #low        "PF01270"="GH8",
                            #low        "PF00759"="GH9",
                            #low        "PF00331"="GH10",
                            #low        "PF00457"="GH11",
                            #low        "PF01670"="GH12",
                            #low        "PF03701"="GH13",
                            #low        "PF01373"="GH14",
                            "PF09137"="Glucodextranase (GH15)",
                            "PF00722"="GH16",
                            #low    "PF00332"="GH17",
                            "PF00704"="GH18",
                            #   "PF00182"="Prokaryote chitinase (GH19)",
                            #low  "PF00728"="GH20",
                            #low    "PF00062"="GH22",
                            #low    "PF00959"="GH24",
                            #low    "PF01183"="GH25",
                            #low    "PF02156"="GH26",
                            "PF16499"="GH27",
                            #,"PF02065"="GH27",
                            "PF00295"="GH28",
                            #low    "PF01120"="GH29",
                            #low    "PF02055"="GH30",
                            "PF01055"="GH31",
                            "PF00251"="GH32",
                            #   "PF02973"="GH33",
                            #       "PF00064"="GH34",
                            #low    "PF01301"="GH35",
                            #low    "PF17167"="GH36",
                            "PF01204"="Trehalase (GH37)",
                            #   "PF18438"="GH38",
                            #low    "PF01229"="GH39",
                            "PF04616"="GH43",
                            #low    "PF02015"="GH45",
                            "PF01374"="GH46",
                            "PF01532"="GH47",
                            "PF03718"="GH49",
                            "PF06964"="GH51",
                            #low    "PF03664"="GH62",
                            "PF03200"="GH63",
                            "PF16483"="GH64",
                            "PF03632"="GH65",
                            "PF13199"="GH66",
                            #low    "PF03648"="GH67",
                            "PF02435"="GH68",
                            #low    "PF03659"="GH71",
                            "PF03198"="GH72",
                            #low    "PF07335"="GH75",
                            "PF03663"="GH76",
                            #low    "PF02446"="GH77",
                            #low    "PF21104"="GH78",
                            #low    "PF03662"="GH79",
                            #low    "PF13647"="GH80",
                            #low    "PF03639"="GH81",
                            #low    "PF21809"="GH84",
                            #low    "PF03644"="GH85",
                            #low    "PF07470"="GH88",
                            "PF07971"="GH92",
                            "PF22124"="GH95",
                            "PF10566"="GH97",
                            "PF08306"="GH98",
                            "PF16317"="GH99",
                            "PF12899"="GH100",
                            "PF12905"="GH101",
                            #low    "PF17132"="GH106",
                            #low    "PF05838"="GH108",
                            #low    "PF21252"="GH109",
                            #low    "PF03537"="GH114",
                            #low    "PF15979"="GH115",
                            #low    "PF21258"="GH120",
                            #low    "PF22680"="GH123",
                            #low    "PF06824"="GH125",
                            "PF11790"="GH128",
                            "PF11308"="GH129",
                            #low    "PF04041"="GH130",
                            "PF21087"="GH134"
                            #"PF01232"="mtdh"
  ),
  limits=rev(GH_vector)) +
  scale_fill_manual(values = c("firebrick", "lightblue"))+ 
  coord_flip() +  theme_nc()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        aspect.ratio = 4) 
+ geom_hline(yintercept=, linetype="dashed", 
             color = "red", linewidth=1)



####tax analysis


#make a column with the tax info
O_totals_1<-O_totals %>% separate(GH_PFAM, into = c("GH_PFAM", "TAX"), sep = "\\s*\\|\\s*")
O_totals_1<-O_totals_1 %>%  filter(!is.na(TAX))

O_totals_1$TAX <- gsub('g__', '', O_totals_1$TAX)


setwd("~/DATA/taxonomizr/")
O_totals_1$taxaId<-getId(O_totals_1$TAX,'accessionTaxa.sql')
taxa_1<-getTaxonomy(O_totals_1$taxaId,'accessionTaxa.sql')
setwd("~/Projects/CryptFunc_general/Ochrolechia_metagenomes/R")
df_1<-bind_cols(O_totals_1,taxa_1)
df_1_clean <- df_1
df_1_clean <- df_1_clean%>% mutate(phylum=replace_na(phylum, "Unknown"))
df_1_clean <- df_1_clean %>% mutate(class=replace_na(class, "Unknown"))

df_1_clean_long <- gather(df_1_clean, sample, val, OFS_1:OFC_3, factor_key=TRUE)

df_1_clean_long$rep = df_1_clean_long$sample 
df_1_clean_long <- df_1_clean_long %>% mutate(across('rep', str_replace, 'OFS_1', 'R1'))
df_1_clean_long <- df_1_clean_long %>% mutate(across('rep', str_replace, 'OFS_2', 'R2'))
df_1_clean_long <- df_1_clean_long %>% mutate(across('rep', str_replace, 'OFS_3', 'R3'))
df_1_clean_long <- df_1_clean_long %>% mutate(across('rep', str_replace, 'OFC_1', 'R1'))
df_1_clean_long <- df_1_clean_long %>% mutate(across('rep', str_replace, 'OFC_2', 'R2'))
df_1_clean_long <- df_1_clean_long %>% mutate(across('rep', str_replace, 'OFC_3', 'R3'))

df_1_clean_long$growth_form = df_1_clean_long$sample 
df_1_clean_long <- df_1_clean_long %>% mutate(across('growth_form', str_replace, 'OFS_1', 'Spines'))
df_1_clean_long <- df_1_clean_long %>% mutate(across('growth_form', str_replace, 'OFS_2', 'Spines'))
df_1_clean_long <- df_1_clean_long %>% mutate(across('growth_form', str_replace, 'OFS_3', 'Spines'))
df_1_clean_long <- df_1_clean_long %>% mutate(across('growth_form', str_replace, 'OFC_1', 'Crust'))
df_1_clean_long <- df_1_clean_long %>% mutate(across('growth_form', str_replace, 'OFC_2', 'Crust'))
df_1_clean_long <- df_1_clean_long %>% mutate(across('growth_form', str_replace, 'OFC_3', 'Crust'))



colors_vec <- c(Actinomycetes = "#c51b7d", Dothideomycetes = "#e9a3c9", Eurotiomycetes = "#fde0ef", Lecanoromycetes = "#f7f7f7", 
                Leotiomycetes = "#e6f5d0", Sordariomycetes = "#a1d76a", Xylonomycetes = "#4d9221", Alphaproteobacteria = "#a6611a",
                Betaproteobacteria = "#dfc27d", Ktedonobacteria = "#80cdc1", Terriglobia = "#018571", Unknown = "grey")

#GH1
GH1_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF00232")
GH1_spines_df<-GH1_df %>% filter(growth_form == "Spines")
GH1_spines_df.summary <- GH1_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH1_spines_df.summary  <- GH1_spines_df.summary  %>% mutate(percentage = (val / sum(val)) * 100)
GH1_spines_df.summary  <- GH1_spines_df.summary  %>% mutate(class = fct_relevel(class, 
                                                                                "Actinomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes", 
                                                                                "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Alphaproteobacteria",
                                                                                "Betaproteobacteria", "Ktedonobacteria", "Terriglobia", "Unknown"))
GH1_spines_plot<-ggplot(GH1_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() +theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec) 

GH1_crust_df<-GH1_df %>% filter(growth_form == "Crust")
GH1_crust_df.summary <- GH1_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH1_crust_df.summary  <- GH1_crust_df.summary  %>% mutate(percentage = (val / sum(val)) * 100)
GH1_crust_df.summary  <- GH1_crust_df.summary  %>% mutate(class = fct_relevel(class, 
                                                                              "Actinomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes", 
                                                                              "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Alphaproteobacteria",
                                                                              "Betaproteobacteria", "Ktedonobacteria", "Terriglobia", "Unknown"))
GH1_crust_plot<-ggplot(GH1_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() +theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

#GH2
GH2_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF00703")
GH2_spines_df<-GH2_df %>% filter(growth_form == "Spines")
GH2_spines_df.summary <- GH2_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH2_spines_df.summary  <- GH2_spines_df.summary  %>% mutate(percentage = (val / sum(val)) * 100)
GH2_spines_df.summary  <- GH2_spines_df.summary  %>% mutate(class = fct_relevel(class, 
                                                                                "Actinomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes", 
                                                                                "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Alphaproteobacteria",
                                                                                "Betaproteobacteria", "Ktedonobacteria", "Terriglobia", "Unknown"))
GH2_spines_plot<-ggplot(GH2_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH2_crust_df<-GH2_df %>% filter(growth_form == "Crust")
GH2_crust_df.summary <- GH2_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH2_crust_df.summary  <- GH2_crust_df.summary  %>% mutate(percentage = (val / sum(val)) * 100)
GH2_crust_df.summary  <- GH2_crust_df.summary  %>% mutate(class = fct_relevel(class, 
                                                                              "Actinomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes", 
                                                                              "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Alphaproteobacteria",
                                                                              "Betaproteobacteria", "Ktedonobacteria", "Terriglobia", "Unknown"))
GH2_crust_plot<-ggplot(GH2_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

#GH3
GH3_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF00933")
GH3_spines_df<-GH3_df %>% filter(growth_form == "Spines")
GH3_spines_df.summary <- GH3_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH3_spines_df.summary  <- GH3_spines_df.summary  %>% mutate(percentage = (val / sum(val)) * 100)
GH3_spines_df.summary  <- GH3_spines_df.summary  %>% mutate(class = fct_relevel(class, 
                                                                                "Actinomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes", 
                                                                                "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Alphaproteobacteria",
                                                                                "Betaproteobacteria", "Ktedonobacteria", "Terriglobia", "Unknown"))
GH3_spines_plot<-ggplot(GH3_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH3_crust_df<-GH3_df %>% filter(growth_form == "Crust")
GH3_crust_df.summary <- GH3_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH3_crust_df.summary  <- GH3_crust_df.summary  %>% mutate(percentage = (val / sum(val)) * 100)
GH3_crust_df.summary  <- GH3_crust_df.summary  %>% mutate(class = fct_relevel(class, 
                                                                              "Actinomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes", 
                                                                              "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Alphaproteobacteria",
                                                                              "Betaproteobacteria", "Ktedonobacteria", "Terriglobia", "Unknown"))
GH3_crust_plot<-ggplot(GH3_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

#GH5
GH5_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF00150")
GH5_spines_df<-GH5_df %>% filter(growth_form == "Spines")
GH5_spines_df.summary <- GH5_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH5_spines_df.summary  <- GH5_spines_df.summary  %>% mutate(percentage = (val / sum(val)) * 100)
GH5_spines_df.summary  <- GH5_spines_df.summary  %>% mutate(class = fct_relevel(class, 
                                                                                "Actinomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes", 
                                                                                "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Alphaproteobacteria",
                                                                                "Betaproteobacteria", "Ktedonobacteria", "Terriglobia", "Unknown"))
GH5_spines_plot<-ggplot(GH5_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH5_crust_df<-GH5_df %>% filter(growth_form == "Crust")
GH5_crust_df.summary <- GH5_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH5_crust_df.summary  <- GH5_crust_df.summary  %>% mutate(percentage = (val / sum(val)) * 100)
GH5_crust_df.summary  <- GH5_crust_df.summary  %>% mutate(class = fct_relevel(class, 
                                                                              "Actinomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes", 
                                                                              "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Alphaproteobacteria",
                                                                              "Betaproteobacteria", "Ktedonobacteria", "Terriglobia", "Unknown"))
GH5_crust_plot<-ggplot(GH5_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

#GH7
GH7_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF00840")
GH7_spines_df<-GH7_df %>% filter(growth_form == "Spines")
GH7_spines_df.summary <- GH7_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH7_spines_df.summary  <- GH7_spines_df.summary  %>% mutate(class = fct_relevel(class, 
                                                                                "Actinomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes", 
                                                                                "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Alphaproteobacteria",
                                                                                "Betaproteobacteria", "Ktedonobacteria", "Terriglobia", "Unknown"))
GH7_spines_plot<-ggplot(GH7_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH7_crust_df<-GH7_df %>% filter(growth_form == "Crust")
GH7_crust_df.summary <- GH7_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH7_crust_df.summary  <- GH7_crust_df.summary  %>% mutate(class = fct_relevel(class, 
                                                                              "Actinomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes", 
                                                                              "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Alphaproteobacteria",
                                                                              "Betaproteobacteria", "Ktedonobacteria", "Terriglobia", "Unknown"))
GH7_crust_plot<-ggplot(GH7_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

#GH15
GH15_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF09137")
GH15_spines_df<-GH15_df %>% filter(growth_form == "Spines")
GH15_spines_df.summary <- GH15_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH15_spines_df.summary  <- GH15_spines_df.summary  %>% mutate(class = fct_relevel(class, 
                                                                                  "Actinomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes", 
                                                                                  "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Alphaproteobacteria",
                                                                                  "Betaproteobacteria", "Ktedonobacteria", "Terriglobia", "Unknown"))
GH15_spines_plot<-ggplot(GH15_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH15_crust_df<-GH15_df %>% filter(growth_form == "Crust")
GH15_crust_df.summary <- GH15_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH15_crust_df.summary  <- GH15_crust_df.summary  %>% mutate(class = fct_relevel(class, 
                                                                                "Actinomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes", 
                                                                                "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Alphaproteobacteria",
                                                                                "Betaproteobacteria", "Ktedonobacteria", "Terriglobia", "Unknown"))
GH15_crust_plot<-ggplot(GH15_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

#GH16
GH16_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF00722")
GH16_spines_df<-GH16_df %>% filter(growth_form == "Spines")
GH16_spines_df.summary <- GH16_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH16_spines_df.summary  <- GH16_spines_df.summary  %>% mutate(percentage = (val / sum(val)) * 100)
GH16_spines_df.summary  <- GH16_spines_df.summary  %>% mutate(class = fct_relevel(class, 
                                                                                  "Actinomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes", 
                                                                                  "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Alphaproteobacteria",
                                                                                  "Betaproteobacteria", "Ktedonobacteria", "Terriglobia", "Unknown"))
GH16_spines_plot<-ggplot(GH16_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH16_crust_df<-GH16_df %>% filter(growth_form == "Crust")
GH16_crust_df.summary <- GH16_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH16_crust_df.summary  <- GH16_crust_df.summary  %>% mutate(percentage = (val / sum(val)) * 100)
GH16_crust_df.summary  <- GH16_crust_df.summary  %>% mutate(class = fct_relevel(class, 
                                                                                "Actinomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes", 
                                                                                "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Alphaproteobacteria",
                                                                                "Betaproteobacteria", "Ktedonobacteria", "Terriglobia", "Unknown"))
GH16_crust_plot<-ggplot(GH16_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

#GH18
GH18_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF00704")
GH18_spines_df<-GH18_df %>% filter(growth_form == "Spines")
GH18_spines_df.summary <- GH18_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH18_spines_df.summary  <- GH18_spines_df.summary  %>% mutate(class = fct_relevel(class, 
                                                                                  "Actinomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes", 
                                                                                  "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Alphaproteobacteria",
                                                                                  "Betaproteobacteria", "Ktedonobacteria", "Terriglobia", "Unknown"))
GH18_spines_plot<-ggplot(GH18_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH18_crust_df<-GH18_df %>% filter(growth_form == "Crust")
GH18_crust_df.summary <- GH18_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH18_crust_df.summary  <- GH18_crust_df.summary  %>% mutate(class = fct_relevel(class, 
                                                                                "Actinomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes", 
                                                                                "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Alphaproteobacteria",
                                                                                "Betaproteobacteria", "Ktedonobacteria", "Terriglobia", "Unknown"))
GH18_crust_plot<-ggplot(GH18_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)


#GH27
GH27_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF16499")
GH27_spines_df<-GH27_df %>% filter(growth_form == "Spines")
GH27_spines_df.summary <- GH27_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH27_spines_df.summary  <- GH27_spines_df.summary  %>% mutate(class = fct_relevel(class, 
                                                                                  "Actinomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes", 
                                                                                  "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Alphaproteobacteria",
                                                                                  "Betaproteobacteria", "Ktedonobacteria", "Terriglobia", "Unknown"))
GH27_spines_plot<-ggplot(GH27_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH27_crust_df<-GH27_df %>% filter(growth_form == "Crust")
GH27_crust_df.summary <- GH27_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH27_crust_df.summary  <- GH27_crust_df.summary  %>% mutate(class = fct_relevel(class, 
                                                                                "Actinomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes", 
                                                                                "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Alphaproteobacteria",
                                                                                "Betaproteobacteria", "Ktedonobacteria", "Terriglobia", "Unknown"))
GH27_crust_plot<-ggplot(GH27_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

#GH28
GH28_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF00295")
GH28_spines_df<-GH28_df %>% filter(growth_form == "Spines")
GH28_spines_df.summary <- GH28_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH28_spines_plot<-ggplot(GH28_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH28_crust_df<-GH28_df %>% filter(growth_form == "Crust")
GH28_crust_df.summary <- GH28_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH28_crust_plot<-ggplot(GH28_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)
#GH31
GH31_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF01055")
GH31_spines_df<-GH31_df %>% filter(growth_form == "Spines")
GH31_spines_df.summary <- GH31_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH31_spines_df.summary  <- GH31_spines_df.summary  %>% mutate(percentage = (val / sum(val)) * 100)
GH31_spines_plot<-ggplot(GH31_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH31_crust_df<-GH31_df %>% filter(growth_form == "Crust")
GH31_crust_df.summary <- GH31_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH31_crust_df.summary  <- GH31_crust_df.summary  %>% mutate(percentage = (val / sum(val)) * 100)
GH31_crust_plot<-ggplot(GH31_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

#GH32
GH32_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF00251")
GH32_spines_df<-GH32_df %>% filter(growth_form == "Spines")
GH32_spines_df.summary <- GH32_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH32_spines_plot<-ggplot(GH32_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH32_crust_df<-GH32_df %>% filter(growth_form == "Crust")
GH32_crust_df.summary <- GH32_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH32_crust_plot<-ggplot(GH32_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)


#GH37
GH37_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF01204")
GH37_spines_df<-GH37_df %>% filter(growth_form == "Spines")
GH37_spines_df.summary <- GH37_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH37_spines_df.summary  <- GH37_spines_df.summary  %>% mutate(percentage = (val / sum(val)) * 100)
GH37_spines_plot<-ggplot(GH37_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH37_crust_df<-GH37_df %>% filter(growth_form == "Crust")
GH37_crust_df.summary <- GH37_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH37_crust_df.summary  <- GH37_crust_df.summary  %>% mutate(percentage = (val / sum(val)) * 100)
GH37_crust_plot<-ggplot(GH37_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)


#GH43
GH43_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF04616")
GH43_spines_df<-GH43_df %>% filter(growth_form == "Spines")
GH43_spines_df.summary <- GH43_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH43_spines_plot<-ggplot(GH43_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH43_crust_df<-GH43_df %>% filter(growth_form == "Crust")
GH43_crust_df.summary <- GH43_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH43_crust_plot<-ggplot(GH43_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

#GH47
GH47_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF01532")
GH47_spines_df<-GH47_df %>% filter(growth_form == "Spines")
GH47_spines_df.summary <- GH47_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH47_spines_plot<-ggplot(GH47_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH47_crust_df<-GH47_df %>% filter(growth_form == "Crust")
GH47_crust_df.summary <- GH47_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH47_crust_plot<-ggplot(GH47_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)


#GH51
GH51_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF06964")
GH51_spines_df<-GH51_df %>% filter(growth_form == "Spines")
GH51_spines_df.summary <- GH51_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH51_spines_plot<-ggplot(GH51_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH51_crust_df<-GH51_df %>% filter(growth_form == "Crust")
GH51_crust_df.summary <- GH51_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH51_crust_plot<-ggplot(GH51_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

#GH63
GH63_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF03200")
GH63_spines_df<-GH63_df %>% filter(growth_form == "Spines")
GH63_spines_df.summary <- GH63_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH63_spines_df.summary  <- GH63_spines_df.summary  %>% mutate(percentage = (val / sum(val)) * 100)
GH63_spines_plot<-ggplot(GH63_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH63_crust_df<-GH63_df %>% filter(growth_form == "Crust")
GH63_crust_df.summary <- GH63_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH63_crust_df.summary  <- GH63_crust_df.summary  %>% mutate(percentage = (val / sum(val)) * 100)
GH63_crust_plot<-ggplot(GH63_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

#GH64
GH64_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF16483")
GH64_spines_df<-GH64_df %>% filter(growth_form == "Spines")
GH64_spines_df.summary <- GH64_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH64_spines_df.summary  <- GH64_spines_df.summary  %>% mutate(percentage = (val / sum(val)) * 100)
GH64_spines_plot<-ggplot(GH64_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH64_crust_df<-GH64_df %>% filter(growth_form == "Crust")
GH64_crust_df.summary <- GH64_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH64_crust_df.summary  <- GH64_crust_df.summary  %>% mutate(percentage = (val / sum(val)) * 100)
GH64_crust_plot<-ggplot(GH64_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

#GH65
GH65_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF03632")
GH65_spines_df<-GH65_df %>% filter(growth_form == "Spines")
GH65_spines_df.summary <- GH65_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH65_spines_plot<-ggplot(GH65_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH65_crust_df<-GH65_df %>% filter(growth_form == "Crust")
GH65_crust_df.summary <- GH65_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH65_crust_plot<-ggplot(GH65_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

#GH72
GH72_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF03198")
GH72_spines_df<-GH72_df %>% filter(growth_form == "Spines")
GH72_spines_df.summary <- GH72_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH72_spines_plot<-ggplot(GH72_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH72_crust_df<-GH72_df %>% filter(growth_form == "Crust")
GH72_crust_df.summary <- GH72_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH72_crust_plot<-ggplot(GH72_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

#GH76
GH76_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF03663")
GH76_spines_df<-GH76_df %>% filter(growth_form == "Spines")
GH76_spines_df.summary <- GH76_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH76_spines_plot<-ggplot(GH76_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH76_crust_df<-GH76_df %>% filter(growth_form == "Crust")
GH76_crust_df.summary <- GH76_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH76_crust_plot<-ggplot(GH76_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

#GH92
GH92_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF07971")
GH92_spines_df<-GH92_df %>% filter(growth_form == "Spines")
GH92_spines_df.summary <- GH92_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH92_spines_plot<-ggplot(GH92_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH92_crust_df<-GH92_df %>% filter(growth_form == "Crust")
GH92_crust_df.summary <- GH92_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH92_crust_plot<-ggplot(GH92_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

#GH128
GH128_df<-df_1_clean_long  %>% filter(GH_PFAM == "PF11790")
GH128_spines_df<-GH128_df %>% filter(growth_form == "Spines")
GH128_spines_df.summary <- GH128_spines_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH128_spines_plot<-ggplot(GH128_spines_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

GH128_crust_df<-GH128_df %>% filter(growth_form == "Crust")
GH128_crust_df.summary <- GH128_crust_df %>% group_by(class) %>% summarise(sd = sd(val, na.rm = TRUE),val = mean(val))
GH128_crust_plot<-ggplot(GH128_crust_df.summary,aes(x="", y=val, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_flip() + scale_y_reverse() + theme_void()+ theme(aspect.ratio = 0.1)+ scale_fill_manual(values = colors_vec)

abundance_grid <-plot_grid(GH1_spines_plot,GH1_crust_plot,
                           GH2_spines_plot,GH2_crust_plot,
                           GH3_spines_plot,GH3_crust_plot,
                           GH5_spines_plot,GH5_crust_plot,
                           GH7_spines_plot,GH7_crust_plot,
                           GH15_spines_plot,GH15_crust_plot,
                           GH16_spines_plot,GH16_crust_plot,
                           GH18_spines_plot,GH18_crust_plot,
                           GH27_spines_plot,GH27_crust_plot,
                           GH28_spines_plot,GH28_crust_plot,
                           GH31_spines_plot,GH31_crust_plot,
                           GH32_spines_plot,GH32_crust_plot,
                           GH37_spines_plot,GH37_crust_plot,
                           GH43_spines_plot,GH43_crust_plot,
                           GH47_spines_plot,GH47_crust_plot,
                           GH51_spines_plot,GH51_crust_plot,
                           GH63_spines_plot,GH63_crust_plot,
                           GH64_spines_plot,GH64_crust_plot,
                           GH65_spines_plot,GH65_crust_plot,
                           GH72_spines_plot,GH72_crust_plot,
                           GH76_spines_plot,GH76_crust_plot,
                           GH92_spines_plot,GH92_crust_plot,
                           GH128_spines_plot,GH128_crust_plot,
                           labels = c("GH1:Spines","GH1:Crust",
                                      "GH2:Spines","GH2:Crust",
                                      "GH3:Spines","GH3:Crust",
                                      "GH5:Spines","GH5:Crust",
                                      "GH7:Spines","GH7:Crust",
                                      "GH15:Spines","GH15:Crust",
                                      "GH16:Spines","GH16:Crust",
                                      "GH18:Spines","GH18:Crust",
                                      "GH27:Spines","GH27:Crust",
                                      "GH28:Spines","GH28:Crust",
                                      "GH31:Spines","GH31:Crust",
                                      "GH32:Spines","GH32:Crust",
                                      "GH37:Spines","GH37:Crust",
                                      "GH43:Spines","GH43:Crust",
                                      "GH47:Spines","GH47:Crust",
                                      "GH51:Spines","GH51:Crust",
                                      "GH63:Spines","GH63:Crust",
                                      "GH64:Spines","GH64:Crust",
                                      "GH65:Spines","GH65:Crust",
                                      "GH72:Spines","GH72:Crust",
                                      "GH76:Spines","GH76:Crust",
                                      "GH92:Spines","GH92:Crust",
                                      "GH128:Spines","GH128:Crust"),
                           label_x = 0,
                           ncol = 2, align = "hv")

plot_grid(cazymes_plot,abundance_grid)
