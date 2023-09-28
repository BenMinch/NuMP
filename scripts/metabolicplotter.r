#Metabolic Plotter
library(tidyverse)

#import data as first argument
args<- commandArgs(trailingOnly=TRUE)
#read in the data
data<- read.csv(args[1])
output_dir<- args[2]
mode<- args[3]
taxonomy<- read.csv(args[4],sep="\t")

dictionary<- read.csv('resources/metabolic_dictionary.tsv',sep="\t")




#split query_id at . and take first part as genome column
data$genome<- strsplit(as.character(data$query_id),split="\\.") %>% sapply(`[[`,1)

#make a list of dictionary genes
dictionary_genes<- dictionary$Gene %>% unique

#filter data so that you only keep rows where id is in the dictionary pfam column
data<- data %>% filter(id %in% dictionary$pfam_name)%>% select(query_id,id,genome)


data$Pathway<- dictionary$Pathway[match(data$id,dictionary$pfam_name)]
data$Gene<- dictionary$Gene[match(data$id,dictionary$pfam_name)]

#split at 2nd to last _ and take all before that as the genome

data$family<- taxonomy$family[match(data$genome,taxonomy$genome)]

#Check the mode
if(mode=="family"){
  #group by family and Gene and get the proportion of members in that family that have the gene, but only count each genome once

  data_family<- data %>% group_by(family,Gene) %>% summarise(Genomes_with_gene=n_distinct(genome))
  data_family2<- data %>% group_by(family) %>% summarise(Genomes_total=n_distinct(genome))
  merged_family<-left_join(data_family,data_family2,by="family")%>% mutate(Proportion=Genomes_with_gene/Genomes_total)
  merged_family$Pathway<- dictionary$Pathway[match(merged_family$Gene,dictionary$Gene)]
  #order by pathway with NCLDV_core at top
  merged_family$Pathway<- factor(merged_family$Pathway,levels=c("NCLDV_core","Glycolysis/Gluconeogenesis","TCA_cycle","Pentose phosphate pathway","Light Harvesting","Beta oxidation","Cytoskeleton","Transporter","DNA processing","Nutrient metabolism","Oxidative stress regulation"))
  #make a bubble plot with family on the x axis and all the genes on the y axis, colored by pathway

    bubble <- ggplot(merged_family, aes(x = family, y = Gene, color = Pathway)) +
    geom_point(aes(size = Proportion)) +
    scale_size_continuous(limits = c(0, 1), range = c(1, 10), breaks = c(0, 0.1, 0.25, 0.5, 0.75, 1)) +
    scale_y_discrete(limits = rev(dictionary_genes)) +
    theme(legend.key = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold"),  # Rotate and make bold
            axis.text.y = element_text(face = "bold"),  # Make y-axis text bold
            axis.title = element_text(face = "bold"),
            panel.border=element_rect(color = 'black',fill=NA,size=1.2),
            ) +
        ggtitle("Family Level Metabolism") +
    scale_color_manual(values = c('NCLDV_core' = '#E20001', 'Glycolysis/Gluconeogenesis' = '#A5E25A', 'TCA_cycle' = '#AA4400', 'Pentose phosphate pathway' = '#CD67E6', 'Light Harvesting' = '#008C00', 'Beta oxidation' = '#000000', 'Cytoskeleton' = '#FD7F81', 'Transporter' = '#ADAD00', 'DNA processing' = '#047AC5', 'Nutrient metabolism' = '#F9B88B', 'Oxidative stress regulation' = '#00FFFE'))
   ggsave(paste0(output_dir,"/bubble_family.pdf"),width=10,height=10)
}

if(mode=='all'){
    data_all<- data%>% group_by(genome,Gene)%>% summarise(Count=n())
    data_all$Pathway<- dictionary$Pathway[match(data_all$Gene,dictionary$Gene)]
    data_all$Pathway<- factor(data_all$Pathway,levels=c("NCLDV_core","Glycolysis/Gluconeogenesis","TCA_cycle","Pentose phosphate pathway","Light Harvesting","Beta oxidation","Cytoskeleton","Transporter","DNA processing","Nutrient metabolism","Oxidative stress regulation"))

    bubble_all <- ggplot(data_all, aes(x = genome, y = Gene, color = Pathway)) +
    geom_point(aes(size = Count)) +
    scale_size_continuous(limits = c(0, max(data_all$Count)), range = c(1, 10), breaks = c(1, 2, 5, 10)) +
    scale_y_discrete(limits = rev(dictionary_genes)) +
    theme(legend.key = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold"),  # Rotate and make bold
            axis.text.y = element_text(face = "bold"),  # Make y-axis text bold
            axis.title = element_text(face = "bold"),
            panel.border=element_rect(color = 'black',fill=NA,size=1.2),
            ) +
        ggtitle("Family Level Metabolism") +
    scale_color_manual(values = c('NCLDV_core' = '#E20001', 'Glycolysis/Gluconeogenesis' = '#A5E25A', 'TCA_cycle' = '#AA4400', 'Pentose phosphate pathway' = '#CD67E6', 'Light Harvesting' = '#008C00', 'Beta oxidation' = '#000000', 'Cytoskeleton' = '#FD7F81', 'Transporter' = '#ADAD00', 'DNA processing' = '#047AC5', 'Nutrient metabolism' = '#F9B88B', 'Oxidative stress regulation' = '#00FFFE'))
   ggsave(paste0(output_dir,"/bubble_genomes.pdf"),width=15,height=10)
}

if(mode=='both'){
    data_all<- data%>% group_by(genome,Gene)%>% summarise(Count=n())
    data_all$Pathway<- dictionary$Pathway[match(data_all$Gene,dictionary$Gene)]
    data_all$Pathway<- factor(data_all$Pathway,levels=c("NCLDV_core","Glycolysis/Gluconeogenesis","TCA_cycle","Pentose phosphate pathway","Light Harvesting","Beta oxidation","Cytoskeleton","Transporter","DNA processing","Nutrient metabolism","Oxidative stress regulation"))

    bubble_all <- ggplot(data_all, aes(x = genome, y = Gene, color = Pathway)) +
    geom_point(aes(size = Count)) +
    scale_size_continuous(limits = c(0, max(data_all$Count)), range = c(1, 10), breaks = c(1, 2, 5, 10)) +
    scale_y_discrete(limits = rev(dictionary_genes)) +
    theme(legend.key = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold"),  # Rotate and make bold
            axis.text.y = element_text(face = "bold"),  # Make y-axis text bold
            axis.title = element_text(face = "bold"),
            panel.border=element_rect(color = 'black',fill=NA,size=1.2),
            ) +
        ggtitle("Family Level Metabolism") +
    scale_color_manual(values = c('NCLDV_core' = '#E20001', 'Glycolysis/Gluconeogenesis' = '#A5E25A', 'TCA_cycle' = '#AA4400', 'Pentose phosphate pathway' = '#CD67E6', 'Light Harvesting' = '#008C00', 'Beta oxidation' = '#000000', 'Cytoskeleton' = '#FD7F81', 'Transporter' = '#ADAD00', 'DNA processing' = '#047AC5', 'Nutrient metabolism' = '#F9B88B', 'Oxidative stress regulation' = '#00FFFE'))
   ggsave(paste0(output_dir,"/bubble_genomes.pdf"),width=15,height=10)

    data_family<- data %>% group_by(family,Gene) %>% summarise(Genomes_with_gene=n_distinct(genome))
  data_family2<- data %>% group_by(family) %>% summarise(Genomes_total=n_distinct(genome))
  merged_family<-left_join(data_family,data_family2,by="family")%>% mutate(Proportion=Genomes_with_gene/Genomes_total)
  merged_family$Pathway<- dictionary$Pathway[match(merged_family$Gene,dictionary$Gene)]
  #order by pathway with NCLDV_core at top
  merged_family$Pathway<- factor(merged_family$Pathway,levels=c("NCLDV_core","Glycolysis/Gluconeogenesis","TCA_cycle","Pentose phosphate pathway","Light Harvesting","Beta oxidation","Cytoskeleton","Transporter","DNA processing","Nutrient metabolism","Oxidative stress regulation"))
  #make a bubble plot with family on the x axis and all the genes on the y axis, colored by pathway

    bubble <- ggplot(merged_family, aes(x = family, y = Gene, color = Pathway)) +
    geom_point(aes(size = Proportion)) +
    scale_size_continuous(limits = c(0, 1), range = c(1, 10), breaks = c(0, 0.1, 0.25, 0.5, 0.75, 1)) +
    scale_y_discrete(limits = rev(dictionary_genes)) +
    theme(legend.key = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold"),  # Rotate and make bold
            axis.text.y = element_text(face = "bold"),  # Make y-axis text bold
            axis.title = element_text(face = "bold"),
            panel.border=element_rect(color = 'black',fill=NA,size=1.2),
            ) +
        ggtitle("Family Level Metabolism") +
    scale_color_manual(values = c('NCLDV_core' = '#E20001', 'Glycolysis/Gluconeogenesis' = '#A5E25A', 'TCA_cycle' = '#AA4400', 'Pentose phosphate pathway' = '#CD67E6', 'Light Harvesting' = '#008C00', 'Beta oxidation' = '#000000', 'Cytoskeleton' = '#FD7F81', 'Transporter' = '#ADAD00', 'DNA processing' = '#047AC5', 'Nutrient metabolism' = '#F9B88B', 'Oxidative stress regulation' = '#00FFFE'))
   ggsave(paste0(output_dir,"/bubble_family.pdf"),width=10,height=10)
}

#write out csv

write.csv(data,paste0(output_dir,"/metabolic_data.csv"),row.names=FALSE)

