### Loading necessary packages and sourcing necessary functions:
library(maftools)
library(readxl)
library(stringr)
source('/Users/zuckerm/Documents/biallelic/analysis01/4-general_biallelic.R')
source('/Users/zuckerm/Documents/biallelic/analysis01/10-clinical.R')
source('/Users/zuckerm/Documents/biallelic/analysis01/8-plotting.R')

#Path where data is located:
dataPath <- '/Users/zuckerm/Documents/biallelic/data'
#path where figures are to be saved
plotPath <- '/Users/zuckerm/Documents/biallelic/analysis01/figures'

###
mutData.tcga <- get(load(paste(dataPath,'/mutationData.tcga.rda',sep = '')))
mutData.impact <- read.maf(paste(dataPath,'/msk_impact_facets_annotated.ccf.maf.gz',sep = ''))@data
clinData.tcga <- read.table(paste(dataPath,'/clinical_data_TCGA.txt',sep = ''),sep='\t',header=T)
clinData.impact <- readxl::read_excel(paste(dataPath,'/data_clinical_sample.oncokb.xlsx',sep = ''))
#TCGA RNA and protein data:
expressionData <- read.table(paste(dataPath,'/expressionDataSet.txt',sep=''),sep='\t',header=T)
proteinData <- read.table(paste(dataPath,'/proteinDataSet.txt',sep=''),sep='\t',header=T)
#Columnames have to be redone since hiphens are automatically turned into periods whenever saved:
colnames(proteinData) <- sapply(colnames(proteinData), function(x){paste(strsplit(x, split = '[.]')[[1]], collapse='-')})
colnames(expressionData) <- sapply(colnames(expressionData), function(x){paste(strsplit(x, split = '[.]')[[1]], collapse='-')})
#Tables of subtypes and TSGs included in analysis:
subtypeTable <- read.table(paste(dataPath,'/subtypeTable.txt',sep=''), sep='\t',header=T,comment.char = '&')
tsgdf <- read.table(paste(dataPath,'/custom_gene_list_autosomes.txt',sep=''),sep='\t')
tsgs <- tsgdf$Hugo_Symbol
#List of impact genes:
impact_genes <- get(load(paste(dataPath,'/impact_genes.rda',sep='')))

#Loading biallelic enrichment results:
biallelicEnrichmentRes <- get(load(paste(dataPath,'/biallelicEnrichment.rda',sep='')))
biallelicEnrichment.vus <- biallelicEnrichmentRes$loh.vus
biallelicEnrichment.driver <- biallelicEnrichmentRes$loh.oncogenic

#Zygosity calls:
filtered <- get(load(paste(dataPath,'/zygosityDataFiltered_current.rda',sep='')))
#Just drivers:
drivers <- get(load(paste(dataPath,'/zygosityDataDrivers_current.rda',sep='')))
filtered.tcga <- read.table(paste(dataPath,'/tcga_zygosity.txt',sep=''),sep='\t',header=T)
tcga_subtype_map <- read.table(paste(dataPath,'/tcga_subtype_map.txt',sep=''),sep='\t')

expressionResPath <- '/Users/zuckerm/Documents/biallelic/analysis01/results/expression/keap1_expression_results.rda'
reactomePathway <- "/Users/zuckerm/Documents/biallelic/data/c2.cp.reactome.v6.2.symbols.gmt.txt"
hallmarksPathway <- "/Users/zuckerm/Documents/biallelic/data/h.all.v2022.1.Hs.symbols.gmt.txt"

#Adding disease subtype column if necessary:
#if(is.null(mutData.impact$disease_subtype)){mutData.impact$disease_subtype <- 
#  clinData.impact$CANCER_TYPE_DETAILED[match(mutData.impact$Tumor_Sample_Barcode, clinData.impact$SAMPLE_ID)]}
#geneids <- get(load('/Users/zuckerm/Documents/Auxiliary/geneids.rda'))
#ens <- rownames(expressionData)
#ens.unversioned <- unname(sapply(ens, function(string){strsplit(string,split='[.]')[[1]][1]}))
#hugo <- sapply(ens.unversioned, function(x){geneids$gene[which(geneids$ensembl == x)]})

#Figure 4: VUSs:
#'stratify' refers to what gene to stratify survival plots by, e.g., EGFR mutation status. If NULL, it doesn't stratify by anything.
generateFigure4 <- function(drivers, filtered, biallelicEnrichment.driver, biallelicEnrichment.vus, mutData.tcga, mainsize=10, stratify=NULL){
  vusdat <- biallelicEnrichment.vus[,c(1:7)]
  colnames(vusdat)[c(3:7)] <- c('OR.vus','p.value.vus','p.value.correcrted.vus','mutloh.vus','mut.vus')
  enrdat <- merge(biallelicEnrichment.driver[,c(1:7)], vusdat, by=c('gene','disease'), )
  enrdat$Mutations <- enrdat$mut + enrdat$mutloh + enrdat$mut.vus + enrdat$mutloh.vus
  colnames(enrdat)[c(1,2,5,10)] <- c('Gene','Disease ','P_Value_Driver','P_Value_VUS')
  enrdat <- enrdat[which(enrdat$Disease != 'all'),]
  enrdat <- enrdat[which(enrdat$Gene %in% tsgs & enrdat$Mutations >= 20 & (enrdat$P_Value_Driver < .05 | enrdat$P_Value_VUS < .05)),]
  enrdat$label <- sapply(enrdat$Disease, function(d){
    str <- strsplit(d, split='[ ]')[[1]]
    l <- length(str)
    if(l == 1){out <- d}else{
      a <- floor(l/2)
      paste(paste(str[1:a], collapse=' '),'\n', paste(str[(a+1):length(str)],collapse=' '), collapse='')
    }
  })
  enrdat$label[which(enrdat$Disease != 'Lung Adenocarcinoma' | enrdat$Gene != 'KEAP1')] <- ''
  f4a <- ggplot() + geom_point(aes(x=-log10(P_Value_VUS),y=-log10(P_Value_Driver), color=Gene), data = enrdat) + theme_bw() +
    geom_hline(yintercept = -log10(.05), linetype = 'dashed') + geom_vline(xintercept = -log10(.05), linetype = 'dashed') + 
    ggtitle('Enrichment in Drivers vs. VUSs') + geom_text(aes(x=-log10(P_Value_VUS),y=-log10(P_Value_Driver),
          label=label),hjust=.7,vjust=-0.5, size=2.5, data=enrdat) + ylab('Log P-Value (Drivers)') + ylab('Log P-Value (VUS)') + 
    theme(plot.title = element_text(size = mainsize))
  
  #Biallelic rate of KEAP1 driver, VUS across diseases
  KEAP1_zygosity <- as.data.frame(filtered[which(filtered$Hugo_Symbol == 'KEAP1' & !filtered$zygosity_call %in% c('Heterozygous - LOH (gene-level)',
             'Indeterminate') & filtered$mutation_class != 'silent'),])
  KEAP1_zygosity$zygosity <- 'wt'
  KEAP1_zygosity$zygosity[which(KEAP1_zygosity$zygosity_call %in% c('Heterozygous','Gain-of-mutant - Mut + copy gain'))] <- 'het'
  KEAP1_zygosity$zygosity[which(KEAP1_zygosity$zygosity_call %in% c('Biallelic - Homdel','Biallelic - Mut + LOH','Biallelic - Mut + fusion','Biallelic - compound'))] <- 'biallelic'
  geneTab.vus <- as.data.frame.matrix(table(KEAP1_zygosity[which(KEAP1_zygosity$mutation_class == 'VUS' & KEAP1_zygosity$zygosity != 'wt'),c('disease_subtype','zygosity')]))
  geneTab.vus$disease <- rownames(geneTab.vus)
  geneTab.vus$Class <- 'VUS'
  geneTab.driver <- as.data.frame.matrix(table(KEAP1_zygosity[which(KEAP1_zygosity$mutation_class == 'Driver'),c('disease_subtype','zygosity')]))
  geneTab.driver$Class <- 'Driver'
  geneTab.driver$disease <- rownames(geneTab.driver)
  geneTab <- rbind(geneTab.driver, geneTab.vus)
  geneTab$mut <- rowSums(geneTab[,c('biallelic','het')])
  geneTab$biallelic_fraction <- geneTab[,'biallelic']/geneTab[,'mut']
  merged <- merge(geneTab.driver, geneTab.vus, by='disease')[,c(1,2:3,5:6)]
  mtab <- rowSums(merged[,2:5])
  names(mtab) <- merged$disease
  geneTab <- geneTab[which(geneTab$disease %in% names(mtab)[which(mtab >= 10)]),]
  geneTab <- unique(cbind(geneTab, t(sapply(1:nrow(geneTab),function(i){prop.test(x=geneTab$biallelic[i],n=geneTab$mut[i])$conf.int}))))
  colnames(geneTab)[(ncol(geneTab) - 1):(ncol(geneTab))] <- c('Lower','Upper')
  geneTab <- geneTab[order(geneTab$biallelic_fraction, decreasing = TRUE),]
  geneTab$disease <- factor(geneTab$disease, levels = unique(geneTab$disease))
  f4b <- ggplot(aes(x=disease,y=biallelic_fraction,fill=Class), data = geneTab) + ylim(0,1) + 
    geom_bar(stat = 'identity',position = 'dodge') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5), plot.title = element_text(size = mainsize)) + 
    scale_fill_manual(values = c('Driver'='gray47','VUS'='gray')) + geom_errorbar(stat = 'identity',data=geneTab, mapping=aes(ymin=Upper, ymax=Lower),
               position_dodge(.75,preserve='single'), width=0.2, size=1) + xlab('Disease') + ylab('Biallelic Fraction') + ggtitle('KEAP1 Across Cancers')
  
  ###CUP Bar plot:
  if(is.null(filtered$loh)){
    filtered$loh <- NA
    filtered$loh[which(filtered$zygosity_call == 'Biallelic - Mut + LOH' | filtered$het_loh)] <- TRUE
    filtered$loh[which((filtered$zygosity_call == 'wt' & filtered$het_loh) | (filtered$zygosity_call %in% 
        c('Biallelic - compound','Gain-of-mutant - Mut + copy gain','Heterozygous','Amplification - CNA','Amplification - (gene-level)')))] <- FALSE
  }
  if(is.null(filtered$loh)){filtered$loh <- filtered$zygosity_call %in% c('Biallelic - Mut + LOH')}
  tab2 <- as.data.frame(table(filtered[which(filtered$Hugo_Symbol == 'KEAP1' & filtered$disease_subtype == 'Cancer of Unknown Primary'),c('loh','mutation_class')]))
  tab2$N <- sapply(1:nrow(tab2), function(i){sum(tab2$Freq[which(tab2$mutation_class == tab2$mutation_class[i])])})
  tab2$Fraction <- tab2$Freq/tab2$N
  tab2 <- cbind(tab2, t(sapply(1:nrow(tab2),function(i){prop.test(x=tab2$Freq[i],n=tab2$N[i])$conf.int})))
  p <- prop.test(x=tab2$Freq[c(2,4)],n=tab2$N[c(2,4)])$p.value
  colnames(tab2)[c(1,2,6,7)] <- c('Biallelic','Class','Lower','Upper')
  tab2$Class <- as.character(tab2$Class)
  tab2$Class <- sapply(1:nrow(tab2),function(i){paste(as.character(tab2$Class[i]),' (n=',tab2$N[i],')',sep='')})
  tab2$Class <- factor(tab2$Class, levels = sort(unique(tab2$Class))[c(1,3,2)])
  abbreviation <- tcga_subtype_map$tcga[which(tcga_subtype_map$impact == 'Cancer of Unknown Primary')]
  f4c <- ggplot(aes(x=Class,y=Fraction,fill=Biallelic), data=tab2) + geom_bar(stat = 'identity', position = 'stack') + 
    theme_bw() + scale_fill_manual(values=c('TRUE'='gray47','FALSE'='gray87')) + 
    ggtitle(paste(abbreviation,' Biallelic Frac. vs. Mut. Class (p=',round(p,digits = 4),')',sep='')) + theme(legend.position="left") + 
    theme(axis.text.x= element_text(size = 8)) + theme(plot.title = element_text(size = mainsize)) + 
    geom_errorbar(stat = 'identity',data=tab2[c(2,4,6),], mapping=aes(ymin=Upper, ymax=Lower), width=0.2, size=1)
  
  #CUP oncoprint:
  f4d <- biOncoprint(data=filtered, gene='KEAP1', disease='Cancer of Unknown Primary', savePath=NULL, impact_genes, cutoff=10)
  f4e <- classSurvival(filtered, clinicalData, disease='Cancer of Unknown Primary', gene='KEAP1')
  
  ###LUAD Bar plot:
  tab2 <- as.data.frame(table(filtered[which(filtered$Hugo_Symbol == 'KEAP1' & filtered$disease_subtype == 'Lung Adenocarcinoma'),c('loh','mutation_class')]))
  tab2$N <- sapply(1:nrow(tab2), function(i){sum(tab2$Freq[which(tab2$mutation_class == tab2$mutation_class[i])])})
  tab2$Fraction <- tab2$Freq/tab2$N
  tab2 <- cbind(tab2, t(sapply(1:nrow(tab2),function(i){prop.test(x=tab2$Freq[i],n=tab2$N[i])$conf.int})))
  p <- prop.test(x=tab2$Freq[c(2,4)],n=tab2$N[c(2,4)])$p.value
  colnames(tab2)[c(1,2,6,7)] <- c('Biallelic','Class','Lower','Upper')
  tab2$Class <- as.character(tab2$Class)
  tab2$Class <- sapply(1:nrow(tab2),function(i){paste(as.character(tab2$Class[i]),' (n=',tab2$N[i],')',sep='')})
  tab2$Class <- factor(tab2$Class, levels = sort(unique(tab2$Class))[c(1,3,2)])
  f4f <- ggplot(aes(x=Class,y=Fraction,fill=Biallelic), data=tab2) + geom_bar(stat = 'identity', position = 'stack') + 
    theme_bw() + scale_fill_manual(values=c('TRUE'='gray47','FALSE'='gray87')) + 
    ggtitle(paste(abbreviation,' Biallelic Frac. vs. Mut. Class (p=',round(p,digits = 4),')',sep='')) + theme(legend.position="left") + 
    theme(axis.text.x= element_text(size = 8)) + theme(plot.title = element_text(size = mainsize)) + 
    geom_errorbar(stat = 'identity',data=tab2[c(2,4,6),], mapping=aes(ymin=Upper, ymax=Lower), width=0.2, size=1)
  
  #LUAD oncoprint:
  f4g <- biOncoprint(data=filtered, gene='KEAP1', disease='Lung Adenocarcinoma', savePath=NULL, impact_genes, cutoff=30)
  
  ###Differential co-mutation: VUS vs. driver.
  ###Differential survival: VUS vs. driver; only drivers
  if(!is.null(stratify)){
    ymutated <- unique(filtered$Tumor_Sample_Barcode[which(filtered$Hugo_Symbol == 'EGFR' & filtered$mutation_class != 'Silent' & 
                                                             filtered$disease_subtype == 'Lung Adenocarcinoma')])
    mutDatA <- filtered[which(filtered$Tumor_Sample_Barcode %in% ymutated),]
    mutDatB <- filtered[which(!filtered$Tumor_Sample_Barcode %in% ymutated),]
    pzgA <- classSurvival(mutDatA, clinicalData, disease='Lung Adenocarcinoma', gene='KEAP1', addendum=paste(stratify, '-mutated'))
    pzgB <- classSurvival(mutDatB, clinicalData, disease='Lung Adenocarcinoma', gene='KEAP1', addendum=paste(stratify, '-WT'))
    f4h <- cowplot::plot_grid(pzgA$plot, pzgB$plot, ncol=1, nrow=2)
  }else{
    f4h <- classSurvival(filtered, clinicalData, disease='Lung Adenocarcinoma', gene='KEAP1')
  }
  
  #Biallelic rates vs. mut. class, TCGA:
  res.tcga.keap1 <- zygosity_and_class_plot(filtered.tcga, gene='KEAP1', disease='LUAD', tcga_map = tcga_subtype_map, dataset='TCGA', mainsize=12)
  f4i <- res.tcga.keap1$plot
  
  #diff expression 
  if(!file.exists(expressionResPath)){
    res <- runDiffExprVUS(expressionData, proteinData, gene='KEAP1', disease='LUAD', mutData.tcga, clinData.tcga,
                          pathwayPath.reactome=reactomePathway,pathwayPath.hallmarks=hallmarksPathway,runGSEA=TRUE)
    if(!is.null(expressionResPath)){save(res, file=expressionResPath)}
  }else{
    res <- get(load(expressionResPath))
  }
  
  #Expression (RNA and protein) plots:
  exprPlots <- generateExpressionPlots2(res$res.vus_v_wt, res$res.vus_v_driver,gene='KEAP1', disease='LUAD', path=NULL, 
                                        counts = res$counts, 'counts.prot'=res$counts.prot)
  f4j <- exprPlots$plots[1]
  f4k <- exprPlots$plots[2]
  f4l <- exprPlots$plots[3]
  f4m <- exprPlots$plots[4]
  
  ###Expression plots section: placement to be worked out later
  #Selection for biallelic in TCGA:
  #res.impact.keap1 <- zygosity_and_class_plot(filtered, gene='KEAP1', disease='Lung Adenocarcinoma', dataset='impact')
  #keap1.luad <- ggarrange(res.impact.keap1$plot, res.tcga.keap1$plot, ncol=2, labels = c('J','K'))
  #dat.keap1.impact <- res.impact.keap1$data
  #dat.keap1.tcga <- res.tcga.keap1$data
  
  #Combining panels:
  layout <- 
    "
  AAAA#BBBB
  CDDDDDDEE
  FGGGGGGHH
  IJKLLLLMM
  "
  f4 <- f4a + f4b + f4c + f4d + f4e$plot + f4f + f4g + f4h$plot + f4i + f4j + f4k + f4l + f4m + 
    plot_layout(design = layout) + plot_annotation(tag_levels = 'A')
  pdf(paste(plotPath,'/f4.pdf',sep=''),width=30,height=30)
  print(f4)
  dev.off() 
  
  #Lollipop plots:
  m1.cup <- read.maf(mutData.impact[which(mutData.impact$Hugo_Symbol == 'KEAP1' & mutData.impact$Tumor_Sample_Barcode %in% 
                                            filtered$Tumor_Sample_Barcode[which(filtered$disease_subtype == 'Cancer of Unknown Primary' &  filtered$mutation_class == 'Driver')]),])
  m2.cup <- read.maf(mutData.impact[which(mutData.impact$Hugo_Symbol == 'KEAP1' & mutData.impact$Tumor_Sample_Barcode %in% 
                                            filtered$Tumor_Sample_Barcode[which(filtered$disease_subtype == 'Cancer of Unknown Primary' & filtered$mutation_class == 'VUS')]),])
  m1.luad <- read.maf(mutData.impact[which(mutData.impact$Hugo_Symbol == 'KEAP1' & mutData.impact$Tumor_Sample_Barcode %in% 
                                             filtered$Tumor_Sample_Barcode[which(filtered$disease_subtype == 'Lung Adenocarcinoma' & filtered$mutation_class == 'Driver')]),])
  m2.luad <- read.maf(mutData.impact[which(mutData.impact$Hugo_Symbol == 'KEAP1' & mutData.impact$Tumor_Sample_Barcode %in% 
                                             filtered$Tumor_Sample_Barcode[which(filtered$disease_subtype == 'Lung Adenocarcinoma' & filtered$mutation_class == 'VUS')]),])
  
  png('/Users/zuckerm/Documents/biallelic/analysis01/figures/f4n.png',width=600, height = 300)
  par(mfrow=c(1,1))
  lollipopPlot2(m1 = m1.cup, m2 = m2.cup,gene = 'KEAP1', m1_name = 'Drivers',m2_name = 'VUSs')
  dev.off()
  png('/Users/zuckerm/Documents/biallelic/analysis01/figures/f4o.png',width=600, height = 300)
  lollipopPlot2(m1 = m1.luad, m2 = m2.luad,gene = 'KEAP1', m1_name = 'Drivers',m2_name = 'VUSs')
  dev.off()
}

generateFigure4(drivers, filtered, biallelicEnrichment.driver, biallelicEnrichment.vus, mutData.tcga, mainsize=10, stratify=NULL)

