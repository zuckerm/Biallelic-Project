library(ggplot2)
library(patchwork)
library(stringr)
library(ggnewscale)
library(RColorBrewer)

#Path where data is located:
dataPath <- '/Users/zuckerm/Documents/biallelic/data'
#path where figures are to be saved
plotPath <- '/Users/zuckerm/Documents/biallelic/analysis01/figures'

#Objects needed: drivers, fga_and_tmb, tsgdf, subtypeTable:
#Load driver-filtered zygosity call file: 
drivers <- get(load(paste(dataPath,'/zygosityDataDrivers.rda',sep='')))
subtypeTable <- read.table(paste(dataPath,'/subtypeTable.txt',sep=''), sep='\t',header=T,comment.char = '&')
fga_and_tmb <- read.table(paste(dataPath,'/fga_and_tmb.txt',sep=''),sep='\t')
tsgdf <- read.table(paste(dataPath,'/custom_gene_list_autosomes.txt',sep=''),sep='\t')

###Plotting functions used in making the panels:
#Driver/VUS zygosity barplot (f1a):
zygosityBarplot <- function(drivers, colorvec, mainsize=12, xsize=6){
  tabf1a <- as.data.frame(table(drivers[which(drivers$mechanism %in% names(colorvec) & drivers$driver_oncogenic),c('mechanism','Gene_Type')]))
  tabf1a$N <- sapply(1:nrow(tabf1a),function(i){sum(tabf1a$Freq[which(tabf1a$Gene_Type == tabf1a$Gene_Type[i])])})
  tabf1a$Fraction <- tabf1a$Freq/tabf1a$N
  colnames(tabf1a)[1] <- c('Category')
  tabf1a$Category <- factor(as.character(tabf1a$Category), levels = names(colorvec))
  tabf1a$Gene_Type <- as.character(tabf1a$Gene_Type)
  nTSGs <- length(unique(drivers$Hugo_Symbol[which(drivers$Gene_Type == 'TSG')]))
  nOGs <- length(unique(drivers$Hugo_Symbol[which(drivers$Gene_Type == 'Oncogene')]))
  tabf1a$Gene_Type[which(tabf1a$Gene_Type == 'TSG')] <- paste('TSG (',nTSGs,')',sep='')
  tabf1a$Gene_Type[which(tabf1a$Gene_Type == 'Oncogene')] <- paste('Oncogene (',nOGs,')',sep='')
  tabf1a$Gene_Type <- factor(tabf1a$Gene_Type, levels = unique(tabf1a$Gene_Type))
  f1a <- ggplot(aes(x=Gene_Type,y=Fraction,fill=Category), data=tabf1a) + 
    geom_bar(stat = 'identity', position = 'stack') + theme_bw() + scale_fill_manual(values=colorvec[which(!names(colorvec) %in% 
                                                                                                             c('LOH','Mut + fusion','Mutation + Gain of Mutant'))]) + xlab('') + 
    ggtitle('Alteration zygosity across drivers') + xlab('') + theme(legend.position="left") + 
    theme(axis.text.x= element_text(size = xsize)) + theme(plot.title = element_text(size = mainsize))
  list('plot'=f1a, 'table'=tabf1a)
}

#Mechanism barplot (f1c):
mechanismBarplot <- function(data, geneset, mainsize=12, filterMutations=NULL,bicats=c('Compound','Mut + Fusion','Homdel','Mut + LOH'), subtypeTable, colorvec){
  dat <- data
  dat$mechanism[which(dat$mechanism %in% c('Mut + fusion'))] <- 'Mut + Fusion'
  tsgs.biallelic <- dat[which(dat$Hugo_Symbol %in% geneset & dat$final_zygosity == 'biallelic'),]
  pancan <- as.data.frame.matrix(table(tsgs.biallelic[,c('Hugo_Symbol','mechanism')]))
  byDisease <- table(tsgs.biallelic[,c('disease_subtype','Hugo_Symbol','mechanism')])
  
  asdf <- as.data.frame(Reduce(rbind, lapply(1:dim(byDisease)[1], function(i){byDisease[i,,]})))
  asdf$Disease <- unlist(lapply(dimnames(byDisease)[[1]], function(x){rep(x,dim(byDisease)[2])}))
  asdf$Gene <- rep(dimnames(byDisease)[[2]],dim(byDisease)[1])
  asdf$n <- rowSums(asdf[,bicats])
  if(!is.null(filterMutations)){
    asdf <- asdf[which(asdf$n >= filterMutations),]
  }
  fracdf <- cbind(asdf[,c('Gene','Disease')], asdf[,bicats]/asdf$n, asdf$n)
  colnames(fracdf)[ncol(fracdf)] <- 'biallelic_events'
  fracdf$pvalue <- sapply(1:nrow(fracdf),function(i){
    test <- chisq.test(rbind(asdf[i,bicats], pancan[which(rownames(pancan) == asdf$Gene[i]),bicats]))
    if(class(test) != 'try-error'){test$p.value}else{NA}
  })
  fracdf$pvalue.corrected <- p.adjust(fracdf$pvalue, method = 'fdr')
  for(i in c(3:ncol(fracdf))){fracdf[,i] <- as.numeric(as.character(fracdf[,i]))}
  fracdf$Significant <- fracdf$pvalue.corrected < .05
  fracdf <- cbind(fracdf, t(sapply(1:nrow(fracdf),function(i){
    index <- which.max(fracdf[i,bicats]) + 2
    if(length(index)==1){c(colnames(fracdf)[index], fracdf[i,index])}else{c(NA,NA)}
  })))
  colnames(fracdf)[(ncol(fracdf)-1):ncol(fracdf)] <- c('Mechanism','Fraction of Biallelic Events')
  fracdf <- merge(fracdf, data.frame('Disease'=subtypeTable$Disease,'Color'=subtypeTable$organ_color), by='Disease')
  for(i in c(3:(2+length(bicats)),length(bicats)+3,length(bicats)+8)){fracdf[,i] <- as.numeric(as.character(fracdf[,i]))}
  
  fracdf2 <- as.data.frame(Reduce(rbind, lapply(1:ncol(pancan),function(i){as.data.frame(pancan[,i]/rowSums(pancan))})))
  fracdf2$Gene <- rep(rownames(pancan),length(bicats))
  fracdf2$Mechanism <- unlist(lapply(colnames(pancan),function(x){rep(x,nrow(pancan))}))
  colnames(fracdf2)[1] <- 'Fraction'
  fracdf$Disease <- factor(fracdf$Disease, levels = unique(fracdf$Disease))
  fracdf$Dominant_Mechanism <- sapply(fracdf$Gene, function(g){
    temp <- fracdf2[which(fracdf2$Gene ==g),]
    temp$Mechanism[which.max(temp$Fraction)]
  })
  fracdf <- fracdf[which(fracdf$pvalue.corrected < .05 & fracdf$Mechanism != fracdf$Dominant_Mechanism),]
  fracdf2 <- fracdf2[which(fracdf2$Gene %in% fracdf$Gene),]
  fracdf$Gene <- factor(fracdf$Gene, levels = geneset)
  fracdf2$Gene <- factor(fracdf2$Gene, levels = geneset)
  
  ###Adding N's:
  dtab <- table(unique(data[,c('Tumor_Sample_Barcode','disease_subtype')])[,2])
  ndf <- data.frame('Disease'=names(dtab), 'n'=as.numeric(dtab))
  ndf$disease_n <- sapply(1:nrow(ndf),function(i){paste(ndf$Disease[i],' (n=',ndf$n[i],')',sep='')})
  fracdf <- merge(fracdf, ndf, by='Disease')
  fracdf$disease_n <- sapply(fracdf$disease_n, function(d){
    str <- strsplit(as.character(d), split = " ")[[1]]
    chars <- length(strsplit(as.character(d), split = '')[[1]])
    if(chars >= 30){
      midpoint <- round(length(str)/2)
      output <- paste(paste(str[1:midpoint],collapse=" "), "\n", paste(str[(midpoint + 1):length(str)],collapse=" "),sep='')
    }else{output <- paste(str, collapse = " ")}
    output
  })
  
  fracdf <- fracdf[order(match(fracdf$Disease, levels(fracdf$Disease))),]
  fracdf$disease_n <- factor(fracdf$disease_n, levels = unique(fracdf$disease_n))
  fracdf <- fracdf[which(fracdf$biallelic_events >= 10),]
  fracdf2 <- fracdf2[which(fracdf2$Gene %in% fracdf$Gene),]
  diseaseColors <- unique(fracdf[,c('disease_n', 'Color')])$Color
  names(diseaseColors) <- unique(fracdf[,c('disease_n', 'Color')])$disease_n
  
  ddf <- rbind(as.matrix(fracdf[,c('Mut + LOH','Gene','Disease','disease_n','n')]),as.matrix(fracdf[,c('Homdel','Gene','Disease','disease_n','n')]),
               as.matrix(fracdf[,c('Mut + Fusion','Gene','Disease','disease_n','n')]),as.matrix(fracdf[,c('Compound','Gene','Disease','disease_n','n')]))
  colnames(ddf)[1] <- 'Fraction'
  ddf <- as.data.frame(ddf)
  ddf$Mechanism <- unlist(lapply(c('Mut + LOH','Homdel','Mut + Fusion','Compound'),function(x){rep(x,nrow(fracdf))}))
  pcdf <- fracdf2
  pcdf$Disease <- 'Pan-Cancer'
  pcdf$disease_n <- pcdf$Disease
  pcdf$n <- NA
  pcdf$disease_n <- sapply(pcdf$disease_n, function(d){
    str <- strsplit(as.character(d), split = " ")[[1]]
    chars <- length(strsplit(as.character(d), split = '')[[1]])
    if(chars >= 30){
      midpoint <- round(length(str)/2)
      output <- paste(paste(str[1:midpoint],collapse=" "), "\n", paste(str[(midpoint + 1):length(str)],collapse=" "),sep='')
    }else{output <- paste(str, collapse = " ")}
    output
  })
  pcdf <- pcdf[order(pcdf$Gene),]
  ddf <- ddf[,c(colnames(pcdf))]
  ddf <- ddf[order(ddf$Gene),]
  df <- rbind(ddf, pcdf)
  df$Fraction <- as.numeric(df$Fraction)
  df <- unique(df)
  df$disease_n <- factor(df$disease_n, levels=names(sort(table(df$disease_n), decreasing = TRUE)))
  diseaseColors <- c(diseaseColors, 'Pan-Cancer'='black')
  df$organ <- subtypeTable$organ[match(df$Disease, subtypeTable$Disease)]
  df$organ_n <- sapply(1:nrow(df),function(i){sum(as.numeric(df$n)[which(df$organ == df$organ[i])])})
  
  ###For making axis labels:
  temp <- merge(unique(df[,c('Disease','n')]), subtypeTable, by='Disease')
  temp$abbrev_n <- sapply(1:nrow(temp),function(i){paste(temp$CTD_ABBR[i], ' (',temp$n[i],')', sep='')})
  temp$abbrev_n[which(temp$Disease == 'Pan-Cancer')] <- 'Pan-Cancer'
  tobind <- temp[,c('Disease','abbrev_n','organ_color')]
  tobind <- rbind(tobind, c('Pan-Cancer','Pan-Cancer','#000000'))
  df <- merge(df, tobind, by='Disease')
  df$n[which(df$Disease == 'Pan-Cancer')] <- -1
  df <- df[order(df$Gene, as.numeric(df$n), decreasing = TRUE),]
  df$abbrev_n <- factor(df$abbrev_n, levels = c(unique(df$abbrev_n)[-which(unique(df$abbrev_n)=='Pan-Cancer')],'Pan-Cancer'))
  df$organ <- factor(df$organ, levels = unique(df$organ))
  df$Gene <- factor(df$Gene, levels = unique(df$Gene))
  #temp <- unique(na.omit(df[,c('organ','abbrev_n','organ_color')]))
  #for(i in 1:3){temp[,i] <- as.character(temp[,i])}
  #temp <- Reduce(rbind, lapply(unique(na.omit(temp$organ)),function(x){rbind(temp[which(temp$organ == x),],c('Pan-Cancer','Pan-Cancer','#000000'))}))
  
  p <- ggplot() + geom_bar(aes(x=abbrev_n,y=Fraction, fill=Mechanism),data=df, position = 'stack',stat = 'identity') + 
    scale_fill_manual(values = colorvec[which(names(colorvec) %in% df$Mechanism)]) + theme(axis.title.x = element_text(size=4)) + 
    theme(axis.text.x = element_text(size=6,angle=45, hjust=1), axis.title.x = element_blank()) +
    #geom_text(data=temp,aes(x=abbrev_n,y=-1.5,colour=organ_color,label=abbrev_n),angle=45,size=4,vjust=1,hjust=1) + theme(legend.position='none') +
    #coord_cartesian(ylim=c(0,1), clip = 'off') + theme(plot.margin = unit(c(1,1,9.7,1.5), "lines")) + 
    facet_wrap(~ Gene, scales = "free_x",ncol = length(unique(df$Gene))) 
  p
  
  list('plot'=p, 'table'=fracdf)
}

#Tileplot (1b):
biallelicRateTileplot <- function(data, xsize=6, mainsize=12, force_gene=TRUE, force_disease=TRUE, path=NULL, colorvec,
                                  justLOH=FALSE, cutoff=0, borders=FALSE, justDrivers=TRUE, fga_and_tmb, tsgdf, subtypeTable, ncolors=4){
  dat <- data
  genes <- tsgdf$Hugo_Symbol
  if(!justLOH){dat$final_zygosity[which(dat$zygosity_call == 'Heterozygous - LOH (gene-level)')] <- 'wt'}
  dat$final_zygosity[which(is.na(dat$final_zygosity))] <- 'wt'
  dat <- as.data.frame(table(dat[,c('disease_subtype','Hugo_Symbol','final_zygosity')]))
  for(i in 1:3){dat[,i] <- as.character(dat[,i])}
  dat$final_zygosity[which(is.na(dat$final_zygosity))] <- 'wt'
  dat <- cbind(dat[which(dat$final_zygosity == 'wt'),c(1,2,4)], dat[which(dat$final_zygosity == 'het'),c(4)], 
               dat[which(dat$final_zygosity == 'biallelic'),4])
  colnames(dat) <- c('Disease','Gene','WT','Het','Biallelic')
  dat$`Biallelic Fraction` <- dat$Biallelic/(dat$Het + dat$Biallelic)
  dat$`Alteration Fraction` <- (dat$Biallelic + dat$Het)/rowSums(dat[,c(3:5)])
  dat <- dat[which(dat$Disease %in% subtypeTable$Disease & dat$Gene %in% genes),]
  if(force_disease){dat$Disease <- factor(dat$Disease, levels = subtypeTable$Disease)}else{
    dat$Disease <- factor(dat$Disease, levels = intersect(subtypeTable$Disease,dat$Disease))}
  if(force_gene){dat$Gene <- factor(dat$Gene, levels = genes)}else{dat$Gene <- factor(dat$Gene, levels = intersect(genes,dat$Gene))}
  dat <- dat[which(dat$`Alteration Fraction` > cutoff & !is.na(dat$`Alteration Fraction`) & !is.na(dat$`Biallelic Fraction`)),]
  grid <- expand.grid(list(unique(dat$Disease), unique(dat$Gene)))
  grid <- cbind(grid, rep(0, nrow(grid)),rep(0, nrow(grid)),rep(0, nrow(grid)),rep(0, nrow(grid)),rep(0, nrow(grid)))
  colnames(grid) <- colnames(dat)
  dat$dg <- sapply(1:nrow(dat),function(i){paste(dat$Disease[i],dat$Gene[i],sep='_')})
  grid$dg <- sapply(1:nrow(grid),function(i){paste(grid$Disease[i],grid$Gene[i],sep='_')})
  grid <- grid[which(!grid$dg %in% dat$dg),]
  dat <- rbind(dat, grid)
  dat$`Biallelic %` <- round(100*(dat$`Biallelic Fraction`))
  
  ###Making table with row for each gene, disease, and mechanism (and one for pancan):
  n <- as.data.frame(table(data$Hugo_Symbol))
  n.drivers <- as.data.frame(table(data$Hugo_Symbol[which(data$mechanism != 'LOH')]))
  n.biallelic <- as.data.frame(table(data$Hugo_Symbol[which(data$final_zygosity == 'biallelic' & data$mechanism != 'LOH')]))
  colnames(n)[1] <- colnames(n.drivers)[1] <- colnames(n.biallelic)[1] <- 'Hugo_Symbol'
  colnames(n)[2] <- 'N'
  colnames(n.drivers)[2] <- 'Drivers'
  colnames(n.biallelic)[2] <- 'Biallelic'
  genedf <- merge(merge(n, n.drivers, by='Hugo_Symbol'), n.biallelic, by='Hugo_Symbol')
  genedf$`Alteration Fraction` <- genedf$Drivers/genedf$N
  genedf$`Biallelic Fraction` <- genedf$Biallelic/genedf$Drivers
  genedf <- merge(genedf, tsgdf[,c('Hugo_Symbol','Pathway')], all = TRUE)
  n.mech <- as.data.frame(table(data[which(data$mechanism != 'LOH' & data$Hugo_Symbol %in% tsgdf$Hugo_Symbol),c('Hugo_Symbol','mechanism')]))
  n.drivers.mech <- as.data.frame(table(data[which(data$mechanism != 'LOH' & data$Hugo_Symbol %in% tsgdf$Hugo_Symbol),c('Hugo_Symbol','mechanism')]))
  colnames(n.mech)[1] <- colnames(n.drivers.mech)[1] <- 'Hugo_Symbol'
  colnames(n.mech)[2] <- colnames(n.drivers.mech)[2] <- 'Mechanism'
  colnames(n.mech)[3] <- 'Instances'
  colnames(n.drivers.mech)[3] <- 'Drivers'
  mechdf <- merge(n.mech, n.drivers.mech, by=c('Hugo_Symbol','Mechanism'))
  mechdf <- merge(mechdf, data.frame('Hugo_Symbol'=genedf$Hugo_Symbol, 'Alterations'=genedf$Drivers), by='Hugo_Symbol')
  mechdf$`Fraction_of_Alterations` <- mechdf$Drivers/mechdf$Alterations
  genedf <- genedf[order(genedf$`Biallelic Fraction`, decreasing = TRUE),]
  mechdf$Hugo_Symbol <- factor(mechdf$Hugo_Symbol, levels = tsgdf$Hugo_Symbol)
  #mechdf$Mechanism <- factor(mechdf$Mechanism, levels = 
  #                              c('Homdel','Mut + LOH','Mut + fusion','Compound','Mutation','Amplification'))
  mechdf <- mechdf[which(mechdf$Mechanism != 'Amplification'),]
  mechdf$Mechanism <- factor(mechdf$Mechanism, levels = 
                               c('Mut + LOH','Mut + fusion','Homdel','Compound','Mutation'))
  dat$`Alteration %` <- round(100*dat$`Alteration Fraction`)
  geneOrder <- as.character(genedf$Hugo_Symbol[order(genedf$`Alteration Fraction`, decreasing = TRUE)])
  genedf$Hugo_Symbol <- factor(genedf$Hugo_Symbol, levels = geneOrder)
  mechdf$Hugo_Symbol <- factor(mechdf$Hugo_Symbol, levels = geneOrder)
  #mechdf$Mechanism <- factor(mechdf$Mechanism, levels = names(colorvec))
  dat$Gene <- factor(dat$Gene, levels = geneOrder)
  #Median biallelic fraction: midpoint for color range
  med <- length(which(data$final_zygosity == 'biallelic'))/length(which(data$zygosity_call %in% 
       c('Biallelic - Homdel','Biallelic - Mut + LOH','Biallelic - Mut + fusion','Biallelic - compound','Gain-of-mutant - Mut + copy gain','Heterozygous')))
  
  ###Adding 'n' values:
  dtab <- table(unique(data[,c('Tumor_Sample_Barcode','disease_subtype')])[,2])
  ndf <- data.frame('Disease'=names(dtab), 'n'=as.numeric(dtab))
  ndf$disease_n <- sapply(1:nrow(ndf),function(i){paste(ndf$Disease[i],' (n=',ndf$n[i],')',sep='')})
  dat <- merge(dat, ndf, by='Disease')
  
  #Adding line breaks:
  dat$disease_n <- sapply(dat$disease_n, function(d){
    str <- strsplit(as.character(d), split = " ")[[1]]
    chars <- length(strsplit(as.character(d), split = '')[[1]])
    if(chars >= 30){
      midpoint <- round(length(str)/2)
      output <- paste(paste(str[1:midpoint],collapse=" "), "\n", paste(str[(midpoint + 1):length(str)],collapse=" "),sep='')
    }else{output <- paste(str, collapse = " ")}
    output
  })
  dat <- dat[order(match(dat$Disease, levels(dat$Disease))),]
  dat$disease_n <- factor(dat$disease_n, levels = unique(dat$disease_n))
  genedf$`Log Alteration Fraction` <- log(genedf$`Alteration Fraction`)
  
  #TMB and FGA plots:
  fga_and_tmb_dat <- unique(merge(fga_and_tmb, dat[,c('Disease','disease_n','n')],by='Disease'))
  fga_and_tmb_dat <- fga_and_tmb_dat[which(fga_and_tmb_dat$Tumor_Sample_Barcode %in% data$Tumor_Sample_Barcode & fga_and_tmb_dat$TMB <= 20),]

  #Pathway info:
  dat <- merge(dat, data.frame('Gene'=tsgdf$Hugo_Symbol, 'Pathway'=tsgdf$Pathway), by='Gene')
  genedf <- merge(genedf, data.frame('Hugo_Symbol'=tsgdf$Hugo_Symbol, 'Pathway'=tsgdf$Pathway), by='Hugo_Symbol')
  genedf <- genedf[which(genedf$Hugo_Symbol %in% genes),]
  mechdf <- merge(mechdf, data.frame('Hugo_Symbol'=tsgdf$Hugo_Symbol, 'Pathway'=tsgdf$Pathway), by='Hugo_Symbol')
  dat <- unique(merge(dat, data.frame('Gene'=tsgdf$Hugo_Symbol, 'Pathway'=tsgdf$Pathway), by='Gene'))
  dat <- unique(merge(dat, data.frame('Disease'=fga_and_tmb$Disease,'organ'=fga_and_tmb$organ), by='Disease'))
  dat <- dat[,-which(colnames(dat) == 'Pathway.x')]
  colnames(dat)[which(colnames(dat) == 'Pathway.y')] <- 'Pathway'
  genedf <- genedf[,-which(colnames(genedf) == 'Pathway.x')]
  colnames(genedf)[which(colnames(genedf) == 'Pathway.y')] <- 'Pathway'

  fga_and_tmb_dat <- merge(fga_and_tmb_dat, subtypeTable[,c('ONCOTREE_CODE','organ_color','CTD_ABBR')],by='ONCOTREE_CODE')
  fga_and_tmb_dat$abbrev_n <- str_c(fga_and_tmb_dat$CTD_ABBR, ' (',fga_and_tmb_dat$n,')')
  temp <- unique(fga_and_tmb_dat[,c('organ','disease_n','organ_color','CTD_ABBR','abbrev_n','n')])
  temp$organ_n <- sapply(1:nrow(temp),function(i){sum(temp$n[which(temp$organ == temp$organ[i])])})
  temp <- temp[order(temp$n, decreasing = TRUE),]
  fga_and_tmb_dat$abbrev_n <- factor(fga_and_tmb_dat$abbrev_n, levels = temp$abbrev_n)
  organOrder <- unique(subtypeTable$organ[order(subtypeTable$organ_order)])
  fga_and_tmb_dat$organ <- factor(fga_and_tmb_dat$organ, levels = organOrder)
  temp$organ <- factor(temp$organ, levels = organOrder)
  dat <- merge(dat, temp[,c('disease_n','abbrev_n')], by='disease_n')
  dat$organ <- factor(dat$organ, levels = organOrder)
  
  #Abbrev_n: x-axis:
  temp <- temp[order(temp$organ, decreasing = FALSE),]
  dat$abbrev_n <- factor(dat$abbrev_n, levels = temp$abbrev_n)
  fga_and_tmb_dat$abbrev_n <- factor(fga_and_tmb_dat$abbrev_n, levels = temp$abbrev_n)
  
  #pathway sorting:
  pathway_n <- sapply(unique(genedf$Pathway), function(pw){sum(genedf$Drivers[which(genedf$Pathway == pw)])})
  names(pathway_n) <- unique(genedf$Pathway)
  pathwayOrder <- unique(names(sort(pathway_n, decreasing = TRUE)))
  pathwayOrder <- c(pathwayOrder[which(pathwayOrder != 'unassigned')],'unassigned')
  genedf$`Alteration %` <- 100*genedf$`Alteration Fraction`
  genedf$`Log Alteration %` <- log(100*genedf$`Alteration Fraction`)
  dat$Pathway <- factor(dat$Pathway, levels = pathwayOrder)
  dat$Gene <- factor(dat$Gene, levels = rev(geneOrder))
  genedf$Pathway <- factor(genedf$Pathway, levels = pathwayOrder)
  genedf$Hugo_Symbol <- factor(genedf$Hugo_Symbol, levels = rev(geneOrder))
  mechdf$Pathway <- factor(mechdf$Pathway, levels = pathwayOrder)
  mechdf$Hugo_Symbol <- factor(mechdf$Hugo_Symbol, levels = rev(geneOrder))
  mechdf$Mechanism <- factor(mechdf$Mechanism, levels = c("Mutation","Mutation + Gain of Mutant","Amplification","Compound",
                                                          "Homdel","Mut + fusion","Mut + LOH"))
  
  #Getting rid of gene/diseases with < 4 mutations:
  dat$Alterations <- (as.numeric(as.character(dat$WT)) + as.numeric(as.character(dat$Het)) + as.numeric(as.character(dat$Biallelic)))*as.numeric(
    as.character(dat$`Alteration Fraction`))
  outputDat <- dat
  dat$`Alteration %`[which(dat$`Alteration %` == 0)] <- ''
  dat$`Biallelic %`[which(dat$`Alteration %` == 0)] <- 0
  dat$`Biallelic %`[which(dat$Alterations < 4)] <- 0
  dat$`Alteration %`[which(dat$Alterations < 4)] <- NA
  dat$`Biallelic %`[which(dat$`Alteration %` == '')] <- 0
  
  #Color gradient:
  colfunc <- colorRampPalette(c("white", "yellow", "red"))
  grad <- colfunc(ncolors)
  
  #Adding number of driver alterations in parentheses; because of facets_grtid, have to add it to all three data frames and 
  #make gene_alts the x-variable instead of gene; also have to make sure order is preserved:
  genedf <- genedf[order(genedf$Hugo_Symbol),]
  genedf$gene_alts <- str_c(genedf$Hugo_Symbol, ' (',genedf$Drivers, ')')
  genedf$gene_alts <- factor(genedf$gene_alts, levels = unique(genedf$gene_alts))
  dat$gene_alts <- factor(genedf$gene_alts[match(dat$Gene, genedf$Hugo_Symbol)], levels = unique(genedf$gene_alts))
  mechdf$gene_alts <- factor(genedf$gene_alts[match(mechdf$Hugo_Symbol, genedf$Hugo_Symbol)], levels = unique(genedf$gene_alts))
  
  p_alt <- ggplot(aes(x=gene_alts, y=1000*`Alteration Fraction`), data=genedf) + scale_y_log10(breaks=c(1,5,10,20,30,50,100,1000)) + geom_bar(stat='identity') + 
    theme_bw() + coord_flip() + xlab('') + theme(axis.text.y = element_text(size=9)) + ggtitle('Pan-Cancer') + 
    theme(axis.text.x = element_text(size = xsize), plot.title  = element_text(size=mainsize), axis.title.x = element_text(size = xsize)) + 
    xlab('10*Alteration Fracion') + 
    facet_grid(Pathway ~ ., scales = "free", space = "free") + 
    theme(strip.text = element_text(
      size = 6))
  p_heatmap <- ggplot(data=dat,aes(x=abbrev_n,y=gene_alts,fill=`Biallelic %`)) + geom_tile() + theme_bw() + 
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(), plot.title = element_text(size=mainsize),
          legend.position = 'right') + 
    geom_text(aes(label = `Alteration %`), color = "black", size = 3) + scale_fill_gradientn(colours = grad,values = seq(from=0, to=1, length=ncolors)) + 
    labs(x='',y='Gene') + ylab('') + ggtitle('Biallelic Fraction Across Cancer Type') + 
    geom_point(data=dat[which(dat$Alterations < 4 & dat$Alterations > 0 | dat$`Alteration %` == ''),], shape=5, aes(x=abbrev_n,y=gene_alts,fill=`Biallelic %`)) + 
    facet_grid(Pathway ~ organ, scales = "free", space = "free")
  p_bi <- ggplot(aes(x=gene_alts,y=Fraction_of_Alterations, fill=Mechanism), data=mechdf) + geom_bar(stat = 'identity',position = 'stack') + 
    coord_flip() + theme_bw() + xlab('TSG') + ylab('Biallelic Fraction') + 
    ggtitle('Pan-Cancer') + theme(axis.text.y = element_blank(), plot.title = element_text(size=mainsize), axis.title.x = element_text(size = xsize)) + 
    scale_fill_manual(values=colorvec) + theme(legend.position = "none") + 
    facet_grid(Pathway ~ ., scales = "free", space = "free") + 
    theme(strip.text = element_text(
      size = 6))
  
  tmb <- ggplot(data=fga_and_tmb_dat) + geom_violin(aes(x=abbrev_n, y=TMB)) + theme(axis.text.x = element_blank(),axis.text.y = element_text(size = xsize),axis.title.y = element_text(size = xsize),
                                                                                    axis.title.x = element_blank()) + xlab('Disease') + 
    geom_text(data=temp,aes(x=abbrev_n,y=-1.5,colour=organ_color,label=abbrev_n),angle=45,size=4,vjust=1,hjust=1) + theme(legend.position='none') +
    coord_cartesian(ylim=c(0,20), clip = 'off') + theme(plot.margin = unit(c(1,1,9.7,1.5), "lines")) + 
    facet_grid(. ~ organ, scales = "free", space = "free")
  
  fga <- ggplot() + geom_violin(data=fga_and_tmb_dat, aes(x=abbrev_n, y=FGA)) + theme(axis.text.x = element_blank(),
                                                                                      axis.text.y = element_text(size = xsize),axis.title.y = element_text(size = xsize)) + 
    facet_grid(. ~ organ, scales = "free", space = "free")
  layout <- 
    "
  ABCCC
  ##DDD
  ##EEE
  "
  p <- p_alt + p_bi + p_heatmap + fga + tmb + plot_layout(design = layout, widths = c(1,1,8), heights = c(8,.5,.5))
  
  list('plot'=p,'table'=dat, 'table2'=mechdf)
}

#Generating figure 1:
generateFigure1 <- function(drivers, xsize=6, mainsize=12, fga_and_tmb, tsgdf, subtypeTable){
  #f1a and b:
  colorvec <- c('Amplification'="#D25047",'Mutation'="#CED1CE",'Mutation + Gain of Mutant'="#EDBE59",'Compound'= "#8CA1D3",'Homdel'="#242860", 
                'Mut + Fusion'="#9E65AA",'Mut + LOH'="#58B9F4")
  f1a <- zygosityBarplot(drivers, colorvec=colorvec)
  f1a_table <- f1a$table
  f1a <- f1a$plot
  
  ###
  #Removing amplifications in RUNX1, GATA3, and FOXA1:
  removeAmps <- which(drivers$Hugo_Symbol %in% c('RUNX1','GATA3','FOXA1') & drivers$zygosity_call %in% 
                        c('Amplification - (gene-level)','Amplification - CNA'))
  if(length(removeAmps) > 0){data <- drivers[-removeAmps,]}
  f1b <- biallelicRateTileplot(data=data, path=NULL,fga_and_tmb=fga_and_tmb, tsgdf=tsgdf, subtypeTable=subtypeTable, colorvec=colorvec)
  f1b_table <- f1b$table
  f1b_table.pancan <- f1b$table2
  f1b <- f1b$plot
  
  #f1c:
  f1c <- mechanismBarplot(data=drivers, geneset=tsgdf$Hugo_Symbol, colorvec=colorvec, mainsize=mainsize, subtypeTable=subtypeTable)
  f1c_table <- f1c$table
  colnames(f1c_table)[c(9,12)] <- c('preferred_mechanism','Dominant_Mechanism_Pancan')
  f1c <- f1c$plot
  top <- ggarrange(f1a, f1b, widths = c(120, 500),labels = c('A','B'), ncol=2, nrow = 1)
  f1 <- ggarrange(top, f1c, heights = c(300, 100), nrow=2, ncol=1, labels = c('','C'))
  
  #Save figure:
  pdf(paste(plotPath,'/f1.pdf',sep=''),width = 25, height = 18)
  print(f1)
  dev.off()
  
  #Write tables
  write.table(f1a_table, file=paste(plotPath,'/f1a_table.txt',sep=''), sep='\t')
  write.table(f1b_table, file=paste(plotPath,'/f1b_table.txt',sep=''), sep='\t')
  write.table(f1c_table, file=paste(plotPath,'/f1c_table.txt',sep=''), sep='\t')
}


