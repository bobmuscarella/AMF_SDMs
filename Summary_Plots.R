
#######################################
#######################################
#######################################
##### TARGET BACKGROUND #####
#######################################
#######################################
#######################################

#######################################
##### BARPLOTS OF SELECTED MODELS #####
#######################################

df <- readRDS('E:/Bob/AMF_SDMs/RESULTS/20180719_targetBG_summary_table.RDS')
dflist <- split(df, df$taxa)

selected_mod <- unlist(lapply(dflist, function(x) x$model.name[x$min.AICc==min(x$min.AICc)]))
rank <- df$rank[seq(1, nrow(df),3)]

selected_mod_prop <- 100 * sweep(table(selected_mod, rank), 2, t(table(rank)), '/')


pdf('E:/Bob/AMF_SDMs/FIGURES/AIC_models_targetBG_barplot.pdf', width=9, height=4)

par(xpd=T, mfrow=c(1,2), mar=c(6,6,4,2))
cols <- RColorBrewer::brewer.pal(3,'Dark2')
b <- barplot(selected_mod_prop, col=cols, las=2, xlim=c(0.25,7.25))
legend('right', horiz = F, legend=c('All','Climate-only','Resource-only'), 
       bty='n', pch=22, pt.bg=cols, pt.cex=2, inset=-.15)
axis(1, labels=F, at=b)
mtext('Proportion of taxa', 2, 2.5, font=2)
mtext('Taxa rank', 1, at=2.5, 4.25, font=2)
mtext('AICc-selected models
(target background)', 3, at=2.5, 1, font=2)


x <- table(selected_mod, rank)
x <- x[,4:1]
cols2 <- (RColorBrewer::brewer.pal(4,'Accent'))
b <- barplot(t(x), col=cols2, las=2)
legend('top', horiz = F, legend=c('Order','Family','Genus','OTU'), 
       bty='n', pch=22, pt.bg=rev(cols2), pt.cex=2)
axis(1, labels=F, at=b)
mtext('Number of taxa', 2, 2.5, font=2)
mtext('Taxa rank', 1, at=mean(b), 4.25, font=2)
mtext('AICc-selected models
(target background)', 3, at=mean(b), 1, font=2)


dev.off()



###########################################
##### BOXPLOTS OF VARIABLE IMPORTANCE #####
###########################################

# Get a subset of the data frame with only the selected models
x <- df[match(paste(names(selected_mod), selected_mod), paste(df$taxa, df$model.name)),]

xa <- x[x$model.name=='all',]
xc <- x[x$model.name=='climate',]
xr <- x[x$model.name=='resources',]

cols <- c(RColorBrewer::brewer.pal(8,'Dark2'), RColorBrewer::brewer.pal(9,'Set1'), RColorBrewer::brewer.pal(8,'Set2'))

pdf('E:/Bob/AMF_SDMs/FIGURES/AIC_models_targetBG_boxplot.pdf', width=9)

varnames <- gsub('pi.','',names(xa)[grepl('pi.', names(xa))])
par(mfcol=c(2,3), oma=c(1,2,2,1), mar=c(6,2,0.25,1))
boxplot(xa[,grepl('pi.', names(xa))], las=2, col=cols, names=varnames)
mtext('All-variable models', 3, 0.5, font=2)
mtext('Permutation importance', 2, 2.5, font=2)
boxplot(xa[,grepl('vi.', names(xa))], las=2, col=cols, names=varnames)
mtext('Variable contribution', 2, 2.5, font=2)

boxplot(xc[,grepl('pi.', names(xc))], las=2, col=cols, names=varnames)
mtext('Climate-only models', 3, 0.5, font=2)
boxplot(xc[,grepl('vi.', names(xc))], las=2, col=cols, names=varnames)

boxplot(xr[,grepl('pi.', names(xr))], las=2, col=cols, names=varnames)
mtext('Resource-only models', 3, 0.5, font=2)
boxplot(xr[,grepl('vi.', names(xr))], las=2, col=cols, names=varnames)


dev.off()





#######################################
#######################################
#######################################
##### RANDOM BACKGROUND #####
#######################################
#######################################
#######################################

#######################################
##### BARPLOTS OF SELECTED MODELS #####
#######################################

df <- readRDS('E:/Bob/AMF_SDMs/RESULTS/20180719_randomBG_summary_table.RDS')
dflist <- split(df, df$taxa)

selected_mod <- unlist(lapply(dflist, function(x) x$model.name[x$min.AICc==min(x$min.AICc)]))
rank <- df$rank[seq(1, nrow(df),3)]

selected_mod_prop <- 100 * sweep(table(selected_mod, rank), 2, t(table(rank)), '/')


pdf('E:/Bob/AMF_SDMs/FIGURES/AIC_models_targetBG_barplot.pdf', width=9, height=4)

par(xpd=T, mfrow=c(1,2), mar=c(6,6,4,2))
cols <- RColorBrewer::brewer.pal(3,'Dark2')
b <- barplot(selected_mod_prop, col=cols, las=2, xlim=c(0.25,7.25))
legend('right', horiz = F, legend=c('All','Climate-only','Resource-only'), 
       bty='n', pch=22, pt.bg=cols, pt.cex=2, inset=-.15)
axis(1, labels=F, at=b)
mtext('Proportion of taxa', 2, 2.5, font=2)
mtext('Taxa rank', 1, at=2.5, 4.25, font=2)
mtext('AICc-selected models
      (target background)', 3, at=2.5, 1, font=2)


x <- table(selected_mod, rank)
x <- x[,4:1]
cols2 <- (RColorBrewer::brewer.pal(4,'Accent'))
b <- barplot(t(x), col=cols2, las=2)
legend('top', horiz = F, legend=c('Order','Family','Genus','OTU'), 
       bty='n', pch=22, pt.bg=rev(cols2), pt.cex=2)
axis(1, labels=F, at=b)
mtext('Number of taxa', 2, 2.5, font=2)
mtext('Taxa rank', 1, at=mean(b), 4.25, font=2)
mtext('AICc-selected models
      (target background)', 3, at=mean(b), 1, font=2)


dev.off()



###########################################
##### BOXPLOTS OF VARIABLE IMPORTANCE #####
###########################################

# Get a subset of the data frame with only the selected models
x <- df[match(paste(names(selected_mod), selected_mod), paste(df$taxa, df$model.name)),]

xa <- x[x$model.name=='all',]
xc <- x[x$model.name=='climate',]
xr <- x[x$model.name=='resources',]

cols <- c(RColorBrewer::brewer.pal(8,'Dark2'), RColorBrewer::brewer.pal(9,'Set1'), RColorBrewer::brewer.pal(8,'Set2'))

pdf('E:/Bob/AMF_SDMs/FIGURES/AIC_models_randomBG_boxplot.pdf', width=9)

varnames <- gsub('pi.','',names(xa)[grepl('pi.', names(xa))])
par(mfcol=c(2,3), oma=c(1,2,2,1), mar=c(6,2,0.25,1))
boxplot(xa[,grepl('pi.', names(xa))], las=2, col=cols, names=varnames)
mtext('All-variable models', 3, 0.5, font=2)
mtext('Permutation importance', 2, 2.5, font=2)
boxplot(xa[,grepl('vi.', names(xa))], las=2, col=cols, names=varnames)
mtext('Variable contribution', 2, 2.5, font=2)

boxplot(xc[,grepl('pi.', names(xc))], las=2, col=cols, names=varnames)
mtext('Climate-only models', 3, 0.5, font=2)
boxplot(xc[,grepl('vi.', names(xc))], las=2, col=cols, names=varnames)

boxplot(xr[,grepl('pi.', names(xr))], las=2, col=cols, names=varnames)
mtext('Resource-only models', 3, 0.5, font=2)
boxplot(xr[,grepl('vi.', names(xr))], las=2, col=cols, names=varnames)


dev.off()


