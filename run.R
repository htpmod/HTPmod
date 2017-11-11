################################################
# Author: Dijun Chen (chendijun2012@gmail.com)
# Update on Mar. 12, 2015
################################################
pks <- c('ade4','reshape2','fmsb','klaR','gplots','outliers','ape','glmnet','randomGLM','randomForest','nnet','rpart','gbm','plotrix','relaimpo','RColorBrewer')
## install.packages(pks)
## source("https://bioconductor.org/biocLite.R")
## biocLite("pcaMethods")


if(!dir.exists('figures')){
	dir.create('figures')
}
if(!dir.exists('tables')){
	dir.create('tables')
}
################ part 1 ################
rm(list=ls())
library(plotrix)
if(file.exists('my.code.R')){
	source('my.code.R')
}
source('models.R')

KN1121 <- read.csv("data/1121KN.data.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)

raw.control <- KN1121[KN1121$Treatment=='control', ]
raw.stress <- KN1121[KN1121$Treatment=='stress', ]
raw.data <- list(Control=raw.control, Stress=raw.stress, "Control+Stress"=KN1121)

treatment.cols <- c("#e4007f", "#f8b62d", "f39800")
names(treatment.cols) <- names(raw.data)

pdf('figures/Fig.1 - KN1121.biomass.correlation.pdf', width=10, heigh=10, pointsize=10)
ref.cors <- data.frame(dataset=names(raw.data), FW=0, DW=0, stringsAsFactors=FALSE)
rownames(ref.cors) <- names(raw.data)
ref.rmsre <- ref.r2 <- ref.cors 
for(dataset in names(raw.data)){
	plot.data <- raw.data[[dataset]][, c('FW', 'DW', 'volume.fluo.iap')]
	pairs(plot.data, main=dataset, lower.panel=panel.smooth, upper.panel=panel.cor) 
	cors <- cor(plot.data)
	ref.cors[dataset, 'FW'] <- cors['FW', 'volume.fluo.iap']
	ref.cors[dataset, 'DW'] <- cors['DW', 'volume.fluo.iap']
	
	for(bio in c('FW', 'DW')){
		dat <- plot.data[, c(bio, 'volume.fluo.iap')]
		colnames(dat) <- c('x', 'y')
		model <- lm(y ~ x, data=dat)
		predicted <- predict(model)
		observed <- dat[,'y']
		ref.r2[dataset, bio] <- summary(lm(predicted ~ observed))$r.squared
		ref.rmsre[dataset, bio] <- rmsre(observed, predicted)
	}
}
dev.off()
write.table(ref.cors, 'tables/Table 1 - biomass.correlation-ref.cors.csv', sep=',', dec='.', na='', row.names=FALSE)
write.table(ref.r2, 'tables/Table 1 - biomass.correlation-ref.r2.csv', sep=',', dec='.', na='', row.names=FALSE)
write.table(ref.rmsre, 'tables/Table 1 - biomass.correlation-ref.rmsre.csv', sep=',', dec='.', na='', row.names=FALSE)

## normalization 
KN1121[, -c(1:5)] <- apply(KN1121[, -c(1:5)], 2, function(x){x/max(abs(x))})

comm.cols <- colnames(KN1121)[c(1:5)]
hist.cols <- colnames(KN1121)[grep('histogram', colnames(KN1121))]
non.hists <- setdiff(colnames(KN1121), hist.cols)
KN1121.data <- KN1121[, non.hists]
KN1121.hist <- KN1121[, c(comm.cols, hist.cols)]
# KN1121.data[, -c(1:5)] <- apply(KN1121.data[, -c(1:5)], 2, function(x){x/max(abs(x))})
# KN1121.hist[, -c(1:5)] <- apply(KN1121.hist[, -c(1:5)], 2, function(x){x/max(abs(x))})

KN1121.control <- KN1121.data[KN1121.data$Treatment=='control', ]
KN1121.stress <- KN1121.data[KN1121.data$Treatment=='stress', ]
# KN1121.control[, -c(1:5)] <- apply(KN1121.control[, -c(1:5)], 2, function(x){x/max(abs(x))})
# KN1121.stress[, -c(1:5)] <- apply(KN1121.stress[, -c(1:5)], 2, function(x){x/max(abs(x))})

barley.list <- list(Control=KN1121.control, Stress=KN1121.stress, "Control+Stress"=KN1121.data)

barley.traits <- colnames(KN1121.data)[-c(1:5)]
trait.index <- 1:length(barley.traits)
names(trait.index) <- sort(barley.traits)

if(FALSE){
	train.x <- KN1121.data[, -c(1:5)]
	train.y <- KN1121.data[, 4]
	train.x2 <- data.frame(apply(train.x, 2, function(x){x/max(x, na.rm=TRUE)}))
	r1 = svm.regression(train.x, train.y)
	r2 = ml.regression(train.x2, train.y)
	
	train <- cbind(y=train.y, train.x)
	formula <- paste('y~', paste(colnames(train.x), collapse='+'))
	## build a multivariate linear model
	lm.model <- lm(formula, data=train)
	
	library(relaimpo, verbose=FALSE)
	re=calc.relimp(lm.model, type=c("car"), rela=TRUE)
}

perform.regression <- function(models, within.val.data.list, traits){
	reg.table <- c()
	reg.results <- list()
	for(dataset in names(within.val.data.list)){
		pdf(paste0("figures/Fig.1 - ", dataset, ".models.pdf"), width=10, heigh=10, pointsize=10)
		op <- par(mfrow=c(4, 4), pty="m")
		for(biomass in c('FW', 'DW')){
			reg.data <- within.val.data.list[[dataset]][, c(biomass, traits)]
			## reg.data[,1] <- log2(reg.data[, 1])
			## reg.data[,-1] <- apply(reg.data[,-1], 2, function(x){x/max(abs(x))})
			for(model in names(models)){
				if(!is.null(models[[model]])){
					reg <- regression.cross.validation(reg.data, model=models[[model]])
					regression.plot(reg, pch=16, main=paste0(model, ": prediction of ", biomass), 
								xlab="Observed values (g)", ylab="Predicted values (g)")
					reg.results[[dataset]][[biomass]][[model]] <- reg
					if(!is.null(reg)){
						## reg.table <- rbind(reg.table, )
					}else{
						cat("********** Error ! *********\n")
					}
				}else{
					plot(1, type="n", xaxt="n", yaxt="n", frame=TRUE, xlab="", ylab="", main="")
					text(1, "Models", cex=4, col=2)
				}
			}
		}
		par(op)
		dev.off()
	}
	results <- list(table=reg.table, data=reg.results)
	return(results)
}

### test all datasets, for a first look
barley.reg.results <- perform.regression(reg.models, barley.list, barley.traits)

model.cols <- c(MLR="#f39800", MARS="#00a0e9", RF="#e4007f", SVR="#1b9e77")

### test the predictive power of models 
pdf('figures/Fig.1 - KN1121.model.evaluation.pdf', width=10, heigh=10, pointsize=10)
eval.table <- data.frame(stringsAsFactors=FALSE)
times <- 10
for(dataset in names(barley.list)){
	input <- barley.list[[dataset]]
	for(biomass in c('FW', 'DW')){
		op <- par(mfrow=c(4, 4), pty="m")
		reg.data <- input[, c(biomass, barley.traits)]
		## reg.data[,-1] <- apply(reg.data[,-1], 2, function(x){x/max(abs(x))})
		R2.table <- data.frame(matrix(0, ncol=length(reg.models), nrow=times))
		PCC.table <- data.frame(matrix(0, ncol=length(reg.models), nrow=times))
		RMSRE.table <- data.frame(matrix(0, ncol=length(reg.models), nrow=times))
		mu.table <- data.frame(matrix(0, ncol=length(reg.models), nrow=times))
		colnames(R2.table) <- colnames(PCC.table) <- colnames(RMSRE.table) <- colnames(mu.table) <- names(reg.models)
		for(n in 1:times){
			cat("Biomass: ", biomass, "\t ", dataset, "\t ", n, "\n")
			for(model in names(reg.models)){
				reg <- regression.cross.validation(reg.data, model=reg.models[[model]])
				R2.table[n, model] <- reg$R2
				PCC.table[n, model] <- reg$PCC
				RMSRE.table[n, model] <- reg$RMSRE
				mu.table[n, model] <- reg$mu
				regression.plot(reg, pch=16, main=paste0(model, ": prediction of ", biomass, " in ", dataset), 
									xlab="Observed values (g)", ylab="Predicted values (g)")
			}
		}
		ylim <- range(R2.table)
		ylim[1] <- ylim[1]-0.1
		ylim <- round(ylim, dig=1)
		barx <- hist.plot(R2.table, xaxt="n", main=paste("Evaluation based on", biomass, " in ", dataset), 
						ylim=ylim, col=model.cols[colnames(R2.table)], ylab=expression(R^2))
		axis(1, at=barx, labels=FALSE)
		text(x=barx, y=par("usr")[3], srt=45, adj=1, labels=colnames(R2.table), xpd=TRUE)
		
		ylim <- range(PCC.table)
		ylim[1] <- ylim[1]-0.1
		ylim <- round(ylim, dig=1)
		barx <- hist.plot(PCC.table, xaxt="n", main=paste("Evaluation based on", biomass, " in ", dataset), 
						ylim=ylim, col=model.cols[colnames(PCC.table)], ylab="PCC")
		axis(1, at=barx, labels=FALSE)
		text(x=barx, y=par("usr")[3], srt=45, adj=1, labels=colnames(PCC.table), xpd=TRUE)
		
		ylim <- range(RMSRE.table)
		ylim[1] <- 0
		barx <- hist.plot(RMSRE.table, xaxt="n", main=paste("Evaluation based on", biomass, " in ", dataset), 
						ylim=ylim, col=model.cols[colnames(RMSRE.table)], ylab="RMSRE")
		axis(1, at=barx, labels=FALSE)
		text(x=barx, y=par("usr")[3], srt=45, adj=1, labels=colnames(RMSRE.table), xpd=TRUE)
		
		ylim <- range(mu.table*1.05)
		barx <- hist.plot(mu.table, xaxt="n", main=paste("Evaluation based on", biomass, " in ", dataset), 
						ylim=ylim, col=model.cols[colnames(mu.table)], ylab=expression(mu))
		axis(1, at=barx, labels=FALSE)
		text(x=barx, y=par("usr")[3], srt=45, adj=1, labels=colnames(mu.table), xpd=TRUE)
		
		for(model in names(reg.models)){
			eval.table <- rbind(eval.table, data.frame(dataset=dataset, biomass=biomass, model=model, R2.avg=mean(R2.table[, model]), R2.std=std.error(R2.table[, model]), PCC.avg=mean(PCC.table[, model]), PCC.std=std.error(PCC.table[, model]), RMSRE.avg=mean(RMSRE.table[, model]), RMSRE.std=std.error(RMSRE.table[, model]), mu.avg=mean(mu.table[, model]), mu.std=std.error(mu.table[, model])))
		}
		par(op)
	}
}
dev.off()

write.table(eval.table, 'tables/Table 1 - evaluation1.data-eval.table.csv', sep=',', dec='.', na='', row.names=FALSE)

pdf('figures/Fig.1 - KN1121.model.evaluation1.pdf', width=10, heigh=10, pointsize=10)
op <- par(mfrow=c(3, 3), pty="m")
for(bio in c('FW', 'DW')){
	sub.data <- subset(eval.table, biomass==bio)
	plot.avg <- sub.data[, c('dataset', 'model', 'R2.avg')]
	plot.std <- sub.data[, c('dataset', 'model', 'R2.std')]
	plot.avg <- reshape(plot.avg, idvar="model", timevar="dataset", direction="wide")
	plot.std <- reshape(plot.std, idvar="model", timevar="dataset", direction="wide")
	colnames(plot.avg) <- gsub('R2.avg.', '', colnames(plot.avg))
	colnames(plot.std) <- gsub('R2.std.', '', colnames(plot.std))
	rownames(plot.avg) <- plot.avg$model
	rownames(plot.std) <- plot.std$model
	plot.avg <- as.matrix(plot.avg[,-1])
	plot.std <- as.matrix(plot.std[,-1])
	ylim <- range(plot.avg)
	ylim[1] <- ylim[1]-0.1
	ylim[2] <- ylim[2]+max(plot.std)+0.03
	ylim <- round(ylim, dig=1)
	if(ylim[2] > 1){ylim[2] <- 1}
	barx <- barplot(plot.avg, ylim=ylim, beside=TRUE, col=model.cols[rownames(plot.avg)], main=paste("Evaluation based on", bio), ylab=expression(R^2), xpd=FALSE, border=NA)
	text(barx, plot.avg+plot.std+abs(par('usr')[4]-par('usr')[3])*0.02, cex=0.8, labels=format(plot.avg, dig=3))
	error.bar(barx, plot.avg, plot.std, lwd=0.8, col=model.cols[rownames(plot.avg)])
	# my.legend("topleft", legend=names(model.cols), ncol=length(model.cols), col=model.cols, text.col=model.cols, pch=15, title="Models")
	cutoff <- apply(barx, 2, range)
	colnames(cutoff) <- colnames(plot.avg)
	for(ds in colnames(plot.avg)){
		lines(cutoff[, ds], rep(ref.r2[ds, bio], 2), lwd=1.5, lty=5)
	}
	
	plot.avg <- sub.data[, c('dataset', 'model', 'PCC.avg')]
	plot.std <- sub.data[, c('dataset', 'model', 'PCC.std')]
	plot.avg <- reshape(plot.avg, idvar="model", timevar="dataset", direction="wide")
	plot.std <- reshape(plot.std, idvar="model", timevar="dataset", direction="wide")
	colnames(plot.avg) <- gsub('PCC.avg.', '', colnames(plot.avg))
	colnames(plot.std) <- gsub('PCC.std.', '', colnames(plot.std))
	rownames(plot.avg) <- plot.avg$model
	rownames(plot.std) <- plot.std$model
	plot.avg <- as.matrix(plot.avg[,-1])
	plot.std <- as.matrix(plot.std[,-1])
	ylim <- range(plot.avg)
	ylim[1] <- ylim[1]-0.1
	ylim[2] <- ylim[2]+max(plot.std)+0.03
	ylim <- round(ylim, dig=1)
	if(ylim[2] > 1){ylim[2] <- 1}
	barx <- barplot(plot.avg, ylim=ylim, beside=TRUE, col=model.cols[rownames(plot.avg)], main=paste("Evaluation based on", bio), ylab="Pearson's correlation", xpd=FALSE, border=NA)
	text(barx, plot.avg+plot.std+abs(par('usr')[4]-par('usr')[3])*0.02, cex=0.8, labels=format(plot.avg, dig=3))
	error.bar(barx, plot.avg, plot.std, lwd=0.8, col=model.cols[rownames(plot.avg)])
	# my.legend("topleft", legend=names(model.cols), ncol=length(model.cols), col=model.cols, text.col=model.cols, pch=15, title="Models")
	cutoff <- apply(barx, 2, range)
	colnames(cutoff) <- colnames(plot.avg)
	for(ds in colnames(plot.avg)){
		lines(cutoff[, ds], rep(ref.cors[ds, bio], 2), lwd=1.5, lty=5)
	}
	
	plot.avg <- sub.data[, c('dataset', 'model', 'RMSRE.avg')]
	plot.std <- sub.data[, c('dataset', 'model', 'RMSRE.std')]
	plot.avg <- reshape(plot.avg, idvar="model", timevar="dataset", direction="wide")
	plot.std <- reshape(plot.std, idvar="model", timevar="dataset", direction="wide")
	colnames(plot.avg) <- gsub('RMSRE.avg.', '', colnames(plot.avg))
	colnames(plot.std) <- gsub('RMSRE.std.', '', colnames(plot.std))
	rownames(plot.avg) <- plot.avg$model
	rownames(plot.std) <- plot.std$model
	plot.avg <- as.matrix(plot.avg[,-1])
	plot.std <- as.matrix(plot.std[,-1])
	ylim <- range(plot.avg)
	ylim[2] <- max(ylim[2], max(ref.rmsre[, bio]), na.rm=TRUE) 
	ylim[1] <- 0
	ylim[2] <- ylim[2]+max(plot.std)*1.1
	barx <- barplot(plot.avg, ylim=ylim, beside=TRUE, col=model.cols[rownames(plot.avg)], main=paste("Evaluation based on", bio), ylab="RMSRE", xpd=FALSE, border=NA)
	text(barx, plot.avg+plot.std+abs(par('usr')[4]-par('usr')[3])*0.02, cex=0.8, labels=format(plot.avg, dig=2))
	error.bar(barx, plot.avg, plot.std, lwd=0.8, col=model.cols[rownames(plot.avg)])
	# my.legend("topleft", legend=names(model.cols), ncol=length(model.cols), col=model.cols, text.col=model.cols, pch=15, title="Models")
    cutoff <- apply(barx, 2, range)
	colnames(cutoff) <- colnames(plot.avg)
	for(ds in colnames(plot.avg)){
		lines(cutoff[, ds], rep(ref.rmsre[ds, bio], 2), lwd=1.5, lty=5)
	}
}
plot(1, type="n", xaxt="n", yaxt="n", frame=FALSE, xlab="", ylab="")
my.legend("topleft", legend=names(model.cols), ncol=length(model.cols), col=model.cols, text.col=model.cols, pch=15, title="Models")
par(op)
dev.off()

########## 
traits.class <- unlist(lapply(barley.traits, trait.class.map))
names(traits.class) <- barley.traits
color.map <- c(Geometric="#3778AD", Color="#006934", FLUO="#F8B62D", NIR="#727171")
model <- 'rf.regression'
category.effect <- data.frame(stringsAsFactors=FALSE)
for(biomass in c('FW', 'DW')){
	for(dataset in names(barley.list)){
		input <- barley.list[[dataset]]
		for(clazz in names(color.map)){
			sub.traits <- barley.traits[which(traits.class==clazz)]
			reg.data <- input[, c(biomass, sub.traits)]
			reg <- regression.cross.validation(reg.data, model=model)
			cat(dataset, "\t", biomass, "\t", clazz, "\t", reg$R2, "\n")
			category.effect <- rbind(category.effect, data.frame(biomass=biomass, dataset=dataset, clazz=clazz, R2=reg$R2))
		}
	}
}

trait.prediction <- data.frame(stringsAsFactors=FALSE)
for(biomass in c('FW', 'DW')){
	for(dataset in names(barley.list)){
		input <- barley.list[[dataset]]
		for(trait in barley.traits){
			reg.data <- input[, c(biomass, trait)]
			reg <- regression.cross.validation(reg.data, model=model)
			cat(dataset, "\t", biomass, "\t", trait, "\t", reg$R2, "\n")
			trait.prediction <- rbind(trait.prediction, data.frame(biomass=biomass, dataset=dataset, trait=trait, R2=reg$R2))
		}
	}
}

trait.relimp <- data.frame(trait=barley.traits, stringsAsFactors=FALSE)
rownames(trait.relimp) <- barley.traits
for(biomass in c('FW', 'DW')){
	for(dataset in names(barley.list)){
		input <- barley.list[[dataset]]
		reg.data <- input[, c(biomass, barley.traits)]
		reg <- regression.cross.validation(reg.data, model=model)
		old.names <- colnames(trait.relimp)
		trait.relimp <- cbind(trait.relimp, reg$relimp[rownames(trait.relimp)])
		colnames(trait.relimp) <- c(old.names, paste(dataset, biomass, sep="-"))
	}
}

write.table(category.effect, 'tables/Table 1 - evaluation2.data-category.effect.csv', sep=',', dec='.', na='', row.names=FALSE)
write.table(trait.prediction, 'tables/Table 1 - evaluation2.data-trait.prediction.csv', sep=',', dec='.', na='', row.names=FALSE)
write.table(trait.relimp, 'tables/Table 1 - evaluation2.data-trait.relimp.csv', sep=',', dec='.', na='', row.names=FALSE)

### test the predictive power of models 
pdf('figures/Fig.1 - KN1121.model.evaluation2.pdf', width=10, heigh=10, pointsize=10)
op <- par(mfrow=c(2, 2), pty="m")
RF.prid.power <- subset(eval.table, model=="RF")[, c("dataset", "biomass", "R2.avg", "R2.std")]
RF.prid.power[, "R2.avg"] <- RF.prid.power[, "R2.avg"] + RF.prid.power[, "R2.std"]
for(bio in c('FW', 'DW')){
	sub.data <- subset(RF.prid.power, biomass==bio)[, c("dataset", "biomass", "R2.avg")]
	sub.data[, 2] <- "All"
	colnames(sub.data) <- c("dataset", "clazz", "R2")
	plot.data <- subset(category.effect, biomass==bio)[, c("dataset", "clazz", "R2")]
	plot.data <- rbind(sub.data, plot.data)
	plot.data <- reshape(plot.data, idvar="clazz", timevar="dataset", direction="wide")
	colnames(plot.data) <- gsub('R2.', '', colnames(plot.data))
	rownames(plot.data) <- plot.data$clazz
	pcols <- color.map[as.character(plot.data$clazz)]
	pcols[is.na(pcols)] <- "#E61673"
	barplot(as.matrix(plot.data[,-1]), ylim=c(0, 1), beside=TRUE, col=pcols, main=paste("Evaluation based on", bio), ylab=expression(R^2), border=NA)
	# my.legend("topleft", legend=names(color.map), ncol=length(color.map), col=color.map, text.col=color.map, pch=15, title="Traits")
}
par(op)

op <- par(mfrow=c(4, 1), pty="m")
for(ds in names(barley.list)){
	for(bio in c('FW', 'DW')){
		R2.data <- subset(trait.prediction, biomass==bio & dataset==ds)[, c("trait", "R2")]
		rownames(R2.data) <- R2.data$trait
		relimp.data <- trait.relimp[, c('trait', paste(ds, bio, sep="-"))]
		
		plot.data <- data.frame(R2=R2.data[barley.traits, 2], relimp=relimp.data[barley.traits, 2])
		rownames(plot.data) <- barley.traits
		plot.data <- plot.data[with(plot.data, order(-relimp)), ]
		barx <- barplot(plot.data[, 'relimp'], border=NA, col=color.map[traits.class[rownames(plot.data)]], ylab="RI (%IncMSE)")
		m <- median(plot.data[, 'relimp'])
		## abline(h=m, col='red', lty=5)
		lines(range(barx), c(m,m), col='red', lty=5)
		title(main=paste0("Relative importance of each trait\n",ds, ", ", bio), xlab='Traits')
		my.legend("topright", legend=names(color.map), col=color.map, text.col=color.map, pch=15, title="Traits")
		
		barx <- barplot(plot.data[, 'R2'], border=NA, col=color.map[traits.class[rownames(plot.data)]], ylab=expression(R^2))
		m <- median(plot.data[, 'R2'])
		## abline(h=m, col='red', lty=5)
		lines(range(barx), c(m,m), col='red', lty=5)
		title(main=paste0("Predictive power of each individual trait\n", ds, ", ", bio), xlab='Traits')
		axis(1, barx, trait.index[rownames(plot.data)], las=1, tick=FALSE, padj=-1)
		my.legend("topright", legend=names(color.map), col=color.map, text.col=color.map, pch=15, title="Traits")
	}
}
par(op)
op <- par(mfrow=c(3, 3), pty="m")
for(ds in names(barley.list)){
	plot.data <- trait.relimp[barley.traits, grep(gsub("\\+", ".", paste0("^", ds, "-")), colnames(trait.relimp))]
	plot(plot.data, col=color.map[traits.class], pch=16, main='Relative importance')
	
	fit <- lm(plot.data[,2]~plot.data[,1]-1)
	abline(fit, col='red', lwd=2)
	res <- residuals(fit)
	names(res) <- rownames(plot.data)
	qtl <- quantile(res,  probs = c(5, 95)/100)
	outlier.id <- names(which(res < qtl[1] | res > qtl[2]))
	points(plot.data[outlier.id, ], pch=16, cex=2, col=color.map[traits.class[outlier.id]])
	text(plot.data[outlier.id, ], outlier.id, pos=1)
	
	p.cor <- format(cor(plot.data[, 1], plot.data[, 2]), dig=4)
	lgd <- substitute(italic(r) == RR, list(RR=p.cor))
	my.legend('topleft', legend=lgd, pch=NA)
}
dev.off()


view.class <- unlist(lapply(barley.traits, trait.view.map))
names(view.class) <- barley.traits
view.map <- c(Side="#4daf4a", Top="#377eb8", Both="#ff7f00")
model <- 'rf.regression'
view.effect <- data.frame(stringsAsFactors=FALSE)
for(biomass in c('FW', 'DW')){
	for(dataset in names(barley.list)){
		input <- barley.list[[dataset]]
		for(clazz in names(view.map)){
			sub.traits <- barley.traits[which(view.class==clazz)]
			reg.data <- input[, c(biomass, sub.traits)]
			reg <- regression.cross.validation(reg.data, model=model)
			cat(dataset, "\t", biomass, "\t", clazz, "\t", reg$R2, "\n")
			view.effect <- rbind(view.effect, data.frame(biomass=biomass, dataset=dataset, clazz=clazz, R2=reg$R2))
		}
	}
}
write.table(view.effect, 'tables/Table 1 - evaluation3.data-view.effect.csv', sep=',', dec='.', na='', row.names=FALSE)


### test the predictive power of models 
pdf('figures/Fig.1 - KN1121.model.evaluation3.pdf', width=10, heigh=10, pointsize=10)
op <- par(mfrow=c(2, 2), pty="m")
RF.prid.power <- subset(eval.table, model=="RF")[, c("dataset", "biomass", "R2.avg", "R2.std")]
RF.prid.power[, "R2.avg"] <- RF.prid.power[, "R2.avg"] + RF.prid.power[, "R2.std"]
for(bio in c('FW', 'DW')){
	sub.data <- subset(RF.prid.power, biomass==bio)[, c("dataset", "biomass", "R2.avg")]
	sub.data[, 2] <- "All"
	colnames(sub.data) <- c("dataset", "clazz", "R2")
	plot.data <- subset(view.effect, biomass==bio)[, c("dataset", "clazz", "R2")]
	plot.data <- rbind(sub.data, plot.data)
	plot.data <- reshape(plot.data, idvar="clazz", timevar="dataset", direction="wide")
	colnames(plot.data) <- gsub('R2.', '', colnames(plot.data))
	rownames(plot.data) <- plot.data$clazz
	pcols <- view.map[as.character(plot.data$clazz)]
	pcols[is.na(pcols)] <- "#525252"
	barplot(as.matrix(plot.data[,-1]), ylim=c(0, 1), beside=TRUE, col=pcols, main=paste("Evaluation based on", bio), ylab=expression(R^2), border=NA)
	my.legend("topleft", legend=names(view.map), ncol=length(view.map), col=view.map, text.col=view.map, pch=15, title="Camera view")
}
par(op)
dev.off()

##########################
## full model
model <- 'rf.regression' ## the 'best' model
all.traits <- colnames(KN1121)[-c(1:5)]
trait.relimp2 <- data.frame(trait=all.traits, stringsAsFactors=FALSE)
rownames(trait.relimp2) <- all.traits
for(biomass in c('FW', 'DW')){
	for(dataset in names(barley.list)){
		if(length(grep('\\+', dataset)) > 0){
			reg.data <- KN1121[, c(biomass, all.traits)]
		}else{
			reg.data <- subset(KN1121, Treatment==tolower(dataset))[, c(biomass, all.traits)]
		}
		reg <- regression.cross.validation(reg.data, model=model)
		cat(dataset, "\t", biomass, "\t", reg$R2, "\n")
		old.names <- colnames(trait.relimp2)
		trait.relimp2 <- cbind(trait.relimp2, reg$relimp[rownames(trait.relimp2)])
		colnames(trait.relimp2) <- c(old.names, paste(dataset, biomass, sep="-"))
	}
}

write.table(trait.relimp2, 'tables/Table 1 - evaluation2.data-trait.relimp2.csv', sep=',', dec='.', na='', row.names=FALSE)


################ part 2 ################
rm(list=ls())
options(stringsAsFactors=FALSE)

library(plotrix)
library(pheatmap)
library(gplots)
require(pcaMethods)
require(ade4)
library(RColorBrewer)

if(file.exists('my.code.R')){
	source('my.code.R')
}
source('models.R')
source('pca.plot.R')

meta.data <- read.csv("data/FW-1121_1130_1137.csv", row.names=1, header=TRUE, sep=',', stringsAsFactors=FALSE)

KN1121 <- read.csv("data/1121KN - regression.data.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)
KN1130 <- read.csv("data/1130KN - regression.data.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)
KN1137 <- read.csv("data/1137KN - regression.data.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)

plant.data <- rbind(cbind(Dataset="KN1121", KN1121), cbind(Dataset="KN1130", KN1130), cbind(Dataset="KN1137", KN1137))
rownames(plant.data) <- plant.data$PlantID
treatments <- sort(unique(plant.data$Treatment))
genotypes <- unique(plant.data$Genotype)
datasets <- unique(plant.data$Dataset)

meta.data$Dataset <- plant.data[rownames(meta.data), "Dataset"]

n.genotypes <- length(genotypes)
geno.col <- rainbow(n.genotypes)
names(geno.col) <- genotypes
geno.lty <- c(1:n.genotypes)
names(geno.lty) <- genotypes
geno.pch <- 0:n.genotypes %% 25
names(geno.pch) <- genotypes

treat.col <- c("#1B9E77", "#A6761D")
names(treat.col) <- treatments
data.col <- c("#377EB8", "#FF7F00", "#984EA3")
names(data.col) <- datasets

## normalization 
plant.scaled.data <- plant.data
plant.scaled.data[, -c(1:5)] <- apply(plant.scaled.data[, -c(1:5)], 2, function(x){x/max(abs(x))})

treatment.pch <- as.integer(as.factor(treatments))
names(treatment.pch) <- treatments
lgd.txt <- expand.grid(D=names(data.col), T=names(treat.col), stringsAsFactors=FALSE)
rownames(lgd.txt) <- paste(lgd.txt$D, lgd.txt$T, sep=', ')
lgd.col <- data.col[lgd.txt$D]
lgd.pch <- treatment.pch[lgd.txt$T]
lgd.lbs <- rownames(lgd.txt)

biomass <- 'FW'
traits <- colnames(plant.data)[-c(1:5)]
cor.matrix <- data.frame(Trait=traits, stringsAsFactors=FALSE)
rownames(cor.matrix) <- traits
for(dataset in datasets){
	sub.data <- subset(plant.data, Dataset==dataset)
	for(trait in traits){
		cor.matrix[trait, dataset] <- cor(sub.data[, biomass], sub.data[, trait])
	}
}
cor.matrix <- cor.matrix[with(cor.matrix, order(KN1121+KN1130+KN1137)), ]


times <- 10
biomass <- 'FW'
model.cols <- c(MLR="#f39800", MARS="#00a0e9", RF="#e4007f", SVR="#1b9e77")
traits <- colnames(plant.scaled.data)[-c(1:5)]
pdf('figures/Fig.2 - whole.evaluation.pdf', width=10, heigh=10, pointsize=10)
op <- par(mfrow=c(4, 4), pty="m")
reg.data <- plant.scaled.data[, c(biomass, traits)]
## reg.data[,-1] <- apply(reg.data[,-1], 2, function(x){x/max(abs(x))})
R2.table <- data.frame(matrix(0, ncol=length(reg.models), nrow=times))
PCC.table <- data.frame(matrix(0, ncol=length(reg.models), nrow=times))
RMSRE.table <- data.frame(matrix(0, ncol=length(reg.models), nrow=times))
mu.table <- data.frame(matrix(0, ncol=length(reg.models), nrow=times))
colnames(R2.table) <- colnames(PCC.table) <- colnames(RMSRE.table) <- colnames(mu.table) <- names(reg.models)
for(n in 1:times){
	cat("Biomass: ", biomass, "\t ", n, "\n")
	for(model in names(reg.models)){
		reg <- regression.cross.validation(reg.data, model=reg.models[[model]])
		R2.table[n, model] <- reg$R2
		PCC.table[n, model] <- reg$PCC
		RMSRE.table[n, model] <- reg$RMSRE
		mu.table[n, model] <- reg$mu
		regression.plot(reg, pch=16, main=paste0(model, ": prediction of ", biomass), 
							xlab="Observed values (g)", ylab="Predicted values (g)")
	}
}
par(op)
op <- par(mfrow=c(2, 2), pty="m")
ylim <- range(R2.table)
ylim[1] <- ylim[1]-0.1
ylim <- round(ylim, dig=1)
barx <- hist.plot(R2.table, xaxt="n", main=paste("Evaluation based on", biomass), 
				ylim=ylim, col=model.cols[colnames(R2.table)], ylab=expression(R^2))
axis(1, at=barx, labels=FALSE)
text(x=barx, y=par("usr")[3], srt=45, adj=1, labels=colnames(R2.table), xpd=TRUE)

ylim <- range(PCC.table)
ylim[1] <- ylim[1]-0.1
ylim <- round(ylim, dig=1)
barx <- hist.plot(PCC.table, xaxt="n", main=paste("Evaluation based on", biomass), 
				ylim=ylim, col=model.cols[colnames(PCC.table)], ylab="PCC")
axis(1, at=barx, labels=FALSE)
text(x=barx, y=par("usr")[3], srt=45, adj=1, labels=colnames(PCC.table), xpd=TRUE)

ylim <- range(RMSRE.table)
ylim[1] <- 0
barx <- hist.plot(RMSRE.table, xaxt="n", main=paste("Evaluation based on", biomass), 
				ylim=ylim, col=model.cols[colnames(RMSRE.table)], ylab="RMSRE")
axis(1, at=barx, labels=FALSE)
text(x=barx, y=par("usr")[3], srt=45, adj=1, labels=colnames(RMSRE.table), xpd=TRUE)

ylim <- range(mu.table*1.05)
barx <- hist.plot(mu.table, xaxt="n", main=paste("Evaluation based on", biomass), 
				ylim=ylim, col=model.cols[colnames(mu.table)], ylab=expression(mu))
axis(1, at=barx, labels=FALSE)
text(x=barx, y=par("usr")[3], srt=45, adj=1, labels=colnames(mu.table), xpd=TRUE)

par(op)
dev.off()


### Cluster based on FW data
pdf("figures/Fig.2 - dataset.summary.pdf", width=10, heigh=10, pointsize=10)
op <- par(mfrow=c(2, 2), pty="m")
FW.genotype.mean <- aggregate(plant.scaled.data$FW, plant.scaled.data[, c('Dataset', 'Treatment', 'Genotype')], mean)
dataset.summary <- reshape(FW.genotype.mean, v.names="x", idvar=c("Dataset", "Treatment"), timevar="Genotype", direction="wide")
rownames(dataset.summary) <- paste(dataset.summary$Dataset, dataset.summary$Treatment, sep=", ")
colnames(dataset.summary) <- gsub('x.', '', colnames(dataset.summary))
## dataset.summary[, -c(1:3)] <- t(apply(dataset.summary[, -c(1:3)], 1, function(x){x/max(abs(x))})) ## scale
plot(hclust(dist(dataset.summary[, -c(1:3)])), col="#487AA1", col.main="#45ADA8", col.lab="#7C8071", col.axis="#F38630", xlab=NA)
plot(plant.scaled.data[, c('FW', 'volume.fluo.iap')], col=data.col[plant.scaled.data$Dataset], pch=16)
boxplot(FW ~ Dataset+Treatment, data=plant.scaled.data[, c('FW', 'Treatment', 'Dataset')], col=data.col, ylab='Fresh weight (g)')

aov(FW ~ Dataset+Treatment, data=plant.scaled.data[, c('FW', 'Treatment', 'Dataset')])

trait.col <- unlist(lapply(rownames(cor.matrix), trait.color.map))
cor.median <- apply(cor.matrix[, datasets], 1, mean)
barplot(cor.median, col=trait.col, border=NA, xlab=NA, ylim=c(-0.5,1), ylab="PCC", panel.first=grid(), space=0.3)

barx <- barplot(as.matrix(t(cor.matrix[, datasets])), border=NA, col=data.col, xlab=NA, ylab="PCC", panel.first=grid(), space=0.3)
points(barx, t(apply(cor.matrix[, datasets], 1, mean)), col="white", pch=16)

par(op)

pairs(t(dataset.summary[, -c(1:3)]), lower.panel=panel.smooth, upper.panel=panel.cor, main='FW')

### PCA based on image data
pca.dudi <- dudi.pca(plant.scaled.data[, -c(1:5)], center=TRUE, scale=TRUE, scan=FALSE, nf=4)
R2 <- pca.dudi$eig/sum(pca.dudi$eig) * 100
loading <- pca.dudi$co
plant.pchs <- as.integer(as.factor(plant.data$Treatment))
plant.cols <- data.col[plant.data$Dataset]
trait.cols <- unlist(lapply(rownames(pca.dudi$co), trait.color.map))
my.pca.plot(pca.dudi, sp=plant.pchs, sc=plant.cols, lc=trait.cols, cex=1, lwd=0.8)
my.legend("topright", legend=names(data.col), pch=15, col=data.col, text.col=data.col, title="Experiment")
my.legend("right", legend=names(treat.col), col=treat.col, text.col=treat.col, pch=as.integer(as.factor(names(treat.col))), title="Treatment")

col.map <- colorRampPalette(brewer.pal(n=9, name="Greys")[1:7])(100) ## bluered(100)

dataset.cors <- cor(t(dataset.summary[, -c(1:3)]), method='p')
heatmap.2(dataset.cors, col=col.map, scale='none', key=TRUE, symkey=FALSE, density.info='none', trace='none', margins=c(15,15))

#### HCA based on image data
genotype.data <- aggregate(plant.data[, -c(1:5)], plant.data[, c('Dataset', 'Treatment', 'Genotype')], mean)
rownames(genotype.data) <- paste(genotype.data[,1], genotype.data[,2], genotype.data[,3], sep="-")
data.norm <- t(apply(genotype.data[, -c(1:3)], 1, function(x){x/max(x, na.rm=TRUE)}))
dataset.cors <- cor(t(data.norm), method='p')

hv1 <- heatmap.2(dataset.cors, col=col.map, scale='none', key=TRUE, symkey=FALSE, density.info='none', trace='none', RowSideColors=data.col[gsub('-.*','', rownames(dataset.cors))], ColSideColors=treat.col[gsub('-.*', '', gsub('^.*?-','', colnames(dataset.cors)))], labRow=NA, labCol=NA) ## geno.col[meta.data[rownames(plant.cors), 'Genotype']]

dataset.cors.new <- dataset.cors[hv1$colInd, hv1$colInd]
dataset.cors.new[lower.tri(dataset.cors.new)] <- NA
image(1:nrow(dataset.cors.new),1:ncol(dataset.cors.new), dataset.cors.new, col=col.map, xlab="Plants", ylab="Plants")

plot(1, type="n", xaxt="n", yaxt="n", frame=FALSE, xlab="", ylab="")
my.legend("left", legend=names(data.col), pch=15, col=data.col, text.col=data.col, title="Experiment")
my.legend("right", legend=names(treat.col), col=treat.col, text.col=treat.col, pch=as.integer(as.factor(names(treat.col))), title="Treatment")
my.legend("center", legend=lgd.lbs, pch=lgd.pch, col=lgd.col, text.col=lgd.col, title="Experiment")
dev.off()


outlier.detection.by.fit <- function(x=x, y=y, alpha=0.05, scale=1){
	if(alpha > 1 | alpha < 0){alpha <- 0.05}
	if(alpha > 0.5){alpha <- 1 - alpha}
	res <- residuals(lm(y~x))
	qtl <- quantile(res,  probs=c(alpha, 1-alpha)) * scale
	oids <- which(res < qtl[1] | res > qtl[2])
	x[oids] <- NA
	y[oids] <- NA
	return(list(x=x, y=y))
}

## cohort prediction
pdf("figures/Fig.2 - cohort.prediction.log.pdf", width=10, heigh=10, pointsize=10)
op <- par(mfrow=c(2, 2), pty="m")
model <- 'rf.regression' ## the 'best' model
biomass <- 'FW'
traits <- colnames(plant.scaled.data)[-c(1:5)]
reg.table <- list()
times <- 10
for(train.exp in datasets){
	for(train.treat in treatments){
		train.id <- paste(train.treat, train.exp, sep=", ")
		train.data <- subset(plant.scaled.data, Dataset==train.exp & Treatment==train.treat)[, c(biomass, traits)]
		for(test.exp in datasets){
			for(test.treat in treatments){
				test.id <- paste(test.treat, test.exp, sep=", ")
				test.data <- subset(plant.scaled.data, Dataset==test.exp & Treatment==test.treat)[, c(biomass, traits)]
				cat(train.id, " --> ", test.id, "\n")
				if(train.id != test.id){
					predicted <- c()
					for(k in 1:times){
						reg <- regression.cohort.validation(train.data, test.data, model=model)
						predicted <- cbind(predicted, reg$predicted)
					}
					predicted <- rowMeans(predicted)
				}else{
					reg <- regression.cross.validation(test.data, model=model, nfold=times)
					predicted <- reg$predicted
				}
				observed <- test.data[names(predicted), biomass]
				reg.table[[train.id]][[test.id]] <- cbind(observed=observed, predicted=predicted)
				plot.data <- evaluation.criteria(observed, predicted)
				plot.data$observed <- observed
				plot.data$predicted <- predicted
				tit <- vector('expression',1)
				tit[1] <- substitute(expression(TR %->% TE), list(TR=train.id, TE=test.id))[2]
				regression.plot(plot.data, mu=TRUE, xlab="Observed values (g)", ylab="Predicted values (g)", main=tit[1])
			}
		}
	}
}
par(op)
dev.off()

pred.power <- data.frame(stringsAsFactors=FALSE)
pdf("figures/Fig.2 - cohort.prediction1.pdf", width=10, heigh=10, pointsize=10)
nr <- nc <- length(reg.table)+1
op <- par(mfrow=c(nr, nc), mar=rep(0.2, 4), oma=rep(2, 4)) # , mgp=c(0, -1.4, 0)
plot(0, 0, type='n', xlab=NA, ylab=NA, axes=FALSE, frame.plot=FALSE)
text(0, 0, labels="Train/Test", cex=2, col='red')
for(test.id in sort(names(reg.table))){
	plot(0, 0, type='n', xlab=NA, ylab=NA, axes=FALSE, frame.plot=FALSE)
	text(0, 0, labels=gsub(", ", "\n", test.id), cex=2)
}
rid <- 1
for(train.id in sort(names(reg.table))){
	cid <- 1
	for(test.id in sort(names(reg.table))){
		## plot(observed, predicted, col=col, xlab=NA, ylab=NA, axes=FALSE, frame.plot=TRUE)
		if(cid == 1){
			plot(0, 0, type='n', xlab=NA, ylab=NA, axes=FALSE, frame.plot=FALSE)
			text(0, 0, labels=gsub(", ", "\n", train.id), cex=2)
		}
		reg.result <- reg.table[[test.id]][[train.id]]
		observed <- reg.result[, 'observed']
		predicted <- reg.result[, 'predicted']
		col <- colorRampPalette(c('#4292c6','#08306b'))(length(predicted))
		col <- col[rank(predicted)]
		without.outliers <- outlier.detection.by.fit(observed, predicted)
		plot(without.outliers$x, without.outliers$y, col=col, xlab=NA, ylab=NA, axes=FALSE, frame.plot=TRUE)
		abline(lm(predicted ~ observed), lwd=1, lty=1, col='grey')
		abline(lm(predicted ~ observed-1), lwd=1, lty=4, col='red')
		eval.result <- evaluation.criteria(observed, predicted)
		o.names <- rownames(pred.power)
		pred.power <- rbind(pred.power, unlist(eval.result))
		colnames(pred.power) <- names(unlist(eval.result))
		rownames(pred.power) <- c(o.names, paste(train.id, "-", test.id))
		lgds <- vector('expression',4)
		lgds[1] <- substitute(expression(italic(r)==PCC), list(PCC=format(eval.result$PCC, dig=4)))[2]
		lgds[2] <- substitute(expression(italic(R)^2==R2), list(R2=format(eval.result$R2, dig=4)))[2]
		lgds[3] <- substitute(expression(paste("RMSRE = ", RMSRE, sep="")), list(RMSRE=format(eval.result$RMSRE, dig=3, scienticif=TRUE)))[2]
		lgds[4] <- substitute(expression(mu==MU), list(MU=format(eval.result$mu, dig=4)))[2]
		legend("topleft", inset=c(-0.05, 0), legend=lgds, box.col=NA, bg=NA, pch=NA)
		if(cid == 1){
		#	axis(2)
		}
		if(rid == nr-1){
			axis(1)
		}
		cid <- 1 + cid
	}
	rid <- 1 + rid
}
par(op)
dev.off()

pdf("figures/Fig.2 - cohort.prediction.log2.pdf", width=10, heigh=10, pointsize=10)
op <- par(mfrow=c(2, 2), pty="m")
model <- 'rf.regression' ## the 'best' model
biomass <- 'FW'
traits <- colnames(plant.scaled.data)[-c(1:5)]
cohort <- list()
times <- 10
for(exp1 in datasets){
	train.data <- subset(plant.scaled.data, Dataset==exp1)[, c(biomass, traits)]
	for(exp2 in datasets){
		test.data <- subset(plant.scaled.data, Dataset==exp2)[, c(biomass, traits)]
		reg <- regression.cohort.validation(train.data, test.data, model=model)
		cat(exp1, " --> ", exp2, "\n")
		if(exp1 != exp2){
			predicted <- c()
			for(k in 1:times){
				reg <- regression.cohort.validation(train.data, test.data, model=model)
				predicted <- cbind(predicted, reg$predicted)
			}
			predicted <- rowMeans(predicted)
		}else{
			reg <- regression.cross.validation(test.data, model=model, nfold=times)
			predicted <- reg$predicted
		}
		observed <- test.data[names(predicted), biomass]
		cohort[[exp1]][[exp2]] <- cbind(observed=observed, predicted=predicted)
		plot.data <- evaluation.criteria(observed, predicted)
		plot.data$observed <- observed
		plot.data$predicted <- predicted
		tit <- vector('expression',1)
		tit[1] <- substitute(expression(TR %->% TE), list(TR=exp1, TE=exp2))[2]
		regression.plot(plot.data, mu=TRUE, xlab="Observed values (g)", ylab="Predicted values (g)", main=tit[1])
	}
}
par(op)
dev.off()


pdf("figures/Fig.2 - cohort.prediction2.pdf", width=10, heigh=10, pointsize=10)
nr <- nc <- length(cohort)
op <- par(mfrow=c(nr, nc), mar=rep(0.2, 4), oma=rep(15, 4)) # , mgp=c(0, -1.4, 0)
rid <- 1
for(exp1 in sort(names(cohort))){
	cid <- 1
	for(exp2 in sort(names(cohort))){
		## plot(observed, predicted, col=col, xlab=NA, ylab=NA, axes=FALSE, frame.plot=TRUE)
		reg.result <- cohort[[exp2]][[exp1]]
		observed <- reg.result[, 'observed']
		predicted <- reg.result[, 'predicted']
		col <- colorRampPalette(c('#4292c6','#08306b'))(length(predicted))
		col <- col[rank(predicted)]
		without.outliers <- outlier.detection.by.fit(observed, predicted)
		plot(without.outliers$x, without.outliers$y, col=treat.col[meta.data[names(predicted), 'Treatment']], pch=treatment.pch[meta.data[names(predicted), 'Treatment']], xlab=NA, ylab=NA, axes=FALSE, frame.plot=TRUE)
		abline(lm(predicted ~ observed), lwd=1, lty=1, col='grey')
		abline(lm(predicted ~ observed-1), lwd=1, lty=4, col='red')
		eval.result <- evaluation.criteria(observed, predicted)
		o.names <- rownames(pred.power)
		pred.power <- rbind(pred.power, unlist(eval.result))
		colnames(pred.power) <- names(unlist(eval.result))
		rownames(pred.power) <- c(o.names, paste(exp1, "-", exp2))
		lgds <- vector('expression',4)
		lgds[1] <- substitute(expression(italic(r)==PCC), list(PCC=format(eval.result$PCC, dig=4)))[2]
		lgds[2] <- substitute(expression(italic(R)^2==R2), list(R2=format(eval.result$R2, dig=4)))[2]
		lgds[3] <- substitute(expression(paste("RMSRE = ", RMSRE, sep="")), list(RMSRE=format(eval.result$RMSRE, dig=3, scienticif=TRUE)))[2]
		lgds[4] <- substitute(expression(mu==MU), list(MU=format(eval.result$mu, dig=4)))[2]
		legend("bottomright", inset=c(-0.05, 0), legend=lgds, box.col=NA, bg=NA, pch=NA)
		if(cid == 1){
			axis(2)
		}
		if(rid == nr){
			axis(1)
		}
		cid <- 1 + cid
	}
	rid <- 1 + rid
}
par(op)
dev.off()

groups <- rownames(pred.power)
names(groups) <- rownames(pred.power)
for(g in rownames(pred.power)){
	gs <- unlist(strsplit(g, " - "))
	if(length(grep(', ', g)) < 1){
		if(gs[1]==gs[2]){
			groups[g] <- 'Within, all plants'
		}else{
			groups[g] <- 'Cross, all plants'
		}
	}else{
		if(gs[1]==gs[2]){
			if(length(grep('control', gs[1])) > 0){
				groups[g] <- 'Within, control'
			}else{
				groups[g] <- 'Within, stress'
			}
		}else{
			d1 <- unlist(strsplit(gs[1], ", "))
			d2 <- unlist(strsplit(gs[2], ", "))
			if(d1[1]==d2[1] & d1[1]=='control'){
				groups[g] <- 'Cross, control'
			}else if(d1[1]==d2[1] & d1[1]=='stress'){
				groups[g] <- 'Cross, stress'
			}else if(d1[1]=='control' & d2[1]=='stress'){
				groups[g] <- 'Control-stress'
			}else if(d1[1]=='stress' & d2[1]=='control'){
				groups[g] <- 'Stress-control'
			}
		}
	}
}

pred.power <- cbind(Group=groups, pred.power)
n.groups <- table(pred.power$Group)

pdf("figures/Fig.2 - cohort.prediction3.pdf", width=10, heigh=10, pointsize=10)
op <- par(mfrow=c(3, 3), pty="m")
bpt <- boxplot(R2~Group, data=pred.power, col='grey', main="R2", ylab="R2", xaxt="n", xlab="", ylim=c(0,1))
axis(1, at=seq_along(n.groups), labels=FALSE)
text(x=seq_along(n.groups), y=par("usr")[3]-0.03, srt=45, adj=1, labels=names(n.groups), xpd=TRUE)
text(1:length(n.groups), bpt$stats[5, ], paste0('n=', n.groups), pos=3)

bpt <- boxplot(PCC~Group, data=pred.power, col='grey', main="PCC", ylab="PCC", xaxt="n", xlab="", ylim=c(0,1))
axis(1, at=seq_along(n.groups), labels=FALSE)
text(x=seq_along(n.groups), y=par("usr")[3]-0.03, srt=45, adj=1, labels=names(n.groups), xpd=TRUE)
text(1:length(n.groups), bpt$stats[5, ], paste0('n=', n.groups), pos=3)

bpt <- boxplot(RMSRE~Group, data=pred.power, col='grey', main="RMSRE", ylab="RMSRE", xaxt="n", xlab="", ylim=c(0,1))
axis(1, at=seq_along(n.groups), labels=FALSE)
text(x=seq_along(n.groups), y=par("usr")[3]-0.03, srt=45, adj=1, labels=names(n.groups), xpd=TRUE)
text(1:length(n.groups), bpt$stats[5, ], paste0('n=', n.groups), pos=3)

bpt <- boxplot(mu~Group, data=pred.power, col='grey', main="mu", ylab="mu", xaxt="n", xlab="")
axis(1, at=seq_along(n.groups), labels=FALSE)
text(x=seq_along(n.groups), y=par("usr")[3]-0.03, srt=45, adj=1, labels=names(n.groups), xpd=TRUE)
text(1:length(n.groups), bpt$stats[5, ], paste0('n=', n.groups), pos=3)

par(op)
dev.off()


#################### relative importance ###############################
trait.relimp <- data.frame(stringsAsFactors=FALSE)
traits <- colnames(plant.scaled.data)[-c(1:5)]
trait.index <- 1:length(traits)
names(trait.index) <- sort(traits)
biomass <- 'FW'
model <- 'rf.regression' ## the 'best' model
for(dataset in datasets){
	for(treatment in treatments){
		reg.data <- subset(plant.scaled.data, Dataset==dataset & Treatment==treatment)[, c(biomass, traits)]
		reg <- regression.cross.validation(reg.data, model=model)
		trait.relimp <- rbind(trait.relimp, data.frame(Trait=traits, Biomass=biomass, Dataset=dataset, Treatment=treatment, Imp=reg$relimp[traits]))
	}
}

#### for DW, only two datasets are applied
two.data <- read.csv("data/Regression-1121_1137.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)
rownames(two.data) <- two.data$PlantID
two.scaled.data <- two.data
two.scaled.data[, -c(1:6)] <- apply(two.scaled.data[, -c(1:6)], 2, function(x){x/max(abs(x))})

biomass <- 'DW'
model <- 'rf.regression' ## the 'best' model
for(dataset in datasets){
	for(treatment in treatments){
		reg.data <- subset(two.scaled.data, Dataset==dataset & Treatment==treatment)[, c(biomass, traits)]
		if(nrow(reg.data) > 0){
			reg <- regression.cross.validation(reg.data, model=model)
			trait.relimp <- rbind(trait.relimp, data.frame(Trait=traits, Biomass=biomass, Dataset=dataset, Treatment=treatment, Imp=reg$relimp[traits]))
		}
	}
}

## for the function of pairs 
panel.p.cor <- function(x, y, digits=4, prefix="", ...){
	usr <- par("usr"); on.exit(par(usr))
	r <- cor(x, y, use="pairwise.complete.obs")
	par(usr=c(0, 1, 0, 1))
	txt <- ifelse(r >= 0, 
				sprintf(paste0("%.", digits, "f"), r), ## format(c(r, 0.123456789), digits=digits)[1]
				txt <- sprintf(paste0("%.", digits-1, "f"), r)
			)
	txt <- paste(prefix, txt, sep="")
	cex.cor <- 0.8/strwidth(txt)
	text(0.5, 0.5, txt, cex=cex.cor * 0.9) ## r
}

## for the function of pairs 
panel.s.cor <- function(x, y, digits=4, prefix="", ...){
	usr <- par("usr"); on.exit(par(usr))
	r <- cor(x, y, use="pairwise.complete.obs", method='s')
	par(usr=c(0, 1, 0, 1))
	txt <- ifelse(r >= 0, 
				sprintf(paste0("%.", digits, "f"), r), ## format(c(r, 0.123456789), digits=digits)[1]
				txt <- sprintf(paste0("%.", digits-1, "f"), r)
			)
	txt <- paste(prefix, txt, sep="")
	cex.cor <- 0.8/strwidth(txt)
	text(0.5, 0.5, txt, cex=cex.cor * 0.9) ## r
}


pdf("figures/Fig.2 - cohort.trait.relimp.pdf", width=10, heigh=10, pointsize=10)
for(biomass in c('FW', 'DW')){
	relimp.data <- subset(trait.relimp, Biomass==biomass)
	plot.data <- reshape(relimp.data, v.names="Imp", idvar=c("Trait", "Treatment"), timevar="Dataset", direction="wide")
	colnames(plot.data) <- gsub('Imp.', '', colnames(plot.data))
	control.data <- subset(plot.data, Treatment=='control')
	rownames(control.data) <- control.data$Trait
	colnames(control.data) <- paste0(colnames(control.data), '\ncontrol')
	stress.data <- subset(plot.data, Treatment=='stress')
	rownames(stress.data) <- stress.data$Trait
	colnames(stress.data) <- paste0(colnames(stress.data), '\nstress')
	plot.data <- cbind(control.data[, -c(1,2,3)], stress.data[, -c(1,2,3)])
	trait.col <- unlist(lapply(rownames(control.data), trait.color.map))
	pairs(plot.data, lower.panel=panel.smooth, upper.panel=panel.s.cor, main=biomass, col=trait.col, pch=16, cex=1.5)
	## pairs(apply(plot.data, 2, rank), lower.panel=panel.smooth, upper.panel=panel.p.cor, main=biomass, col=trait.col, pch=16, cex=2)
}
for(biomass in c('FW', 'DW')){
	for(treatment in treatments){
		relimp.data <- subset(trait.relimp, Biomass==biomass & Treatment==treatment)
		plot.data <- reshape(relimp.data, v.names="Imp", idvar=c("Trait"), timevar="Dataset", direction="wide")
		colnames(plot.data) <- gsub('Imp.', '', colnames(plot.data))
		trait.col <- unlist(lapply(plot.data[, 'Trait'], trait.color.map))
		if(ncol(plot.data) > 5){
			pairs(plot.data[, -c(1,2,3)], lower.panel=panel.smooth, upper.panel=panel.p.cor, main=paste0(biomass, ', ', treatment), col=trait.col, pch=16, cex=2)
		}else{
			op <- par(mfrow=c(2, 2), pty="m")
			p.cor <- format(cor(plot.data[, 4], plot.data[, 5]), dig=4)
			lgd <- substitute(italic(r) == RR, list(RR=p.cor))
			plot(plot.data[, 4], plot.data[, 5], cex=2, main=paste0(biomass, ', ', treatment), col=trait.col, pch=16, xlab=colnames(plot.data)[4], ylab=colnames(plot.data)[5])
			my.legend('topleft', legend=lgd, pch=NA)
			par(op)
		}
	}
}
dev.off()

save.image('barley.regression.RData')
