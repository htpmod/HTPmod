################################################
# Author: Dijun Chen (chend@ipk-gatersleben.de)
# Update on Mar. 12, 2015
################################################
rm(list=ls())

if(file.exists('my.code.R')){
	source('my.code.R')
}

options(stringsAsFactors=FALSE)

manual.data <- read.delim("../barley/manual.data.1121KN.txt", header=TRUE, stringsAsFactors=FALSE)
rownames(manual.data) <- manual.data$PlantID

## loading the image data
image.data <- read.csv("../barley/report_1121KN_v2.csv", header=TRUE, sep=";", check.names=FALSE, stringsAsFactors=FALSE)
## change "Plant ID" to "PlantID"
colnames(image.data)[grep("Plant ID", colnames(image.data))] <- "PlantID"
## change "Day (Int)" to "DAS"
colnames(image.data)[grep("Day \\(Int\\)", colnames(image.data))] <- "DAS"
## change "Day (Float)" to "Realtime"
colnames(image.data)[grep("Day \\(Float\\)", colnames(image.data))] <- "Realtime"

common.columns <- c('PlantID', 'Genotype', 'Treatment', 'DAS', 'Realtime')
trait.columns <- colnames(image.data)[-(1:38)]
trait.columns <- trait.columns[-grep('stddev|skewness|wue', trait.columns)]

###################################################
system.time(
	trait.reproducibility <- trait.reproducibility.calculation(image.data, trait.columns)
)

pdf("Figure 1 - 1121KN.reproducible.traits.pdf", width=10, heigh=10, pointsize=10)
op <- par(mfrow=c(2, 2), pty="m")
rep.control.traits <- trait.reproducibility.plot(trait.reproducibility, 'control', cutoff.p=1, cutoff.cor=0.6)
rep.stress.traits <- trait.reproducibility.plot(trait.reproducibility, 'stress', cutoff.p=1, cutoff.cor=0.6)
par(op)
dev.off()

rep.traits <- intersect(rep.control.traits, rep.stress.traits)

###################################################
last.data <- subset(image.data, DAS==max(DAS))

last.data <- last.data[, c(common.columns, trait.columns)]
last.data[last.data == 0] <- NA

dim(last.data)
## remove columns or rows with less than 20% missing values
last.data <- last.data[, colSums(is.na(last.data)) < nrow(last.data)/5]
last.data <- last.data[rowSums(is.na(last.data)) < ncol(last.data)/5, ]

last.traits <- sort(intersect(trait.columns, colnames(last.data)))
genotypes <- sort(unique(last.data$Genotype))
treatments <- sort(unique(last.data$Treatment))

## outlier detection and filled by average values
result <- input <- last.data
for(genotype in genotypes){
	if(length(treatments) > 0){
		for(treatment in treatments){
			plant.rows <- which(input$Genotype %in% genotype & 
								input$Treatment %in% treatment)
			plant.data <- input[plant.rows, last.traits]
			result[plant.rows, last.traits] <- grubbs.test.fill.outliers(plant.data)
		}
	}else{
		plant.rows <- which(input$Genotype %in% genotype)
		plant.data <- input[plant.rows, last.traits]
		result[plant.rows, last.traits] <- grubbs.test.fill.outliers(plant.data)
	}
}
rownames(result) <- result$PlantID
sum(is.na(result))

## fill NAs by average values
null.index <- which(is.na(result), arr.ind=TRUE)
for(ind in 1:nrow(null.index)){
	genotype <- result[null.index[ind, 'row'], 'Genotype']
	treatment <- result[null.index[ind, 'row'], 'Treatment']
	geno.data <- result[result$Genotype==genotype & result$Treatment==treatment, null.index[ind, 'col']]
	if(all(is.na(geno.data))){
		result[null.index[ind, 'row'], null.index[ind, 'col']] <- mean(geno.data, na.rm=TRUE)
	}
	result[null.index[ind, 'row'], null.index[ind, 'col']] <- mean(geno.data, na.rm=TRUE)
}

if(sum(is.na(result)) > 0){
	rm.cols <- unique(which(is.na(result), arr.ind=TRUE)[, 'col'])
	result <- result[, -rm.cols]
}

sum(is.na(result))
traits <- sort(intersect(intersect(last.traits, colnames(result)), rep.traits))
write.table(traits, '1121KN.traits.out', sep='\t', dec='.', na='', row.names=FALSE, col.names=FALSE)

#################################################
filtered.traits <- scan("1121KN.filtered.traits.txt", what=character()) #####
# result$Treatment == manual.data[rownames(result), 'Treatment']
# result$Genotype == manual.data[rownames(result), 'Genotype']
filtered.data <- cbind(manual.data[rownames(result), c("PlantID", "Genotype", "Treatment", "FW", "DW")], abs(result[, filtered.traits]))
write.table(filtered.data, '1121KN.data.csv', sep=';', dec='.', na='', row.names=FALSE)

save.image('barley.input.RData')


#######################################################################
rm(list=ls())

library(plotrix)
if(file.exists('my.code.R')){
	source('my.code.R')
}
source('models.R')

KN1121 <- read.csv("1121KN.data.csv", header=TRUE, sep=";", stringsAsFactors=FALSE) # check.names=FALSE, 

raw.control <- KN1121[KN1121$Treatment=='control', ]
raw.stress <- KN1121[KN1121$Treatment=='stress', ]
raw.data <- list(Control=raw.control, Stress=raw.stress, "Control+Stress"=KN1121)

treatment.cols <- c("#e4007f", "#f8b62d", "f39800")
names(treatment.cols) <- names(raw.data)

pdf('Figure 1 - KN1121.biomass.correlation.pdf', width=10, heigh=10, pointsize=10)
ref.cors <- data.frame(dataset=names(raw.data), FW=0, DW=0, stringsAsFactors=FALSE)
rownames(ref.cors) <- names(raw.data)
for(dataset in names(raw.data)){
	plot.data <- raw.data[[dataset]][, c('FW', 'DW', 'volume.fluo.iap')]
	pairs(plot.data, main=dataset, lower.panel=panel.smooth, upper.panel=panel.cor) 
	cors <- cor(plot.data)
	ref.cors[dataset, 'FW'] <- cors['FW', 'volume.fluo.iap']
	ref.cors[dataset, 'DW'] <- cors['DW', 'volume.fluo.iap']
}
dev.off()
write.table(ref.cors, 'Table 1 - biomass.correlation-ref.cors.csv', sep=';', dec='.', na='', row.names=FALSE)

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
		pdf(paste0("Figure 1 - ", dataset, ".models.pdf"), width=10, heigh=10, pointsize=10)
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
pdf('Figure 1 - KN1121.model.evaluation.pdf', width=10, heigh=10, pointsize=10)
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

write.table(eval.table, 'Table 1 - evaluation1.data-eval.table.csv', sep=';', dec='.', na='', row.names=FALSE)

pdf('Figure 1 - KN1121.model.evaluation1.pdf', width=10, heigh=10, pointsize=10)
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
	ylim[1] <- 0
	ylim[2] <- ylim[2]+max(plot.std)*1.1
	barx <- barplot(plot.avg, ylim=ylim, beside=TRUE, col=model.cols[rownames(plot.avg)], main=paste("Evaluation based on", bio), ylab="RMSRE", xpd=FALSE, border=NA)
	text(barx, plot.avg+plot.std+abs(par('usr')[4]-par('usr')[3])*0.02, cex=0.8, labels=format(plot.avg, dig=2))
	error.bar(barx, plot.avg, plot.std, lwd=0.8, col=model.cols[rownames(plot.avg)])
	# my.legend("topleft", legend=names(model.cols), ncol=length(model.cols), col=model.cols, text.col=model.cols, pch=15, title="Models")
}
plot(1, type="n", xaxt="n", yaxt="n", frame=FALSE, xlab="", ylab="")
my.legend("topleft", legend=names(model.cols), ncol=length(model.cols), col=model.cols, text.col=model.cols, pch=15, title="Models")
par(op)
dev.off()


########## 
model <- 'rf.regression' ## the 'best' model
input <- KN1121.data

trait.effect <- data.frame(stringsAsFactors=FALSE)
for(biomass in c('FW', 'DW')){
	current <- barley.traits
	ind <- 1
	while(length(current) > 1){
		R2s <- rep(0, length(current))
		names(R2s) <- current
		for(trait in current){
			sub.traits <- setdiff(current, trait)
			reg.data <- input[, c(biomass, sub.traits)]
			reg <- regression.cross.validation(reg.data, model=model)
			R2s[trait] <- reg$R2
		}
		min.eff.trait <- names(which.max(R2s)[1])
		trait.effect <- rbind(trait.effect, 
						data.frame(biomass=biomass, index=ind, 
						avg=mean(R2s), std=std.error(R2s), 
						trait=min.eff.trait)
						)
		current <- setdiff(current, min.eff.trait)
		ind <- 1 + ind
	}
}
## write.table(trait.effect, 'Table 1 - trait.effect.data.csv', sep=';', dec='.', na='', row.names=FALSE)
trait.effect <- read.table('Table 1 - trait.effect.data.csv', sep=';', header=TRUE, stringsAsFactors=FALSE)

biomass.cols <- c("#ea5514", "#c9a063")
names(biomass.cols) <- c('FW', 'DW')
pdf('Figure 1 - KN1121.trait.effect.pdf', width=10, heigh=10, pointsize=10)
op <- par(mfrow=c(2, 1), pty="m")
for(bio in c('FW', 'DW')){
	sub.data <- subset(trait.effect, biomass==bio)
	barx <- barplot(sub.data$avg, ylim=c(0.9, 1), col=biomass.cols[bio], main=bio, ylab=expression(R^2), xpd=FALSE, border=NA)
	text(barx, sub.data$avg+sub.data$std+abs(par('usr')[4]-par('usr')[3])*0.02, cex=0.8, labels=trait.index[sub.data$trait])
	error.bar(barx, sub.data$avg, sub.data$std, lwd=0.8, col=biomass.cols[bio])
}
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

write.table(category.effect, 'Table 1 - evaluation2.data-category.effect.csv', sep=';', dec='.', na='', row.names=FALSE)
write.table(trait.prediction, 'Table 1 - evaluation2.data-trait.prediction.csv', sep=';', dec='.', na='', row.names=FALSE)
write.table(trait.relimp, 'Table 1 - evaluation2.data-trait.relimp.csv', sep=';', dec='.', na='', row.names=FALSE)

### test the predictive power of models 
pdf('Figure 1 - KN1121.model.evaluation2.pdf', width=10, heigh=10, pointsize=10)
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
write.table(view.effect, 'Table 1 - evaluation3.data-view.effect.csv', sep=';', dec='.', na='', row.names=FALSE)


### test the predictive power of models 
pdf('Figure 1 - KN1121.model.evaluation3.pdf', width=10, heigh=10, pointsize=10)
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

write.table(trait.relimp2, 'Table 1 - evaluation2.data-trait.relimp2.csv', sep=';', dec='.', na='', row.names=FALSE)


save.image('barley.regression.RData')
