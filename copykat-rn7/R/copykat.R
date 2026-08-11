# copykat (rat fork, genome="rn7") — vendored from copykat 1.1.0 commit b795ff7.
# Modifications: rn7 annotation branch, annotateGenes.rn7(), mm10 gene-space block widened to rn7.

copykat <- 
function(rawmat=rawdata, id.type="S", cell.line="no", ngene.chr=5,LOW.DR=0.05, UP.DR=0.1, win.size=25, norm.cell.names="", KS.cut=0.1, sam.name="", distance="euclidean", output.seg="FALSE", plot.genes="TRUE", genome="hg20", n.cores=1){

start_time <- Sys.time()
  set.seed(1234)
  sample.name <- paste(sam.name,"_copykat_", sep="")

  print("running copykat v1.1.0")

  print("step1: read and filter data ...")
  print(paste(nrow(rawmat), " genes, ", ncol(rawmat), " cells in raw data", sep=""))

  genes.raw <- apply(rawmat, 2, function(x)(sum(x>0)))

  if(sum(genes.raw> 200)==0) stop("none cells have more than 200 genes")
  if(sum(genes.raw<100)>1){
    rawmat <- rawmat[, -which(genes.raw< 200)]
    print(paste("filtered out ", sum(genes.raw<=200), " cells with less than 200 genes; remaining ", ncol(rawmat), " cells", sep=""))
  }
  ##
  der<- apply(rawmat,1,function(x)(sum(x>0)))/ncol(rawmat)

  if(sum(der>LOW.DR)>=1){
    rawmat <- rawmat[which(der > LOW.DR), ]; print(paste(nrow(rawmat)," genes past LOW.DR filtering", sep=""))
  }

  WNS1 <- "data quality is ok"
  if(nrow(rawmat) < 7000){
    WNS1 <- "low data quality"
    UP.DR<- LOW.DR
    print("WARNING: low data quality; assigned LOW.DR to UP.DR...")
  }

  print("step 2: annotations gene coordinates ...")
  if(genome=="hg20"){
  anno.mat <- annotateGenes.hg20(mat = rawmat, ID.type = id.type) #SYMBOL or ENSEMBLE
  } else if(genome=="mm10"){
  anno.mat <- annotateGenes.mm10(mat = rawmat, ID.type = id.type) #SYMBOL or ENSEMBLE
  dim(rawmat)
  } else if(genome=="rn7"){
  anno.mat <- annotateGenes.rn7(mat = rawmat, ID.type = id.type) #SYMBOL or ENSEMBLE (rat mRatBN7.2)
  dim(rawmat)
  }
  anno.mat <- anno.mat[order(as.numeric(anno.mat$abspos), decreasing = FALSE),]

# print(paste(nrow(anno.mat)," genes annotated", sep=""))

  ### module 3 removing genes that are involved in cell cycling

  if(genome=="hg20"){
  HLAs <- anno.mat$hgnc_symbol[grep("^HLA-", anno.mat$hgnc_symbol)]
  toRev <- which(anno.mat$hgnc_symbol %in% c(as.vector(cyclegenes[[1]]), HLAs))
  if(length(toRev)>0){
    anno.mat <- anno.mat[-toRev, ]
  }
  }
#  print(paste(nrow(anno.mat)," genes after rm cell cycle genes", sep=""))
  ### secondary filtering
  ToRemov2 <- NULL
  for(i in 8:ncol(anno.mat)){
    cell <- cbind(anno.mat$chromosome_name, anno.mat[,i])
    cell <- cell[cell[,2]!=0,]
    if(length(as.numeric(cell))< 5){
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    } else if(length(rle(cell[,1])$length)<length(unique((anno.mat$chromosome_name)))|min(rle(cell[,1])$length)< ngene.chr){
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    }
    i<- i+1
  }

  if(length(ToRemov2)==(ncol(anno.mat)-7)) stop("all cells are filtered")
  if(length(ToRemov2)>0){
    anno.mat <-anno.mat[, -which(colnames(anno.mat) %in% ToRemov2)]
  }

  # print(paste("filtered out ", length(ToRemov2), " cells with less than ",ngene.chr, " genes per chr", sep=""))
  rawmat3 <- data.matrix(anno.mat[, 8:ncol(anno.mat)])
  norm.mat<- log(sqrt(rawmat3)+sqrt(rawmat3+1))
  norm.mat<- apply(norm.mat,2,function(x)(x <- x-mean(x)))
  colnames(norm.mat) <-  colnames(rawmat3)

  #print(paste("A total of ", ncol(norm.mat), " cells, ", nrow(norm.mat), " genes after preprocessing", sep=""))

  ##smooth data
  print("step 3: smoothing data with dlm ...")
  dlm.sm <- function(c){
    model <- dlm::dlmModPoly(order=1, dV=0.16, dW=0.001)
    x <- dlm::dlmSmooth(norm.mat[, c], model)$s
    x<- x[2:length(x)]
    x <- x-mean(x)
  }

  test.mc <-parallel::mclapply(1:ncol(norm.mat), dlm.sm, mc.cores = n.cores)
  norm.mat.smooth <- matrix(unlist(test.mc), ncol = ncol(norm.mat), byrow = FALSE)

  colnames(norm.mat.smooth) <- colnames(norm.mat)

  print("step 4: measuring baselines ...")
  if (cell.line=="yes"){
  	print("running pure cell line mode")
  	    relt <- baseline.synthetic(norm.mat=norm.mat.smooth, min.cells=10, n.cores=n.cores)
		norm.mat.relat <- relt$expr.relat
		CL <- relt$cl
        WNS <- "run with cell line mode"
    	preN <- NULL

      } else if(length(norm.cell.names)>1){

        #print(paste(length(norm.cell.names), "normal cells provided", sep=""))
         NNN <- length(colnames(norm.mat.smooth)[which(colnames(norm.mat.smooth) %in% norm.cell.names)])
         print(paste(NNN, " known normal cells found in dataset", sep=""))

         if (NNN==0) stop("known normal cells provided; however none existing in testing dataset")
         print("run with known normal...")

         basel <- apply(norm.mat.smooth[, which(colnames(norm.mat.smooth) %in% norm.cell.names)],1,median); print("baseline is from known input")

          d <- parallelDist::parDist(t(norm.mat.smooth),threads =n.cores, method="euclidean") ##use smooth and segmented data to detect intra-normal cells

          km <- 6
          fit <- hclust(d, method="ward.D2")
           CL <- cutree(fit, km)

           while(!all(table(CL)>5)){
          km <- km -1
          CL <- cutree(fit, k=km)
         if(km==2){
         break
         }
         }

        WNS <- "run with known normal"
        preN <- norm.cell.names
         ##relative expression using pred.normal cells
      	norm.mat.relat <- norm.mat.smooth-basel

        }else {
         basa <- baseline.norm.cl(norm.mat.smooth=norm.mat.smooth, min.cells=5, n.cores=n.cores)
          basel <- basa$basel
          WNS <- basa$WNS
          preN <- basa$preN
          CL <- basa$cl
          if (WNS =="unclassified.prediction"){

                    basa <- baseline.GMM(CNA.mat=norm.mat.smooth, max.normal=5, mu.cut=0.05, Nfraq.cut=0.99,RE.before=basa,n.cores=n.cores)
                    basel <-basa$basel
                    WNS <- basa$WNS

                    preN <- basa$preN

              }
          ##relative expression using pred.normal cells
             norm.mat.relat <- norm.mat.smooth-basel

             }

  ###use a smaller set of genes to perform segmentation
  DR2 <- apply(rawmat3,1,function(x)(sum(x>0)))/ncol(rawmat3)
  ##relative expression using pred.normal cells
  norm.mat.relat <- norm.mat.relat[which(DR2>=UP.DR),]

  ###filter cells
  anno.mat2 <- anno.mat[which(DR2>=UP.DR), ]

  ToRemov3 <- NULL
  for(i in 8:ncol(anno.mat2)){
    cell <- cbind(anno.mat2$chromosome_name, anno.mat2[,i])
    cell <- cell[cell[,2]!=0,]
    if(length(as.numeric(cell))< 5){
      rm <- colnames(anno.mat2)[i]
      ToRemov3 <- c(ToRemov3, rm)
    } else if(length(rle(cell[,1])$length)<length(unique((anno.mat$chromosome_name)))|min(rle(cell[,1])$length)< ngene.chr){
      rm <- colnames(anno.mat2)[i]
      ToRemov3 <- c(ToRemov3, rm)
    }
    i<- i+1
  }

  if(length(ToRemov3)==ncol(norm.mat.relat)) stop ("all cells are filtered")

  if(length(ToRemov3)>0){
    norm.mat.relat <-norm.mat.relat[, -which(colnames(norm.mat.relat) %in% ToRemov3)]
   #print(paste("filtered out ", length(ToRemov3), " cells with less than ",ngene.chr, " genes per chr", sep=""))
  }

  #print(paste("final segmentation: ", nrow(norm.mat.relat), " genes; ", ncol(norm.mat.relat), " cells", sep=""))

  CL <- CL[which(names(CL) %in% colnames(norm.mat.relat))]
  CL <- CL[order(match(names(CL), colnames(norm.mat.relat)))]

  print("step 5: segmentation...")
  results <- CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = KS.cut, n.cores=n.cores)

  if(length(results$breaks)<25){
    print("too few breakpoints detected; decreased KS.cut to 50%")
    results <- CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = 0.5*KS.cut, n.cores=n.cores)
  }

  if(length(results$breaks)<25){
    print("too few breakpoints detected; decreased KS.cut to 75%")
    results <- CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = 0.5*0.5*KS.cut, n.cores=n.cores)
  }

  if(length(results$breaks)<25) stop ("too few segments; try to decrease KS.cut; or improve data")

  colnames(results$logCNA) <- colnames(norm.mat.relat)
  results.com <- apply(results$logCNA,2, function(x)(x <- x-mean(x)))
  RNA.copycat <- cbind(anno.mat2[, 1:7], results.com)

  write.table(RNA.copycat, paste(sample.name, "CNA_raw_results_gene_by_cell.txt", sep=""), sep="\t", row.names = FALSE, quote = F)

  if(genome=="hg20"){
  print("step 6: convert to genomic bins...") ###need multi-core
  Aj <- convert.all.bins.hg20(DNA.mat = DNA.hg20, RNA.mat=RNA.copycat, n.cores = n.cores)

  uber.mat.adj <- data.matrix(Aj$RNA.adj[, 4:ncol(Aj$RNA.adj)])

  print("step 7: adjust baseline ...")

    if(cell.line=="yes"){

               mat.adj <- data.matrix(Aj$RNA.adj[, 4:ncol(Aj$RNA.adj)])
               write.table(cbind(Aj$RNA.adj[, 1:3], mat.adj), paste(sample.name, "CNA_results.txt", sep=""), sep="\t", row.names = FALSE, quote = F)

                if(distance=="euclidean"){
                 hcc <- hclust(parallelDist::parDist(t(mat.adj),threads =n.cores, method = distance), method = "ward.D")
                  }else {
                 hcc <- hclust(as.dist(1-cor(mat.adj, method = distance)), method = "ward.D")
                   }


                  saveRDS(hcc, file = paste(sample.name,"clustering_results.rds",sep=""))

                   #plot heatmap
                   print("step 8: ploting heatmap ...")
                  my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

                   chr <- as.numeric(Aj$DNA.adj$chrom) %% 2+1
                   rbPal1 <- colorRampPalette(c('black','grey'))
                   CHR <- rbPal1(2)[as.numeric(chr)]
                   chr1 <- cbind(CHR,CHR)


                   if (ncol(mat.adj)< 3000){
                   h <- 10
                   } else {
                   h <- 15
                     }

                  col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))
                  #library(parallelDist)

                   if(distance=="euclidean"){
                          jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
                          heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
                          ColSideColors=chr1,Colv=NA, Rowv=TRUE,
                          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
                          keysize=1, density.info="none", trace="none",
                          cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
                          symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
                          dev.off()
                          ### add a step to plot out gene by cell matrix
                          if(plot.genes=="TRUE"){

                          rownames(results.com) <- anno.mat2$hgnc_symbol
                          chrg <- as.numeric(anno.mat2$chrom) %% 2+1
                          rbPal1g <- colorRampPalette(c('black','grey'))
                          CHRg <- rbPal1(2)[as.numeric(chrg)]
                          chr1g <- cbind(CHRg,CHRg)

                          pdf(paste(sample.name,"with_genes_heatmap.pdf",sep=""), height=h*2.5, width=40)
                          heatmap.3(t(results.com),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
                          ColSideColors=chr1g,Colv=NA, Rowv=TRUE,
                          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
                          keysize=1, density.info="none", trace="none",
                          cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
                          symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
                          dev.off()
                           }
                         #end of ploting gene by cell matrix

                } else {
                          jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
                          heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
                          ColSideColors=chr1,Colv=NA, Rowv=TRUE,
                          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
                          keysize=1, density.info="none", trace="none",
                          cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
                          symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
                          dev.off()
                           ### add a step to plot out gene by cell matrix
             if(plot.genes=="TRUE"){

                          rownames(results.com) <- anno.mat2$hgnc_symbol
                          chrg <- as.numeric(anno.mat2$chrom) %% 2+1
                          rbPal1g <- colorRampPalette(c('black','grey'))
                          CHRg <- rbPal1(2)[as.numeric(chrg)]
                          chr1g <- cbind(CHRg,CHRg)

                          pdf(paste(sample.name,"with_genes_heatmap.pdf",sep=""), height=h*2.5, width=40)
                          heatmap.3(t(results.com),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
                          ColSideColors=chr1g,Colv=NA, Rowv=TRUE,
                          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
                          keysize=1, density.info="none", trace="none",
                          cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
                          symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
                          dev.off()
                          }
                         #end of ploting gene by cell matrix
                          }

                          end_time<- Sys.time()
                          print(end_time -start_time)

                         reslts <- list(cbind(Aj$RNA.adj[, 1:3], mat.adj), hcc)
                         names(reslts) <- c("CNAmat","hclustering")
                         return(reslts)
    } else {
      ########## cell line mode ends here ####################

      #removed baseline adjustment
        if(distance=="euclidean"){
        hcc <- hclust(parallelDist::parDist(t(uber.mat.adj),threads =n.cores, method = distance), method = "ward.D")
        }else {
        hcc <- hclust(as.dist(1-cor(uber.mat.adj, method = distance)), method = "ward.D")
        }
        hc.umap <- cutree(hcc,2)
        names(hc.umap) <- colnames(results.com)

        cl.ID <- NULL
        for(i in 1:max(hc.umap)){
        cli <- names(hc.umap)[which(hc.umap==i)]
        pid <- length(intersect(cli, preN))/length(cli)
        cl.ID <- c(cl.ID, pid)
        i<- i+1
         }

        com.pred <- names(hc.umap)
        com.pred[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "diploid"
        com.pred[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "aneuploid"
        names(com.pred) <- names(hc.umap)

  ################removed baseline adjustment
        results.com.rat <- uber.mat.adj-apply(uber.mat.adj[,which(com.pred=="diploid")], 1, mean)
        results.com.rat <- apply(results.com.rat,2,function(x)(x <- x-mean(x)))
        results.com.rat.norm <- results.com.rat[,which(com.pred=="diploid")]; dim(results.com.rat.norm)

        cf.h <- apply(results.com.rat.norm, 1, sd)
        base <- apply(results.com.rat.norm, 1, mean)

        adjN <- function(j){
        a <- results.com.rat[, j]
        a[abs(a-base) <= 0.25*cf.h] <- mean(a)
        a
        }


        mc.adjN <-  parallel::mclapply(1:ncol(results.com.rat),adjN, mc.cores = n.cores)
        adj.results <- matrix(unlist(mc.adjN), ncol = ncol(results.com.rat), byrow = FALSE)
        colnames(adj.results) <- colnames(results.com.rat)

        #rang <- 0.5*(max(adj.results)-min(adj.results))
        #mat.adj <- adj.results/rang
        mat.adj <- t(t(adj.results)-apply(adj.results,2,mean))

        print("step 8: final prediction ...")

        if(distance=="euclidean"){
         hcc <- hclust(parallelDist::parDist(t(mat.adj),threads =n.cores, method = distance), method = "ward.D")
         }else {
         hcc <- hclust(as.dist(1-cor(mat.adj, method = distance)), method = "ward.D")
         }

         hc.umap <- cutree(hcc,2)
         names(hc.umap) <- colnames(results.com)

        saveRDS(hcc, file = paste(sample.name,"clustering_results.rds",sep=""))

        cl.ID <- NULL
        for(i in 1:max(hc.umap)){
        cli <- names(hc.umap)[which(hc.umap==i)]
        pid <- length(intersect(cli, preN))/length(cli)
        cl.ID <- c(cl.ID, pid)
        i<- i+1
         }

        com.preN <- names(hc.umap)
        com.preN[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "diploid"
        com.preN[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "aneuploid"
        names(com.preN) <- names(hc.umap)

        if(WNS=="unclassified.prediction"){
        com.preN[which(com.preN == "diploid")] <- "c1:diploid:low.conf"
        com.preN[which(com.preN == "aneuploid")] <- "c2:aneuploid:low.conf"
        }

      print("step 9: saving results...")

  ##add back filtered cells as not defined in prediction results
  '%!in%' <- function(x,y)!('%in%'(x,y))

  ndef <- colnames(rawmat)[which(colnames(rawmat) %!in% names(com.preN))]
  if(length(ndef)>0){
    res <- data.frame(cbind(c(names(com.preN),ndef), c(com.preN, rep("not.defined",length(ndef)))))
    colnames(res) <- c("cell.names", "copykat.pred")
  } else {
    res <- data.frame(cbind(names(com.preN), com.preN))
    colnames(res) <- c("cell.names", "copykat.pred")
  }
  ##end
  write.table(res, paste(sample.name, "prediction.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)

  ####save copycat CNA
  write.table(cbind(Aj$RNA.adj[, 1:3], mat.adj), paste(sample.name, "CNA_results.txt", sep=""), sep="\t", row.names = FALSE, quote = F)

  ####%%%%%%%%%%%%%%%%%next heatmaps, subpopulations and tSNE overlay
  print("step 10: ploting heatmap ...")
  my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

  chr <- as.numeric(Aj$DNA.adj$chrom) %% 2+1
  rbPal1 <- colorRampPalette(c('black','grey'))
  CHR <- rbPal1(2)[as.numeric(chr)]
  chr1 <- cbind(CHR,CHR)

  rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
  compreN_pred <- rbPal5(2)[as.numeric(factor(com.preN))]

  cells <- rbind(compreN_pred,compreN_pred)

  if (ncol(mat.adj)< 3000){
    h <- 10
  } else {
    h <- 15
  }

  col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

  if(distance=="euclidean"){
  jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
   heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
            ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))

  legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=1)
  dev.off()

  ### add a step to plot out gene by cell matrix
  if(plot.genes=="TRUE"){
    dim(results.com)
    rownames(results.com) <- anno.mat2$hgnc_symbol
    chrg <- as.numeric(anno.mat2$chrom) %% 2+1
    rbPal1g <- colorRampPalette(c('black','grey'))
    CHRg <- rbPal1(2)[as.numeric(chrg)]
    chr1g <- cbind(CHRg,CHRg)

    pdf(paste(sample.name,"with_genes_heatmap.pdf",sep=""), height=h*2.5, width=40)
    heatmap.3(t(results.com),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
              ColSideColors=chr1g,RowSideColors=cells,Colv=NA, Rowv=TRUE,
              notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
              keysize=1, density.info="none", trace="none",
              cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
              symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
    dev.off()
  }
  #end of ploting gene by cell matrix



  } else {
    jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
    heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
                 ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
              notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
              keysize=1, density.info="none", trace="none",
              cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
              symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))

    legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=1)

    dev.off()
    ### add a step to plot out gene by cell matrix
    if(plot.genes=="TRUE"){
      dim(results.com)
      rownames(results.com) <- anno.mat2$hgnc_symbol
      chrg <- as.numeric(anno.mat2$chrom) %% 2+1
      rbPal1g <- colorRampPalette(c('black','grey'))
      CHRg <- rbPal1(2)[as.numeric(chrg)]
      chr1g <- cbind(CHRg,CHRg)

      pdf(paste(sample.name,"with_genes_heatmap.pdf",sep=""), height=h*2.5, width=40)
      heatmap.3(t(results.com),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
                ColSideColors=chr1g,RowSideColors=cells, Colv=NA, Rowv=TRUE,
                notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
                keysize=1, density.info="none", trace="none",
                cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
                symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
      dev.off()
    }
    #end of ploting gene by cell matrix
  }

 if(output.seg=="TRUE"){
  print("generating seg files for IGV viewer")

  thisRatio <- cbind(Aj$RNA.adj[, 1:3], mat.adj)
  Short <- NULL
  chr <- rle(thisRatio$chrom)[[2]]

  for(c in 4:ncol(thisRatio))
  {
    for (x in 1:length(chr)){
      thisRatio.sub <- thisRatio[which(thisRatio$chrom==chr[x]), ]
      seg.mean.sub <- rle(thisRatio.sub[,c])[[2]]

      rle.length.sub <- rle(thisRatio.sub[,c])[[1]]

      num.mark.sub <- seq(1,length(rle.length.sub),1)
      loc.start.sub <-seq(1,length(rle.length.sub),1)
      loc.end.sub <- seq(1,length(rle.length.sub),1)

      len <-0
      j <-1

      for (j in 1: length(rle.length.sub)){
        num.mark.sub[j] <- rle.length.sub[j]
        loc.start.sub[j] <- thisRatio.sub$chrompos[len+1]
        len <- num.mark.sub[j]+len
        loc.end.sub[j] <- thisRatio.sub$chrompos[len]
        j <- j+1
      }

      ID <- rep(colnames(thisRatio[c]), times=length(rle.length.sub))
      chrom <- rep(chr[x], times=length(rle.length.sub))
      Short.sub <- cbind(ID,chrom,loc.start.sub,loc.end.sub,num.mark.sub,seg.mean.sub)
      Short <- rbind(Short, Short.sub)
      x <- x+1
    }
    c<- c+1
  }

  colnames(Short) <- c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")
  head(Short)
  write.table(Short, paste(sample.name, "CNA_results.seg", sep=""), row.names = FALSE, quote=FALSE, sep="\t")

}
  end_time<- Sys.time()
  print(end_time -start_time)

  reslts <- list(res, cbind(Aj$RNA.adj[, 1:3], mat.adj), hcc)
  names(reslts) <- c("prediction", "CNAmat","hclustering")
  return(reslts)
}

  }

  if(genome=="mm10" || genome=="rn7") {
    uber.mat.adj <- data.matrix(results.com)
    dim(uber.mat.adj)
    if(distance=="euclidean"){
      hcc <- hclust(parallelDist::parDist(t(uber.mat.adj),threads =n.cores, method = distance), method = "ward.D")
    }else {
      hcc <- hclust(as.dist(1-cor(uber.mat.adj, method = distance)), method = "ward.D")
    }
    hc.umap <- cutree(hcc,2)
    names(hc.umap) <- colnames(results.com)

    cl.ID <- NULL
    for(i in 1:max(hc.umap)){
      cli <- names(hc.umap)[which(hc.umap==i)]
      pid <- length(intersect(cli, preN))/length(cli)
      cl.ID <- c(cl.ID, pid)
      i<- i+1
    }

    com.pred <- names(hc.umap)
    com.pred[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "diploid"
    com.pred[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "aneuploid"
    names(com.pred) <- names(hc.umap)

    ################removed baseline adjustment
    results.com.rat <- uber.mat.adj-apply(uber.mat.adj[,which(com.pred=="diploid")], 1, mean)

    results.com.rat <- apply(results.com.rat,2,function(x)(x <- x-mean(x)))
    results.com.rat.norm <- results.com.rat[,which(com.pred=="diploid")]; dim(results.com.rat.norm)

    cf.h <- apply(results.com.rat.norm, 1, sd)
    base <- apply(results.com.rat.norm, 1, mean)

    adjN <- function(j){
      a <- results.com.rat[, j]
      a[abs(a-base) <= 0.25*cf.h] <- mean(a)
      a
    }


    mc.adjN <-  parallel::mclapply(1:ncol(results.com.rat),adjN, mc.cores = n.cores)
    adj.results <- matrix(unlist(mc.adjN), ncol = ncol(results.com.rat), byrow = FALSE)
    colnames(adj.results) <- colnames(results.com.rat)

    #rang <- 0.5*(max(adj.results)-min(adj.results))
    #mat.adj <- adj.results/rang
    mat.adj <- t(t(adj.results)-apply(adj.results,2,mean))

    print("step 8: final prediction ...")

    if(distance=="euclidean"){
      hcc <- hclust(parallelDist::parDist(t(mat.adj),threads =n.cores, method = distance), method = "ward.D")
    }else {
      hcc <- hclust(as.dist(1-cor(mat.adj, method = distance)), method = "ward.D")
    }

    hc.umap <- cutree(hcc,2)
    names(hc.umap) <- colnames(results.com)

    saveRDS(hcc, file = paste(sample.name,"clustering_results.rds",sep=""))

    cl.ID <- NULL
    for(i in 1:max(hc.umap)){
      cli <- names(hc.umap)[which(hc.umap==i)]
      pid <- length(intersect(cli, preN))/length(cli)
      cl.ID <- c(cl.ID, pid)
      i<- i+1
    }

    com.preN <- names(hc.umap)
    com.preN[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "diploid"
    com.preN[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "aneuploid"
    names(com.preN) <- names(hc.umap)

    if(WNS=="unclassified.prediction"){
      com.preN[which(com.preN == "diploid")] <- "c1:diploid:low.conf"
      com.preN[which(com.preN == "aneuploid")] <- "c2:aneuploid:low.conf"
    }

    print("step 9: saving results...")

    ##add back filtered cells as not defined in prediction results
    '%!in%' <- function(x,y)!('%in%'(x,y))
    ndef <- colnames(rawmat)[which(colnames(rawmat) %!in% names(com.preN))]
    if(length(ndef)>0){
      res <- data.frame(cbind(c(names(com.preN),ndef), c(com.preN, rep("not.defined",length(ndef)))))
      colnames(res) <- c("cell.names", "copykat.pred")
    } else {
      res <- data.frame(cbind(names(com.preN), com.preN))
      colnames(res) <- c("cell.names", "copykat.pred")
    }
    ##end
    write.table(res, paste(sample.name, "prediction.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)

    ####save copycat CNA
    write.table(cbind(anno.mat2[, 1:7], mat.adj), paste(sample.name, "CNA_results.txt", sep=""), sep="\t", row.names = FALSE, quote = F)

    ####%%%%%%%%%%%%%%%%%next heatmaps, subpopulations and tSNE overlay
    print("step 10: ploting heatmap ...")
    my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

    rownames(mat.adj) <- anno.mat2$mgi_symbol
    chrg <- as.numeric(anno.mat2$chromosome_name) %% 2+1
    rle(as.numeric(anno.mat2$chromosome_name))
    rbPal1g <- colorRampPalette(c('black','grey'))
    CHRg <- rbPal1g(2)[as.numeric(chrg)]
    chr1g <- cbind(CHRg,CHRg)


    rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
    compreN_pred <- rbPal5(2)[as.numeric(factor(com.preN))]

    cells <- rbind(compreN_pred,compreN_pred)

    if (ncol(mat.adj)< 3000){
      h <- 10
    } else {
      h <- 15
    }

    col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

    if(distance=="euclidean"){

        pdf(paste(sample.name,"with_genes_heatmap.pdf",sep=""), height=h*2.5, width=40)
        heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
                  ColSideColors=chr1g,RowSideColors=cells,Colv=NA, Rowv=TRUE,
                  notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
                  keysize=1, density.info="none", trace="none",
                  cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
                  symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
        dev.off()


    } else {

        pdf(paste(sample.name,"with_genes_heatmap1.pdf",sep=""), height=h*2.5, width=40)
        heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
                  ColSideColors=chr1g,RowSideColors=cells,Colv=NA, Rowv=TRUE,
                  notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
                  keysize=1, density.info="none", trace="none",
                  cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
                  symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
        dev.off()

      #end of ploting gene by cell matrix
    }

    if(output.seg=="TRUE"){
      print("generating seg files for IGV viewer")

      thisRatio <- cbind(anno.mat2[, c(2,3,1)], mat.adj)
      Short <- NULL
      chr <- rle(thisRatio$chromosome_name)[[2]]

      for(c in 4:ncol(thisRatio))
      {
        for (x in 1:length(chr)){
          thisRatio.sub <- thisRatio[which(thisRatio$chromosome_name==chr[x]), ]
          seg.mean.sub <- rle(thisRatio.sub[,c])[[2]]

          rle.length.sub <- rle(thisRatio.sub[,c])[[1]]

          num.mark.sub <- seq(1,length(rle.length.sub),1)
          loc.start.sub <-seq(1,length(rle.length.sub),1)
          loc.end.sub <- seq(1,length(rle.length.sub),1)

          len <-0
          j <-1

          for (j in 1: length(rle.length.sub)){
            num.mark.sub[j] <- rle.length.sub[j]
            loc.start.sub[j] <- thisRatio.sub$start_position[len+1]
            len <- num.mark.sub[j]+len
            loc.end.sub[j] <- thisRatio.sub$start_position[len]
            j <- j+1
          }

          ID <- rep(colnames(thisRatio[c]), times=length(rle.length.sub))
          chrom <- rep(chr[x], times=length(rle.length.sub))
          Short.sub <- cbind(ID,chrom,loc.start.sub,loc.end.sub,num.mark.sub,seg.mean.sub)
          Short <- rbind(Short, Short.sub)
          x <- x+1
        }
        c<- c+1
      }

      colnames(Short) <- c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")

      write.table(Short, paste(sample.name, "CNA_results.seg", sep=""), row.names = FALSE, quote=FALSE, sep="\t")

    }
    end_time<- Sys.time()
    print(end_time -start_time)

    reslts <- list(res, cbind(anno.mat2[, 1:7], mat.adj), hcc)
    names(reslts) <- c("prediction", "CNAmat","hclustering")
    return(reslts)

  }
}

annotateGenes.hg20 <- 
function(mat, ID.type="S"){
  print("start annotation ...")

  if(substring(ID.type,1,1) %in% c("E", "e")){
    shar <- intersect(rownames(mat), full.anno$ensembl_gene_id)
    mat <- mat[which(rownames(mat) %in% shar),]
    anno <- full.anno[which(as.vector(full.anno$ensembl_gene_id) %in% shar),]
    anno <- anno[!duplicated(anno$hgnc_symbol),]
    anno <- anno[order(match(anno$ensembl_gene_id, rownames(mat))),]
    data <- cbind(anno, mat)

  }else if(substring(ID.type,1,1) %in% c("S", "s")) {

    shar <- intersect(rownames(mat), full.anno$hgnc_symbol)
    mat <- mat[which(rownames(mat) %in% shar),]
    anno <- full.anno[which(as.vector(full.anno$hgnc_symbol) %in% shar),]
    anno <- anno[!duplicated(anno$hgnc_symbol),]
    anno <- anno[order(match(anno$hgnc_symbol, rownames(mat))),]
    data <- cbind(anno, mat)
  }
}


annotateGenes.mm10 <- 
function(mat, ID.type="S", full.anno=full.anno.mm10){
  print("start annotation ...")

  if(substring(ID.type,1,1) %in% c("E", "e")){
    shar <- intersect(rownames(mat), full.anno$ensembl_gene_id)
    mat <- mat[which(rownames(mat) %in% shar),]
    anno <- full.anno[which(as.vector(full.anno$ensembl_gene_id) %in% shar),]
    anno <- anno[!duplicated(anno$mgi_symbol),]
    anno <- anno[order(match(anno$ensembl_gene_id, rownames(mat))),]
    data <- cbind(anno, mat)

  }else if(substring(ID.type,1,1) %in% c("S", "s")) {

    shar <- intersect(rownames(mat), full.anno$mgi_symbol)
    rownames(mat)[1:10]
    full.anno$mgi_symbol[]

    mat <- mat[which(rownames(mat) %in% shar),]
    anno <- full.anno[which(as.vector(full.anno$mgi_symbol) %in% shar),]
    anno <- anno[!duplicated(anno$mgi_symbol),]
    anno <- anno[order(match(anno$mgi_symbol, rownames(mat))),]
    data <- cbind(anno, mat)
  }
}


baseline.GMM <- 
function(CNA.mat, max.normal=5, mu.cut=0.05, Nfraq.cut=0.99, RE.before=basa, n.cores=1){

     N.normal <-NULL
     for(m in 1:ncol(CNA.mat)){

      print(paste("cell: ", m, sep=""))
      sam <- CNA.mat[, m]
      sg <- max(c(0.05, 0.5*sd(sam)));
      GM3 <- mixtools::normalmixEM(sam, lambda = rep(1,3)/3, mu = c(-0.2, 0, 0.2),sigma = sg,arbvar=FALSE,ECM=FALSE,maxit=500);#maxrestarts=10; arbmean=TRUE; arbvar=TRUE;epsilon=0.01

      ###decide normal or tumor
      s <- sum(abs(GM3$mu)<=mu.cut)

      if(s>=1){
        frq <- sum(GM3$lambda[which(abs(GM3$mu)<=mu.cut)])
     #   print(paste("N.fraq ", frq, sep=""))
     #   print(paste("sigma: ", GM3$sigma[1], sep=""))
        if(frq> Nfraq.cut){
          pred <- "diploid"
        }else{pred<-"aneuploid"}

      }else {pred <- "aneuploid"}
    #  print(paste("pred: ", pred, sep=""))
       N.normal<- c(N.normal,pred)

      if(sum(N.normal=="diploid")>=max.normal){break}
      m<- m+1
    }

    names(N.normal) <- colnames(CNA.mat)[1:length(N.normal)]
    preN <- names(N.normal)[which(N.normal=="diploid")]

    d <- parallelDist::parDist(t(CNA.mat), threads = n.cores) ##use smooth and segmented data to detect intra-normal cells
    km <- 6
    fit <- hclust(d, method="ward.D2")
    ct <- cutree(fit, k=km)


    if(length(preN) >2){
      WNS <- ""
      basel <- apply(CNA.mat[, which(colnames(CNA.mat) %in% preN)], 1, mean)

      RE <- list(basel, WNS, preN, ct)
      names(RE) <- c("basel", "WNS", "preN", "cl")
      return(RE)
    }else{
      return(RE.before) ##found this bug
    }

}


baseline.norm.cl <- 
function(norm.mat.smooth, min.cells=5, n.cores=n.cores){

  d <- parallelDist::parDist(t(norm.mat.smooth), threads = n.cores) ##use smooth and segmented data to detect intra-normal cells
  km <- 6
  fit <- hclust(d, method="ward.D2")
  ct <- cutree(fit, k=km)

  while(!all(table(ct)>min.cells)){
    km <- km -1
    ct <- cutree(fit, k=km)
    if(km==2){
      break
    }
  }

  SDM <-NULL
  SSD <-NULL
  for(i in min(ct):max(ct)){

    data.c <- apply(norm.mat.smooth[, which(ct==i)],1, median)
    sx <- max(c(0.05, 0.5*sd(data.c)))
    GM3 <- mixtools::normalmixEM(data.c, lambda = rep(1,3)/3, mu = c(-0.2, 0, 0.2), sigma = sx,arbvar=FALSE,ECM=FALSE,maxit=5000)
    SDM <- c(SDM, GM3$sigma[1])
    SSD <- c(SSD, sd(data.c))
       i <- i+1
      }

  wn <- mean(cluster::silhouette(cutree(fit, k=2), d)[, "sil_width"])

  ####
 PDt <- pf(max(SDM)^2/min(SDM)^2, nrow(norm.mat.smooth), nrow(norm.mat.smooth), lower.tail = FALSE)
  #PDt <- dt((min(SDM)-max(SDM))/mad(SDM),df=km-1)

 # print(c("low sigma pvalue:", PDt))
  #print(c("low sd pvalue:", dt((min(SSD)-max(SSD))/mad(SSD),df=km-1)))

  if(wn <= 0.15|(!all(table(ct)>min.cells))| PDt > 0.05){
    WNS <- "unclassified.prediction"
    print("low confidence in classification")
  }else {
    WNS <- ""
  }
    basel <- apply(norm.mat.smooth[, which(ct %in% which(SDM==min(SDM)))], 1, median)
    preN <- colnames(norm.mat.smooth)[which(ct %in% which(SDM==min(SDM)))]

  ### return both baseline and warning message
  RE <- list(basel, WNS, preN, ct)
  names(RE) <- c("basel", "WNS", "preN", "cl")
  return(RE)
    }


baseline.synthetic <- 
function(norm.mat=norm.mat, min.cells=10, n.cores){ 

 d <- parallelDist::parDist(t(norm.mat), threads = n.cores) ##use smooth and segmented data to detect intra-normal cells
  km <- 6
  fit <- hclust(d, method="ward.D2")
  ct <- cutree(fit, k=km)

  while(!all(table(ct)>min.cells)){
    km <- km -1
    ct <- cutree(fit, k=km)
    if(km==2){
      break
    }
  }

  
  expr.relat <- NULL
  syn <- NULL
  for(i in min(ct):max(ct)){
    data.c1 <- norm.mat[, which(ct==i)]
    sd1 <- apply(data.c1,1,sd)
    set.seed(123)
    syn.norm <- sapply(sd1,function(x)(x<- rnorm(1,mean = 0,sd=x)))
    relat1 <- data.c1 -syn.norm
    expr.relat <- rbind(expr.relat, t(relat1))
    syn <- cbind(syn,syn.norm)
    i <- i+1
  }

  reslt <- list(data.frame(t(expr.relat)), data.frame(syn), ct)
  names(reslt) <- c("expr.relat","syn.normal", "cl")
  
  return(reslt)
}


CNA.MCMC <- 
function(clu,fttmat, bins, cut.cor, n.cores){
  CON<- NULL
  for(i in min(clu):max(clu)){
    data.c <- apply(fttmat[, which(clu==i)],1, median)
    CON <- cbind(CON, data.c)
    i <- i+1
  }

  norm.mat.sm <- exp(CON)
  n <- nrow(norm.mat.sm)

  BR <- NULL

  for(c in 1:ncol(norm.mat.sm)){

    breks <- c(seq(1, as.integer(n/bins-1)*bins, bins),n)
    bre <- NULL

    for (i in 1:(length(breks)-2)){
      #i<-42
      a1<-  max(mean(norm.mat.sm[breks[i]:breks[i+1],c]), 0.001)
      posterior1 <-MCMCpack::MCpoissongamma(norm.mat.sm[breks[i]:breks[i+1],c], a1, 1, mc=1000)


      a2 <- max(mean(norm.mat.sm[(breks[i+1]+1):breks[i+2],c]), 0.001)
      posterior2 <-MCMCpack::MCpoissongamma(norm.mat.sm[(breks[i+1]+1):breks[i+2],c], a2, 1, mc=1000)

      if (ks.test(posterior1,posterior2)$statistic[[1]] > cut.cor){
        bre <- c(bre, breks[i+1])
       }

        i<- i+1
        }

    breks <- sort(unique(c(1, bre, n)))
    BR <- sort(unique(c(BR, breks)))
    c<-c+1
  }

  #print(paste(length(BR), " breakpoints", sep=""))

  ###CNA
  norm.mat.sm <- exp(fttmat)

  seg <- function(z){
      x<-numeric(n)
      for (i in 1:(length(BR)-1)){
        a<- max(mean(norm.mat.sm[BR[i]:BR[i+1],z]), 0.001)
        posterior1 <-MCMCpack::MCpoissongamma(norm.mat.sm[BR[i]:BR[i+1],z], a, 1, mc=1000)
        x[BR[i]:BR[i+1]]<-mean(posterior1)
         i<- i+1
      }
      x<-log(x)

  }

  seg.test <- parallel::mclapply(1:ncol(norm.mat.sm), seg, mc.cores = n.cores)
  logCNA <- matrix(unlist(seg.test), ncol = ncol(norm.mat.sm), byrow = FALSE)

  res <- list(logCNA, BR)
  names(res) <- c("logCNA","breaks")
  return(res)
}


convert.all.bins.hg20 <- 
function(DNA.mat, RNA.mat, n.cores){
##make list obj for each window
         DNA <- DNA.mat[-which(DNA.mat$chrom==24),]; dim(DNA)
         end <- DNA$chrompos
         start <- c(0, end[-length(end)])
          ls.all <- list()
          for(i in 1:nrow(DNA)){
          sub.anno <- full.anno[which(full.anno$chromosome_name==DNA$chrom[i]),]
          cent.gene <- 0.5*(sub.anno$start_position+sub.anno$end_position)
          x <- sub.anno$hgnc_symbol[which(cent.gene<=end[i] & cent.gene>= start[i])]
          if(length(x)==0){x <- "NA"}
          ls.all[[i]] <- x
          i<- i+1
          }

          ##convert gene to bin
          RNA <- RNA.mat[, 8:ncol(RNA.mat)]

         ###adj
          R.ADJ <- function(i){
            shr <- intersect(ls.all[[i]], RNA.mat$hgnc_symbol)
            if(length(shr)>0){
              Aj <- apply(RNA[which(RNA.mat$hgnc_symbol %in% shr), ], 2, median)
            }

          }

          test.adj <- parallel::mclapply(1:nrow(DNA), R.ADJ, mc.cores = n.cores)
          RNA.aj <- matrix(unlist(test.adj), ncol = ncol(RNA), byrow = TRUE)
          colnames(RNA.aj) <- colnames(RNA.mat)[8:ncol(RNA.mat)]

          mind <- which(test.adj=="NULL")

           if(length(mind)>1){
           ind <- 1:nrow(DNA)
           Rw <- ind[-which(ind %in% mind)]

           FK.again <- function(i){
            fkI <- abs(Rw-mind[i])
            fk <-  RNA.aj[which(fkI==min(fkI))[1], ]
          }

          tt.FK <-  parallel::mclapply(1:length(mind), FK.again, mc.cores = n.cores)
          FK <- matrix(unlist(tt.FK), ncol = ncol(RNA), byrow = TRUE)
          colnames(FK) <- colnames(RNA.mat)[8:ncol(RNA.mat)]

           RNA.aj <- cbind(Rw, RNA.aj)
           FK <- cbind(mind, FK)
           RNA.co <- data.frame(rbind(RNA.aj, FK))

          RNA.com <- RNA.co[order(RNA.co$Rw, decreasing = FALSE), ]
          RNA.adj <- cbind(DNA[, 1:3], RNA.com[, 2:ncol(RNA.com)])
         } else{
           RNA.adj <- cbind(DNA[, 1:3], RNA.aj)
         }

        reslt <- list(DNA, RNA.adj)
        names(reslt) <- c("DNA.adj", "RNA.adj")
        return(reslt)

        }


heatmap.3 <- 
function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      ColSideColorsSize = 1,
                      RowSideColorsSize = 1,
                      KeyValueName="Value",...){

    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
          return(TRUE)
      if (is.list(x))
          return(all(sapply(x, invalid)))
      else if (is.vector(x))
          return(all(is.na(x)))
      else return(FALSE)
    }

    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)

        if (!missing(ColSideColors)) {
           #if (!is.matrix(ColSideColors))
           #stop("'ColSideColors' must be a matrix")
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
                stop("'ColSideColors' must be a matrix of nrow(x) rows")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            #lhei <- c(lhei[1], 0.2, lhei[2])
             lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
        }

        if (!missing(RowSideColors)) {
            #if (!is.matrix(RowSideColors))
            #stop("'RowSideColors' must be a matrix")
            if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
                stop("'RowSideColors' must be a matrix of ncol(x) columns")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
            #lwid <- c(lwid[1], 0.2, lwid[2])
            lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }

    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)){
                par(mar = c(margins[1], 0, 0, 0.5))
                image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[,rowInd, drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
            if (length(rownames(RowSideColors)) > 0) {
                axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
            }
        }
    }

    if (!missing(ColSideColors)) {

        if (!is.matrix(ColSideColors)){
            par(mar = c(0.5, 0, 0, margins[2]))
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            par(mar = c(0.5, 0, 0, margins[2]))
            csc = ColSideColors[colInd, , drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            image(csc, col = as.vector(csc.colors), axes = FALSE)
            if (length(colnames(ColSideColors)) > 0) {
                axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
            }
        }
    }

    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
        cex.axis = cexCol)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
        cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }

        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, KeyValueName, line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}


annotateGenes.rn7 <- 
function(mat, ID.type="S", full.anno=full.anno.rn7){
  print("start annotation ...")

  if(substring(ID.type,1,1) %in% c("E", "e")){
    shar <- intersect(rownames(mat), full.anno$ensembl_gene_id)
    mat <- mat[which(rownames(mat) %in% shar),]
    anno <- full.anno[which(as.vector(full.anno$ensembl_gene_id) %in% shar),]
    anno <- anno[!duplicated(anno$mgi_symbol),]
    anno <- anno[order(match(anno$ensembl_gene_id, rownames(mat))),]
    data <- cbind(anno, mat)

  }else if(substring(ID.type,1,1) %in% c("S", "s")) {

    shar <- intersect(rownames(mat), full.anno$mgi_symbol)
    rownames(mat)[1:10]
    full.anno$mgi_symbol[]

    mat <- mat[which(rownames(mat) %in% shar),]
    anno <- full.anno[which(as.vector(full.anno$mgi_symbol) %in% shar),]
    anno <- anno[!duplicated(anno$mgi_symbol),]
    anno <- anno[order(match(anno$mgi_symbol, rownames(mat))),]
    data <- cbind(anno, mat)
  }
}

