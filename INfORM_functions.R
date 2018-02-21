
suppressMessages(library(igraph))
suppressMessages(library(TopKLists))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(GO.db))
suppressMessages(library(plyr))
suppressMessages(library(GSEABase))
suppressMessages(library(GOSemSim))
suppressMessages(library(treemap))
suppressMessages(library(abind))
suppressMessages(library(minet))
suppressMessages(library(foreach))
suppressMessages(library(parallel))
suppressMessages(library(doParallel))
suppressMessages(library(ggplot2))
suppressMessages(library(visNetwork))
suppressMessages(library(radarchart))
suppressMessages(library(WriteXLS))
suppressMessages(library(gplots))

net_attr <- c("betweenness", "cc", "degree", "eccentricity", "closeness", "eigenvector")

methods <- c("clr","aracne","mrnet","mrnetb") # i
est.opt <- c("pearson","spearman","kendall","mi.empirical","mi.mm","mi.shrink","mi.sg")# j
disc.opt <- c("none","equalfreq","equalwidth","globalequalwidth") # k
shiny_app_vars <- NULL
shiny_app_env <- NULL

combineList <- function(...){
	myParam <- list(...)
	#c(myParam)
	myParam
}

utils::globalVariables(names=c("GO.db"))

#Get lock to avoid current access by different threads 
getLock <- function(tmpDir, logFile){
        lockDir <- paste0(tmpDir, "/", gsub(".txt", "", logFile), "_lock/")
        lockChk <- FALSE
        while(lockChk==FALSE){
                Sys.sleep(1)
                lockChk <- dir.create(lockDir, showWarnings=FALSE)
        }
}

#Delete the lock directory and remove lock 
unLock <- function(tmpDir, logFile){
        lockDir <- paste0(tmpDir, "/", gsub(".txt", "", logFile), "_lock/")
        unlink(c(lockDir), recursive=TRUE)
}

#' Infer correlation matrix from gene expression table by using mutual information.
#'
#' calculate_correlation_matrix uses the MINET package to create correlation matrix by mutual information method. User can specify
#' the inference algorithms, correlation calculation methods and discretization methods to create mutiple combinations of parameters 
#' for multiple runs of minet(). Multiple inferences are unified by taking median to create a consensus matrix.
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom minet minet
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach getDoParWorkers getDoParName %dopar%
#' @importFrom plyr aaply laply
#' @importFrom stats quantile median
#' @importFrom utils capture.output
#'
#' @param gx_table Gene expression table as a data frame.
#' @param iMethods Vector of valid inference algorithms for MINET package.
#' @param iEst Vector of valid correlation methods for MINET package.
#' @param iDisc Vector of valid discretization methods for MINET package.
#' @param summ_by Measure for summarizing the co-expression score of gene pairs from a set of netword to create a consensus network; options:median, mean, max; default:median.
#' @param ncores Number of cores for running instances of MINET in parallel default:2.
#' @param debug_output Print help and status messages to help debug the running of the function default:FALSE.
#' @param updateProgress Shiny application can request for update of progress from this function default:NULL.
#' @return A binary symmetrix matrix representing the median of mutual information correlation computed across various MINET combinations
#' @examples
#' \dontrun{
#' calculate_correlation_matrix(gx_table=gene_expression.df,
#' iMethods=c("clr","aracne","mrnet","mrnetb"),
#' iEst=c("pearson","spearman","kendall","mi.empirical","mi.mm","mi.shrink","mi.sg"),
#' iDisc=c("none","equalfreq","equalwidth","globalequalwidth"),
#' summ_by="median",
#' ncores=2,
#' debug_output=TRUE,
#' updateProgress=NULL)
#' }
#' @keywords internal
#' @export
calculate_correlation_matrix <- function(gx_table, iMethods, iEst, iDisc, summ_by="median", ncores=2, debug_output=FALSE, updateProgress=NULL){
	parList <- list()
	out.GX.list<- list()
	out.tmp.list<- list()
	out.list<- list()
	tGX <- t(gx_table)
	stdGX <- scale(tGX)

	if(is.null(iMethods)){
		stop("Please Select atLeast One Inference Algorithm!")
	}

	if(is.null(iEst)){
		stop("Please Select At Least One Correlation!")
	}

	if(is.null(iDisc)){
		stop("Please Select At Least One Discretization Method!")
	}

	cntr <- 0
	for(i in 1:length(iMethods)){
		#print(methods2[i])
		for(j in 1:length(iEst)){
		#print(paste("-",est.opt2[j]))
			for(k in 1:length(iDisc)){
                                if(grepl("mi\\..*", est.opt[j]) && disc.opt[k]=="none"){
                                        next
                                }
				cntr <- cntr + 1
				parList[[cntr]] <- list()
				parList[[cntr]][["mt"]] <- iMethods[i]
				parList[[cntr]][["est"]] <- iEst[j]
				parList[[cntr]][["disc"]] <- iDisc[k]
			}
		}
	}
	#print(paste0("List of MINET Combinations: ", parList))

	if(ncores > parallel::detectCores()){
		ncores <- parallel::detectCores()
	}

	cl <- parallel::makeCluster(ncores, outfile="")
	doParallel::registerDoParallel(cl)
	print(paste("Total Cores: ", parallel::detectCores(), sep=""))
	print(paste("Total Workers: ", foreach::getDoParWorkers(), sep=""))
	print(paste("DoPar Name: ",  foreach::getDoParName(), sep=""))

	print(paste0("Before For Each, ", "Number of Combinations:", length(parList)))
	#utils::capture.output(print(paste0("Starting Parallel Cluster For ", length(parList), " Combinations ", timestamp())), file="minet-log.txt", append=TRUE)
	#utils::capture.output(print(paste0("Starting Parallel Cluster For ", length(parList), " Combinations ", timestamp())), file="minet-to-run.txt", append=TRUE)
	#utils::capture.output(print(paste0("Starting Parallel Cluster For ", length(parList), " Combinations ", timestamp())), file="minet-completed.txt", append=TRUE)
	utils::capture.output(print(paste0("Starting Parallel Cluster For ", length(parList), " Combinations ", timestamp())), file="minet-log.txt")
	utils::capture.output(print(paste0("Starting Parallel Cluster For ", length(parList), " Combinations ", timestamp())), file="minet-to-run.txt")
	utils::capture.output(print(paste0("Starting Parallel Cluster For ", length(parList), " Combinations ", timestamp())), file="minet-completed.txt")
	#out.tmp.list <- foreach(i=1:length(parList)) %do% {

        #tmpDir <- tempdir()
        #if(is.null(env)){
        #    env <- environment()
        #}
        logDir <- "./"
        env <- environment()
	parallel::clusterExport(cl, list("getLock", "unLock", "logDir"), envir=env)

	out.tmp.list <- foreach::foreach(i=1:length(parList)) %dopar% {
		#capture.output(print("Start For Each"), file="minet-log.txt", append=T)
		mt <- parList[[i]]$mt
		est <- parList[[i]]$est
		disc <- parList[[i]]$disc
		#print(parList[[i]]$mt)

		#if((est == "mi.empirical" | est == "mi.mm" | est == "mi.shrink" | est == "mi.sg") & disc == "none") {
		#	miMat <- -1
		#	miMatName <- "np"
		#}
		#else{
                if(debug_output==TRUE){
                        getLock(tmpDir=logDir, logFile="minet-log.txt")
                        utils::capture.output(print(paste("----",mt,est,disc, sep="__")), file="minet-log.txt", append=TRUE)
                        unLock(tmpDir=logDir, logFile="minet-log.txt")
                }

                getLock(tmpDir=logDir, logFile="minet-to-run.txt")
                utils::capture.output(print(paste0("Iteration-", i, ": ", mt, "-", est, "-", disc)), file="minet-to-run.txt", append=TRUE)
                unLock(tmpDir=logDir, logFile="minet-to-run.txt")
                #print("Before Updating text")
                #if (is.function(updateProgress)){
                #	text <- paste("MINET: ", mt, "-", est, "-", disc, "-", sep="")
                #	print("Updating text")
                #	updateProgress(detail = text)
                #}

                ptm <- proc.time()
                miMat <- minet::minet(stdGX, method=mt, estimator=est, disc=disc)

                getLock(tmpDir=logDir, logFile="minet-completed.txt")
                utils::capture.output(print(paste0("Iteration-", i, ", ", mt, "-", est, "-", disc, ": ", "MINET Execution Time - ", round(proc.time() - ptm)[3], " sec")), file="minet-completed.txt", append=TRUE)
                #capture.output(print(proc.time() - ptm), file="minet-log.txt", append=T)
                unLock(tmpDir=logDir, logFile="minet-completed.txt")


                miMatName <- paste(mt,est,disc,sep="__")
		#}
		out.list <- list("mat"=miMat, "name"=miMatName)
	}
	utils::capture.output(print("For Each Finished, Stopping Cluster..."), file="minet-log.txt", append=TRUE)
	parallel::stopCluster(cl)

	for(i in 1:length(out.tmp.list)){
	  out.GX.list[[i]] <- out.tmp.list[[i]]$mat
	  names(out.GX.list)[[i]] <- out.tmp.list[[i]]$name
	}

	#if(debug_output==TRUE)
	#print(paste("Got minet list of lists --- ", str(out.GX.list)))

        #removing ones having no values
        #llGX <- out.GX.list[-(which(names(out.GX.list)=="np"))]
        idx <- which(names(out.GX.list)=="np")
        if(length(idx)>0){
            llGX <- out.GX.list[-idx]
        }else{
            llGX <- out.GX.list
        }

        if(length(llGX)>1){
            print(paste0("Making ", summ_by, " Matrix..."))
            arr1<-abind::abind(llGX,along=3)
            if(summ_by=="median"){
                matGX <- apply(arr1,c(1,2),median)
            }else if(summ_by=="mean"){
                matGX <- apply(arr1,c(1,2),mean)
            }else if(summ_by=="max"){
                matGX <- apply(arr1,c(1,2),max)
            }
        }else{
            cat("Only one matrix, not computing ", summ_by, " matrix!\n")
            matGX <- llGX[[1]]
        }

        print("Return Matrix")
        return(matGX)
}

#' Create edge ranked inference matrix by combining information from different inference algorithms with the help of internal function calculate_correlation_matrix().
#'
#' get_ranked_consensus_matrix uses the internal function calculate_correlation_matrix() to get a single consensus matrix per inference algorithm. User can specify
#' the inference algorithms, correlation calculation methods and discretization methods, a combination of parameters will be created per inference algorithm to run calculate_correlation_matrix(), 
#' this will generate a consensus matrix per inference algorithm. The consesus matrices from different inference algorithms are used to create a single binary matrix by rank based selection of edges.
#'
#' @importFrom TopKLists Borda
#'
#' @param gx_table Gene expression table as a data frame.
#' @param iMethods Vector of valid inference algorithms for MINET package.
#' @param iEst Vector of valid correlation methods for MINET package.
#' @param iDisc Vector of valid discretization methods for MINET package.
#' @param ncores Number of cores for running instances of MINET in parallel; default:2.
#' @param matList List of co-expression matrices provided by user; default:NULL.
#' @param mat_weights Type of scores in the user provided matrices; default:rank.
#' @param ensemble_strategy Strategy to use for ensemble "minet" for minet generated matrices, or "user" for only user provided matrices, or "minet+user" to combine minet generated matrices and user provided matrices; default:minet.
#' @param debug_output Print help and status messages to help debug the running of the function; default:FALSE.
#' @param updateProgress Shiny application can request for update of progress from this function; default:NULL.
#' @return A symmetrix matrix with edge ranks representing the edge rank based consensus from different inference algorithms.
#' @examples
#' \dontrun{
#' get_ranked_consensus_matrix(gx_table=gene_expression.df,
#' iMethods=c("clr","aracne","mrnet","mrnetb"),
#' iEst=c("pearson","spearman","kendall","mi.empirical","mi.mm","mi.shrink","mi.sg"),
#' iDisc=c("none","equalfreq","equalwidth","globalequalwidth"),
#' ncores=12,
#' matList=list(mat1, mat2, matN),
#' mat_weights="rank",
#' ensemble_strategy="minet",
#' debug_output=TRUE,
#' updateProgress=TRUE)
#' }
#' @keywords internal
#' @export
#get_ranked_consensus_binary_matrix <- function(gx_table, iMethods, iEst, iDisc, ncores=2, debug_output=FALSE, updateProgress=NULL){
get_ranked_consensus_matrix <- function(gx_table=NULL, iMethods=NULL, iEst=NULL, iDisc=NULL, summ_by="median", score_type="median", ncores=2, matList=NULL, mat_weights="rank", ensemble_strategy="minet", debug_output=FALSE, updateProgress=NULL){
	mat_ll <- list()
    	ranked_edges_ll <- list()

        if(!score_type %in% c("mean", "median", "geo.mean", "l2norm")){
                print(paste0("score_type : '", score_type, "' is not a valid Borda score!"))
                return(NULL)
        }

        if(length(grep("minet" ,ensemble_strategy))>0){
                mthdCount <- 1
                totalMthds <- length(iMethods)+1
                for(mthd in iMethods){
                        if (is.function(updateProgress)) {
                                text <- paste0("'", mthd, "' consensus by ", summ_by)
                                value <- mthdCount / totalMthds
                                updateProgress(detail = text, value = value)
                        }
                        mthdCount <- mthdCount + 1

                        print(paste0("Calculate correlation matrix for method : ", mthd))
                        mat_ll[[mthd]] <- calculate_correlation_matrix(gx_table=gx_table, iMethods=mthd, iEst=iEst, iDisc=iDisc, summ_by=summ_by, ncores=ncores)

                        print(paste0("Get ranked edges for method : ", mthd))
                        mat_ll[[mthd]][lower.tri(mat_ll[[mthd]], diag=TRUE)] <- NA
                                
                        edge_df <- as.data.frame(as.table(mat_ll[[mthd]]))
                        print("Minet edge_df before:")
                        print(str(edge_df))
                        edge_df <- edge_df[-which(is.na(edge_df$Freq)),]
                        edge_df <- data.frame(edge=paste0(edge_df$Var1,";",edge_df$Var2), weight=edge_df$Freq, stringsAsFactors=FALSE)
                        print("Minet edge_df after:")
                        print(str(edge_df))

                        ranked_edges_ll[[mthd]] <- edge_df[order(edge_df$weight, decreasing=TRUE), "edge"]
                }
        }

        if(length(grep("user", ensemble_strategy))>0){
                matList <- as.list(matList)
                for(i in c(1:length(matList))){
                        itmName <- paste0("user_", i)
                        matList[[i]][lower.tri(matList[[i]], diag=TRUE)] <- NA
                                
                        edge_df <- as.data.frame(as.table(matList[[i]]))
                        print("User edge_df before:")
                        print(str(edge_df))
                        edge_df <- edge_df[-which(is.na(edge_df$Freq)),]
                        edge_df <- data.frame(edge=paste0(edge_df$Var1,";",edge_df$Var2), weight=edge_df$Freq, stringsAsFactors=FALSE)
                        hasNulls <- FALSE
                        if(mat_weights=="rank"){
                                nullIdx <- which(edge_df$weight==0)
                                if(length(nullIdx)>0){
                                        edge_df_null <- edge_df[nullIdx,]
                                        edge_df <- edge_df[-nullIdx,]
                                        hasNulls <- TRUE
                                }
                        }
                        print("User edge_df after:")
                        print(str(edge_df))

                        ordOption <- ifelse(mat_weights=="rank", FALSE, TRUE)
                        print(paste0("ordOption : ", ordOption))
                        edge_df <- edge_df[order(edge_df$weight,decreasing=ordOption),]
                        if(hasNulls){
                                edge_df <- as.data.frame(rbind(edge_df, edge_df_null), stringsAsFactors=FALSE)
                        }
                        #ranked_edges_ll[[itmName]] <- edge_df[order(edge_df$weight, decreasing=ordOption), "edge"]
                        ranked_edges_ll[[itmName]] <- edge_df$edge
                }
        }
        print("ranked_edges_ll:")
        print(length(ranked_edges_ll))
        print(names(ranked_edges_ll))
        print(lapply(ranked_edges_ll, length))

        if(is.function(updateProgress)) {
                updateProgress(detail = "Consensus Binary", value = 1)
        }

	print("Perform Borda on list of list of ranked edges.")
	borda_res <- TopKLists::Borda(ranked_edges_ll)
        #Borda.plot(borda_res)

	print("Get a consensus binary matrix by selecting the most significant ranked edges from median rank of Borda result.")
        if(length(mat_ll)>0){
	        rank_mat <- mat_ll[[1]]
        }else{
                rank_mat <- matList[[1]]
        }
	rank_mat[,] <- 0
        print("rank_mat")
        print(dim(rank_mat))

	ra_list <- borda_res$TopK[[score_type]]
        #Kendall.plot(ranked_edges_ll, median_list)

	#if(edge_selection_strategy=="default"){
	#	#input_genes <- dim(gx_table)[1]
	#	input_genes <- nrow(rank_mat)
	#	genes <- NULL
	#	total_genes <- 0
	#	cutoffIdx <- NULL
	#	for(i in c(1:length(median_list))){
	#		if(total_genes<input_genes){
	#			local_genes <- strsplit(median_list[i], ";")[[1]]
	#			rank_mat[local_genes[1],local_genes[2]] <- i
	#			rank_mat[local_genes[2],local_genes[1]] <- i
	#			genes[local_genes[1]] <- 1
	#			genes[local_genes[2]] <- 1
	#			total_genes <- length(genes)
	#		}else{
	#			cutoffIdx <- i-1
	#			break
	#		}
	#	}
	#}else if(edge_selection_strategy=="top"){
	#	cutOff <- round((as.numeric(topN)*length(median_list))/100)
	#	for(i in c(1:length(median_list))){
	#		if(i<=cutOff){
	#			local_genes <- strsplit(median_list[i], ";")[[1]]
	#			rank_mat[local_genes[1],local_genes[2]] <- i
	#			rank_mat[local_genes[2],local_genes[1]] <- i
	#		}else{
	#			break
	#		}
	#	}
	#}

        print("ra_list:")
        print(str(ra_list))
        print(length(ra_list))
        print(head(ra_list))
        for(i in c(1:length(ra_list))){
                local_genes <- strsplit(ra_list[i], ";")[[1]]
                rank_mat[local_genes[1],local_genes[2]] <- i
                rank_mat[local_genes[2],local_genes[1]] <- i
        }

	print("Rank matrix computed, returning!")
	return(rank_mat)
}

#' Get edge rank list and binary inference matrix from edge rank matrix computed by get_ranked_consensus_matrix().
#'
#' parse_edge_rank_matrix parses the edge rank matrix created by using the internal function get_ranked_consensus_matrix_matrix() to get a ranked edge list and a binary matrix.
#'
#' @importFrom TopKLists Borda
#'
#' @param edge_rank_matrix A symmetrix matrix with edge ranks as weight.
#' @param edge_selection_strategy Strategy to select edges for ensemble, "default" selects ranked edges untill all nodes have degree>=1, "top" swtiches to top N percentage ranked genes specfied by topN parameter.
#' @param mat_weights Type of scores in the ranked matrix; default:rank.
#' @param topN Top N percentage ranked edges to create ensemble if the edge_selection_strategy is "top"; default:10
#' @param debug_output Print help and status messages to help debug the running of the function; default:FALSE.
#' @param updateProgress Shiny application can request for update of progress from this function; default:NULL.
#' @return A list containing a vector of consensus edge ranks and a binary symmetrix matrix representing the edge rank matrix.
#' @examples
#' \dontrun{
#' parse_edge_rank_matrix <- function(edge_rank_matrix, 
#' edge_selection_strategy="default", 
#' mat_weights="rank", 
#' topN=10, 
#' debug_output=FALSE, 
#' updateProgress=NULL)
#' }
#' @keywords internal
#' @export
parse_edge_rank_matrix <- function(edge_rank_matrix, edge_selection_strategy="default", mat_weights="rank", topN=10, debug_output=FALSE, updateProgress=NULL){
        bin_mat <- edge_rank_matrix
        #idx <- which(bin_mat>0)
        #if(length(idx)>0){
        #        print("Creating binary matrix...")
        #        bin_mat[which(bin_mat>0)] <- 1
        #}
        bin_mat[,] <- 0

        print("Getting edge list ordered by rank...")
        rank_matrix <- edge_rank_matrix
        rank_matrix[lower.tri(rank_matrix, diag=TRUE)] <- NA
        edge_df <- as.data.frame(as.table(rank_matrix))
        edge_df <- edge_df[-which(is.na(edge_df$Freq)),]
        edge_df <- data.frame(edge=paste0(edge_df$Var1,";",edge_df$Var2), weight=edge_df$Freq, stringsAsFactors=FALSE)
        edge_df <- edge_df[which(edge_df$weight>0),]
        ordOption <- ifelse(mat_weights=="rank", FALSE, TRUE)
        print(paste0("ordOption : ", ordOption))
        print(class(edge_df$weight))
        edge_rank <- edge_df$edge[order(edge_df$weight, decreasing=ordOption)]
        print(paste0("edge rank before:", length(edge_rank)))
        print(head(edge_rank))

	if(edge_selection_strategy=="default"){
		#input_genes <- nrow(rank_matrix)
                input_genes <- length(unique(unlist(strsplit(edge_rank, ";"))))
                print(paste0("input genes : ", input_genes))
		genes <- NULL
		total_genes <- 0
		cutoffIdx <- NULL
		for(i in c(1:length(edge_rank))){
			if(total_genes<input_genes){
				local_genes <- strsplit(edge_rank[i], ";")[[1]]
				bin_mat[local_genes[1],local_genes[2]] <- 1
				bin_mat[local_genes[2],local_genes[1]] <- 1
				genes[local_genes[1]] <- 1
				genes[local_genes[2]] <- 1
				total_genes <- length(genes)
				cutoffIdx <- i
			}else{
                                print(paste0("Cutoff before break:", cutoffIdx))
				break
			}
		}
	}else if(edge_selection_strategy=="top"){
		cutoffIdx <- round((as.numeric(topN)*length(edge_rank))/100)
		for(i in c(1:length(edge_rank))){
			if(i<=cutoffIdx){
				local_genes <- strsplit(edge_rank[i], ";")[[1]]
				bin_mat[local_genes[1],local_genes[2]] <- 1
				bin_mat[local_genes[2],local_genes[1]] <- 1
			}else{
				break
			}
		}
	}
        print(paste0("Cutoff:", cutoffIdx))
        edge_rank <- edge_rank[c(1:cutoffIdx)]
        print(paste0("edge rank after:", length(edge_rank)))

        res_ll <- list(bin_mat=bin_mat, edge_rank=edge_rank)

	print("Returning!")
	return(res_ll)
}

#' Create iGraph object from a symmetrical adjacency matrix, annotate it with centrality attributes and return the annotated iGraph.
#'
#' @importFrom igraph graph.adjacency vertex_attr betweenness closeness degree eigen_centrality transitivity edge_attr list.vertex.attributes
#'
#' @param adj_mat Binary consensus matrix computed by using get_ranked_consensus_binary_matrix() function.
#' @return Annotated igraph object representing the binary consensus matrix computed by get_ranked_consensus_binary_matrix() function.
#' @examples
#' \dontrun{
#' get_iGraph(adj_mat=binary.inference.matrix)
#' }
#' @keywords internal
#' @export
get_iGraph <- function(adj_mat){
	if(is.null(adj_mat)){
                stop("Input Adjacency Matrix is NULL!")
        }

        if(!isSymmetric(adj_mat)){
                print("Matrix is not symmetric. Getting symmetric matrix!")
                #adj_mat <- get_symmetric_matrix(adj_mat)
                adj_mat <- pmax(adj_mat, t(adj_mat))
        }

        iG <- igraph::graph.adjacency(adj_mat, mode="undirected", weighted=NULL)

        igraph::vertex_attr(iG, name="betweenness") <- as.vector(igraph::betweenness(iG))
        igraph::vertex_attr(iG, name="closeness") <- as.vector(igraph::closeness(iG, normalized=TRUE))
        igraph::vertex_attr(iG, name="degree") <- as.vector(igraph::degree(iG))
        igraph::vertex_attr(iG, name="eigenvector") <- as.vector(igraph::eigen_centrality(iG)$vector)
        igraph::vertex_attr(iG, name="cc") <- igraph::transitivity(iG, type="local", isolates="zero")

        igraph::vertex_attr(iG, name="color") <- "lightgray"
        igraph::vertex_attr(iG, name="highlightcolor") <- "darkgray"
        igraph::edge_attr(iG, name="color") <- "lightgray"
        igraph::edge_attr(iG, name="highlightcolor") <- "darkgray"

        print(igraph::list.vertex.attributes(iG))
        iG
}

#' Annotate iGraph object with centrality attributes and return the annotated iGraph.
#'
#' @importFrom igraph vertex_attr betweenness closeness degree eigen_centrality transitivity edge_attr list.vertex.attributes list.edge.attributes
#'
#' @param iG igraph object created from the binary consensus matrix computed by using get_ranked_consensus_binary_matrix() function.
#' @return Annotated igraph object representing the binary consensus matrix computed by get_ranked_consensus_binary_matrix() function.
#' @examples
#' \dontrun{
#' annotate_iGraph(iG=inferred.iGraph)
#' }
#' @keywords internal
#' @export
annotate_iGraph <- function(iG){
        igraph::vertex_attr(iG, name="betweenness") <- as.vector(igraph::betweenness(iG))
        igraph::vertex_attr(iG, name="closeness") <- as.vector(igraph::closeness(iG, normalized=TRUE))
        igraph::vertex_attr(iG, name="degree") <- as.vector(igraph::degree(iG))
        igraph::vertex_attr(iG, name="eigenvector") <- as.vector(igraph::eigen_centrality(iG)$vector)
        igraph::vertex_attr(iG, name="cc") <- igraph::transitivity(iG, type="local", isolates="zero")

        vertex_attr_list <- igraph::list.vertex.attributes(iG)
        edge_attr_list <- igraph::list.edge.attributes(iG)

        if(!("color" %in% vertex_attr_list)){
                igraph::vertex_attr(iIG, name="color") <- "lightgray"
        }
        if(!("highlightcolor" %in% vertex_attr_list)){
                igraph::vertex_attr(iIG, name="highlightcolor") <- "darkgray"
        }

        if(!("color" %in% edge_attr_list)){
                igraph::edge_attr(iIG, name="color") <- "lightgray"
        }
        if(!("highlightcolor" %in% edge_attr_list)){
                igraph::edge_attr(iIG, name="highlightcolor") <- "darkgray"
        }

        print(igraph::list.vertex.attributes(iG))
        return(iG)
}

#' Set general vertex colors and also colors to highlight the identified important vertices.
#'
#' @importFrom igraph vertex_attr V
#' @importFrom stats quantile
#'
#' @param iGraph igraph object created from the binary consensus matrix computed by using get_ranked_consensus_binary_matrix() function.
#' @param gx_data_table Gene expression table as a data frame, must have gene symbols as rownames and sample names as colnames.
#' @param dgx_table Differential Gene expression table as a data frame, must have gene symbols as rownames but should not have colnames.
#' @param pos_cor_color Valid color name or hex code for genes positively associated with differential gene expression. 
#' @param pos_cor_highlight_color Valid color name or hex code for highlighted genes positively associated with differential gene expression.
#' @param neg_cor_color Valid color name or hex code for genes negatively associated with differential gene expression.
#' @param neg_cor_highlight_color Valid color name or hex code for highlighted genes negatively associated with differential gene expression.
#' @param pos_perc Percentile cutoff for positive association of genes with differential gene expression score, default:0.95.
#' @param neg_perc Percentile cutoff for negative association of genes with differential gene expression score, default:0.05.
#' @return igraph object with user provided vertex color for normal and highlighted vertices.
#' @examples
#' \dontrun{
#' set_vertex_color(iGraph=inferred.igraph,
#' gx_data_table=gene_expression.df,
#' dgx_table=differential_gene_expression.df,
#' pos_cor_color=color.pos, pos_cor_highlight_color=color.pos.highlight,
#' neg_cor_color=color.neg,
#' neg_cor_highlight_color=color.neg.highlight,
#' pos_perc=pos.threshold,
#' neg_perc=neg.threshold
#' )
#' }
#' @keywords internal
#' @export
set_vertex_color <- function(iGraph, gx_data_table, dgx_table, pos_cor_color="salmon", pos_cor_highlight_color="red", neg_cor_color="lightblue", neg_cor_highlight_color="royalblue", pos_perc=0.95, neg_perc=0.05){
	#Intitialize color vectors
        col_length <- vcount(iGraph)
        color_vector <- rep("lightgrey", col_length)
        highlight_color_vector <- rep("darkgrey", col_length)
	vNames <- igraph::V(iGraph)$name
	if("score" %in% list.vertex.attributes(iGraph))
	{
		posCor <- vNames[which(igraph::V(iGraph)$score>=stats::quantile(igraph::V(iGraph)$score, pos_perc))]
		negCor <- vNames[which(igraph::V(iGraph)$score<=stats::quantile(igraph::V(iGraph)$score, neg_perc))]
	}
	else{
		dgx_table <- as.matrix(dgx_table)
		dgx_table <- dgx_table[vNames,]
		scoreColIdx <- NULL
		if("score" %in% colnames(dgx_table)){
			scoreColIdx <- which(colnames(dgx_table)=="score")
		}else if(any(c("lfc", "logfc", "log.fc", "log-fc") %in% tolower(colnames(dgx_table)))){
			scoreColIdx <- which(tolower(colnames(dgx_table))%in%c("lfc", "logfc", "log.fc", "log-fc"))
		}else{
			if(ncol(dgx_table)>=2){
				scoreColIdx <- 2
			}else{
				scoreColIdx <- 1
			}
		}
		posCor <- rownames(dgx_table)[which(dgx_table[,scoreColIdx]>= stats::quantile(as.vector(dgx_table[,scoreColIdx]), pos_perc))]
		negCor <- rownames(dgx_table)[which(dgx_table[,scoreColIdx]<= stats::quantile(as.vector(dgx_table[,scoreColIdx]), neg_perc))]
	}

        print("For Loop for populating vertex color vector...")
        #for(i in c(which(is.element(rownames(gx_data_table),posCor)))){
        for(i in c(which(is.element(vNames,posCor)))){
                color_vector[i] <- pos_cor_color
                highlight_color_vector[i] <- pos_cor_highlight_color
        }

        #for(i in c(which(is.element(rownames(gx_data_table),negCor)))){
        for(i in c(which(is.element(vNames,negCor)))){
                color_vector[i] <- neg_cor_color
                highlight_color_vector[i] <- neg_cor_highlight_color
        }

        igraph::vertex_attr(iGraph, name="color") <- color_vector
        igraph::vertex_attr(iGraph, name="highlightcolor") <- highlight_color_vector
        return(iGraph)
}

#' Set general edge colors and also colors to highlight the identified important edges.
#'
#' set_edge_color assigns color and highlight color to the edges. User provided gene expression table is used to compute correlation 
#' between the genes, edges with positive value of correlation are assigned specific colors and edges with negative value of correlation 
#' are assigned specific colors. The highlight colors are used to faciliate sub-network extraction and highlighting the seed genes used for 
#' extracting the sub-networks.
#'
#' @importFrom igraph get.edges edge_attr E
#' @importFrom stats cor
#'
#' @param iGraph igraph object created from the binary consensus matrix computed by using get_ranked_consensus_binary_matrix() function.
#' @param gx_data_table Gene expression table as a data frame, must have gene symbols as rownames and sample names as colnames.
#' @param pos_cor_color Valid color name or hex code for edges showing positive correlation between genes. 
#' @param pos_cor_highlight_color Valid color name or hex code for highlighted edges showing positive correlation between genes.
#' @param neg_cor_color Valid color name or hex code for edges showing negative correlation between genes.
#' @param neg_cor_highlight_color Valid color name or hex code for highlighted edges showing negative correlation between genes.
#' @return igraph object with user provided vertex color for normal and highlighted vertices.
#' @examples
#' \dontrun{
#' set_edge_color(iGraph=inferred.igraph,
#' gx_data_table=gene_expression.df,
#' pos_cor_color=color.pos,
#' pos_cor_highlight_color=color.pos.highlight,
#' neg_cor_color=color.neg,
#' neg_cor_highlight_color=color.neg.highlight
#' )
#' }
#' @keywords internal
#' @export
set_edge_color <- function(iGraph, gx_data_table, pos_cor_color="salmon", pos_cor_highlight_color="red", neg_cor_color="lightblue", neg_cor_highlight_color="royalblue"){
	tGX <- t(gx_data_table)
	corGX <- cor(tGX)

	matColor <- ifelse(corGX<=0, neg_cor_color, pos_cor_color)
	matHighlightColor <- ifelse(corGX<=0, neg_cor_highlight_color, pos_cor_highlight_color)

	#edgeIDs <- igraph::get.edges(iGraph, igraph::E(iGraph))
	edgeIDs <- igraph::ends(iGraph, igraph::E(iGraph), names=TRUE)
        #print("edgeIDs")
        #print(edgeIDs)
        cols <- apply(edgeIDs, 1, function(x) {
		matColor[x[1],x[2]]
	})
        #print(cols)

	igraph::edge_attr(iGraph, name="color") <- apply(edgeIDs, 1, function(x) {
		matColor[x[1],x[2]]
	})

	igraph::edge_attr(iGraph, name="highlightcolor") <- apply(edgeIDs, 1, function(x) {
		matHighlightColor[x[1],x[2]]
	})
	return(iGraph)
}

#' Get top ranked candidates computed by Borda on the provided list of attributes and the cutoff.
#'
#' get_ranked_gene_list uses the annotation associated with each vertex to rank them and generated a list of vertices names ordered by rank. 
#' By default the annotations calculated by INfORM are "betweenness", "cc", "degree", "eccentricity", "closeness" & "eigenvector". The vertices are 
#' ranked by each annotation separately nad then the ranks are unified by means of Borda(). User must provide the annotations to use in ranking scheme, 
#' if the user has custom annotations such as "score" for differential gene expression then it can also be used for ranking.
#'
#' @importFrom igraph get.vertex.attribute
#' @importFrom TopKLists Borda
#' @importFrom utils head
#'
#' @param iGraph igraph object created from the binary consensus matrix computed by using get_ranked_consensus_binary_matrix() function.
#' @param rank_list_attr Vector of network attributes/scores to use for generating a combined gene rank
#' @param debug_output Print help and status messages to help debug the running of the function default:FALSE.
#' @return vector of genes ordered by rank, based on the selected attributes associated to each gene.
#' @examples
#' \dontrun{
#' get_ranked_gene_list(iGraph=inferred.igraph,
#' rank_list_attr=c("betweenness",
#' "cc",
#' "degree",
#' "eccentricity",
#' "closeness",
#' "eigenvector",
#' "score"),
#' debug_output=FALSE
#' )
#' }
#' @keywords internal
#' @export
get_ranked_gene_list <- function(iGraph, rank_list_attr, debug_output=FALSE){
        if(is.null(rank_list_attr)){
                stop("List of attributes for ranking is NULL!")
        }
        
        attrOrdMat <- list()
        for(a in 1:length(rank_list_attr)){
		val_ord=TRUE

                if(debug_output==TRUE)
                print(paste("Using attribute: ", rank_list_attr[a], sep=""))

                attrValueVector <- igraph::get.vertex.attribute(iGraph, rank_list_attr[a])

                if(rank_list_attr[a] == "score")
                attrValueVector <- abs(attrValueVector)

                attrOrdList <- cbind(igraph::V(iGraph)$name, attrValueVector)[order(attrValueVector, decreasing=val_ord)]

                if(debug_output==TRUE)
                print(utils::head(attrOrdList))

                attrOrdMat[[a]] <- attrOrdList
        }

        attrBorda <- TopKLists::Borda(attrOrdMat)

        if(debug_output==TRUE)
        print(utils::head(attrBorda$TopK, n=10))

	return(attrBorda$TopK$median)
}

#' Get subgraph from the main igraph.
#'
#' Extract a subgraph from the main graph by using a list of genes to identify and select genes interconnected in a hub with them.
#'
#' @importFrom igraph neighborhood V induced.subgraph
#'
#' @param iGraph igraph object representing the main graph from which the subgraph should be extracted.
#' @param selected_vertices Vector of genes to be used for identifying interconnected hub.
#' @param vertex_level define the area of inclusion around the selected_vertices default:1.
#' @param shortest_paths Use shortest paths between selected_vertices to identify interconnected hubs default:FALSE.
#' @return igraph object representing the extracted subgraph
#' @examples
#' \dontrun{
#' get_sub_graph(iGraph=main.graph, selected_vertices, vertex_level=1, shortest_paths=TRUE)
#' }
#' @keywords internal
#' @export
get_sub_graph <- function(iGraph, selected_vertices, vertex_level=1, shortest_paths=FALSE)
{
	tmpVertices <- NULL
	#print("Getting nodes in the neighborhood...")
	theHood <- igraph::neighborhood(graph=iGraph, nodes=selected_vertices, order=vertex_level)
	for(i in 1:length(theHood))
	{
		tmpVertices <- c(tmpVertices, theHood[[i]])
	}
	
	#print("Getting Shortest Paths!")
	for(i in 1:length(selected_vertices))
	{
		fromVertex <- selected_vertices[i]
		#print("Pre shortes_path function!")
		#pathObject <- shortest_paths(graph=iGraph, from=fromVertex, to=selected_vertices)
		pathObject <- shortest_paths(graph=iGraph, from=fromVertex, to=selected_vertices[-i])
		#print("Post shortes_path function!")
		pathVertices <- pathObject$vpath
		for(j in 1:length(pathVertices))
		{
			tmpVertices <- c(tmpVertices, pathVertices[[j]])
		}
	}
	#print("Got All vertices!")

	vertices4SubGraph <- unique(tmpVertices)
	subiGraph <- induced.subgraph(graph=iGraph, vids=vertices4SubGraph)
	print(list.vertex.attributes(subiGraph))
	
	#vertices2Highlight <- which(get.vertex.attribute(subiGraph,"name") %in% selected_vertices)
	#igraph::V(subiGraph)$color[vertices2Highlight] <- igraph::V(subiGraph)$highlightcolor[vertices2Highlight]
	subiGraph 
}

#' Get modules from the main igraph.
#'
#' Extract modules from the main graph by using a specified community detection algorithm from igraph.
#'
#' @importFrom igraph cluster_walktrap cluster_spinglass cluster_louvain cluster_fast_greedy
#'
#' @param iGraph igraph object representing the main graph from which the subgraph should be extracted.
#' @param method community detection method from igraph.
#' @return communities object cotaining communities identified from igraph object by using a specific method
#' @examples
#' \dontrun{
#' get_modules(iGraph=main.graph, method="walktrap")
#' }
#' @keywords internal
#' @export
get_modules <- function(iGraph, method="walktrap")
{
	switch(method,
	"walktrap" = {
		optimalStep <- c(2:10)[which.max(sapply(c(2:10), function(x){igraph::modularity(igraph::cluster_walktrap(iGraph, step=x))}))]
		igraph::cluster_walktrap(iGraph, step=optimalStep)
	},
	"spinglass" = igraph::cluster_spinglass(iGraph),
	"louvain" = igraph::cluster_louvain(iGraph),
	"greedy" = igraph::cluster_fast_greedy(iGraph)
       ) 
}

#' Annotate modules identified from the main igraph.
#'
#' Annotate modules and return list-of-list containing igraph objects, annotation enrichment and ranks.
#'
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @importFrom AnnotationDbi select
#' @importFrom igraph communities induced.subgraph get.edgelist
#'
#' @param iGraph igraph object representing the main graph from which the subgraph should be extracted.
#' @param modules Modules detected from the main igraph with specific community detection method.
#' @param rl Ordered ranked list of genes obtained by taking median rank from Borda over all scores.
#' @param rl.c Ordered ranked list of genes obtained by taking median rank from Borda over centrality scores.
#' @param rl.pv Ordered ranked list of genes by P.Value obtained from differential expression analysis.
#' @param rl.lfc Ordered ranked list of genes by Log FC obtained from differential expression analysis.
#' @param rl.edge Ordered ranked list of edges obtained by taking median rank from Borda over ranked edges from matrices inferred by different algorithms.
#' @param annDB Organism specific annotation library default:'org.Hs.eg.db'.
#' @param p_value P.Value cutoff for the enriched GO terms default:0.05.
#' @param min_mod_size Minimum module size to select modules for annotation and to be added in results; DEFAULT:10.
#' @param prefix string to prefix module names; DEFAULT:NULL.
#' @return list-of-list containing igraph objects, annotation enrichment and median values for vertice and edges by different rank lists in each module
#' @examples
#' \dontrun{
#' annotate_modules(iGraph=main.graph,
#' modules,
#' rl,
#' rl.c,
#' rl.pv=NULL,
#' rl.lfc=NULL,
#' rl.edge=NULL,
#' annDB="org.Hs.eg.db",
#' p_value=0.05,
#' min_mod_size=10,
#' prefix=NULL
#' )
#' }
#' @keywords internal
#' @export
annotate_modules <- function(iGraph, modules, rl, rl.c, rl.pv=NULL, rl.lfc=NULL, rl.edge=NULL, annDB="org.Hs.eg.db", p_value=0.05, min_mod_size=10, prefix=NULL){
	res_ll <- list()
	res_ll[["names"]] <- igraph::communities(modules)

	#Filter modules by size	
	len_vec <- as.vector(unlist(lapply(res_ll[["names"]], length)))
	idx <- which(len_vec<min_mod_size)
	if(length(idx)>0){
		res_ll[["names"]] <- res_ll[["names"]][-idx]
	}

        if(!is.null(prefix)){
                names(res_ll[["names"]]) <- paste0(prefix, "_", c(1:length(res_ll[["names"]])))
        }
	res_ll[["ig"]] <- lapply(res_ll[["names"]], function(x){igraph::induced.subgraph(iGraph, x)})
	res_ll[["ae"]] <- lapply(
		res_ll[["names"]], function(x){
			res <- annotation_enrichment(genelist=x, annDB=annDB)
			#res <- res[res$ease<=0.05,]
			res <- res[res$ease<=p_value,]
			if(nrow(res)==0)
			return(NULL)
		
			selectAnnDF <- AnnotationDbi::select(GO.db, keys=as.vector(res$ID), columns=c("TERM", "ONTOLOGY"), keytype="GOID")
			#res <- cbind(selectAnnDF, res[,c(3,2)])
			#colnames(res)[c(4,5)] <- c("Score", "Genes")
			res <- cbind(selectAnnDF, res[,c("pValue", "ease", "User_Genes", "Symbols")])
			colnames(res)[c(4:ncol(res))] <- c("P_Value", "EASE_Score", "Gene_Count", "Genes")
			res
		}
	)
	res_ll[["ranks"]] <- lapply(
		res_ll[["names"]],
		function(x){
			tmp_list <- list()
			tmp_list[["rank_combined"]] <- round(1-(median(which(rl %in% x))/length(rl)),2)
			tmp_list[["rank_centrality"]] <- round(1-(median(which(rl.c %in% x))/length(rl.c)),2)
                        if(!is.null(rl.pv)){
			        tmp_list[["rank_pv"]] <- round(1-(median(which(rl.pv %in% x))/length(rl.pv)),2)
                        }
                        if(!is.null(rl.lfc)){
			        tmp_list[["rank_lfc"]] <- round(1-(median(which(rl.lfc %in% x))/length(rl.lfc)),2)
                        }
                        if(!is.null(rl.edge)){
                                el <- apply(igraph::get.edgelist(igraph::induced.subgraph(iGraph, x)), 1, paste, collapse=";")
                                tmp_list[["rank_edge"]] <- round(1-(median(which(rl.edge %in% el))/length(rl.edge)),2)
                        }
			tmp_list
		}
	)
	return(res_ll)
}

#' Get table containing score for genes in the iGraph by th euser provided ranked lists.
#'
#' Get table containing score for genes in the iGraph by th euser provided ranked lists.
#'
#' @importFrom igraph communities induced.subgraph get.edgelist
#'
#' @param iGraph igraph object representing the main graph from which the subgraph should be extracted.
#' @param gene_ll Named list of character vectors containing gene symbols, representing the genes present in modules.
#' @param rl Ordered ranked list of genes obtained by taking median rank from Borda over all scores.
#' @param rl.c Ordered ranked list of genes obtained by taking median rank from Borda over centrality scores.
#' @param rl.pv Ordered ranked list of genes by P.Value obtained from differential expression analysis.
#' @param rl.lfc Ordered ranked list of genes by Log FC obtained from differential expression analysis.
#' @param rl.edge Ordered ranked list of edges obtained by taking median rank from Borda over ranked edges from matrices inferred by different algorithms.
#' @return Data frame containing median(ranks), size and GO within similarity.
#' @examples
#' \dontrun{
#' get_rank_table(iGraph=main.graph,
#' gene_ll,
#' rl,
#' rl.c,
#' rl.pv=NULL,
#' rl.lfc=NULL,
#' rl.edge=NULL)
#' }
#' @keywords internal
#' @export
get_rank_table <- function(iGraph, gene_ll, rl, rl.c, rl.pv, rl.lfc, rl.edge){
	#res_ll <- list()
	res_ll <- lapply(
		gene_ll,
		function(x){
			res <- c()
			res["rank_combined"] <- round(1-(median(which(rl %in% x))/length(rl)),2)
			res["rank_centrality"] <- round(1-(median(which(rl.c %in% x))/length(rl.c)),2)
                        if(!is.null(rl.pv)){
			        res["rank_pv"] <- round(1-(median(which(rl.pv %in% x))/length(rl.pv)),2)
                        }
                        if(!is.null(rl.lfc)){
			        res["rank_lfc"] <- round(1-(median(which(rl.lfc %in% x))/length(rl.lfc)),2)
                        }
                        if(!is.null(rl.edge)){
                                el <- apply(igraph::get.edgelist(igraph::induced.subgraph(iGraph, x)), 1, paste, collapse=";")
                                res["rank_edge"] <- round(1-(median(which(rl.edge %in% el))/length(rl.edge)),2)
                        }
                        res["size"] <- round(length(x)/igraph::vcount(iGraph),2)
			res
		}
	)
        resDF <- as.data.frame(t(as.data.frame(res_ll)))

	#resDF[,"GO"] <- rep(NA, nrow(resDF))
	#tmp_vec <- unlist(lapply(ae_ll, function(x) get_jaccard_sim(x, x, IC_ll=IC_ll)))
	#idx <- which(rownames(resDF) %in% names(tmp_vec))
	#if(length(idx)>0){
	#	resDF[idx,"GO"] <- tmp_vec
	#}

	return(resDF)
}

#' Convert annotated modules list-of-list to data frame containing rank and size.
#'
#' Converts the list-of-list object of annotated modules and returns a data frame with annotation information, ready to be used for plotting.
#'
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @importFrom AnnotationDbi select
#' @importFrom igraph vcount
#'
#' @param iGraph igraph object representing the main graph from which the subgraph should be extracted.
#' @param annotated_modules_ll list-of-list containing igraph objects, annotation enrichment and ranks for each module.
#' @param IC_ll list-of-list containing GOSemSimDATA objects for BP, MF and CC.
#' @param prefix string to prefix module names; DEFAULT:NULL.
#' @return Data frame containing median(ranks), size and GO within similarity.
#' @examples
#' \dontrun{
#' as.rank.table(iGraph=main.graph, annotated_modules_ll, IC_ll=NULL, prefix=NULL)
#' }
#' @keywords internal
#' @export
as.rank.table <- function(iGraph, annotated_modules_ll, IC_ll=NULL, prefix=NULL){
	resDF <- as.data.frame(t(as.data.frame(lapply(annotated_modules_ll[["ranks"]], unlist))))
        if(!is.null(prefix)){
                rownames(resDF) <- paste0(prefix, "_", c(1:nrow(resDF)))
        }
	resDF[,"size"] <- unlist(lapply(annotated_modules_ll[["names"]], function(x){round(length(x)/igraph::vcount(iGraph),2)}))

	#resDF[,"GO"] <- rep(NA, nrow(resDF))
	#tmp_vec <- unlist(lapply(annotated_modules_ll[["ae"]], function(x) get_jaccard_sim(x, x, IC_ll=IC_ll)))
        #if(!is.null(prefix)){
        #        names(tmp_vec) <- paste0(prefix, "_", names(tmp_vec))
        #}
	#idx <- which(rownames(resDF) %in% names(tmp_vec))
	#if(length(idx)>0){
	#	resDF[idx,"GO"] <- tmp_vec
	#}

	return(resDF)
}

#' Get Jaccard Index.
#'
#' Get Jaccard index for two sets (vectors) of identifiers.
#'
#' @param set1 Set of identifiers provided as a vector.
#' @param set2 Set of identifiers provided as a vector.
#' @return Numerical value between 0 and 1 representing the Jaccard similarity coefficient.
#' @examples
#' \dontrun{
#' get_jaccard(set1, set2)
#' }
#' @keywords internal
#' @export
get_jaccard <- function(set1, set2){
	val <- length(intersect(set1, set2))/length(union(set1, set2))
	return(val)
}

#' Get Jaccard Index with Sematic Support.
#'
#' Get Jaccard index for two sets (data frame) of annotation enrichment results with the support of semantic similarity between GO terms to perform intersection.
#'
#' @importFrom GOSemSim godata mgoSim
#' @importFrom AnnotationDbi select
#' @importFrom stats setNames
#'
#' @param set1 Set of identifiers provided as a vector.
#' @param set2 Set of identifiers provided as a vector.
#' @param simThr Threshold of semantic similarity score to consider items as similar; DEFAULT:0.6.
#' @param IC_ll list-of-list containing GOSemSimDATA objects for BP, MF and CC.
#' @param annDB Organism specific annotation library default:'org.Hs.eg.db'.
#' @param simMeasure Sematic similarity measure option from GOSemSim package default:"Rel".
#' @return Numerical value between 0 and 1 representing the Jaccard similarity coefficient.
#' @examples
#' \dontrun{
#' get_jaccard_sim(set1, set2, simThr=0.6, IC_ll=NULL, annDB="org.Hs.eg.db", simMeasure="Rel")
#' }
#' @keywords internal
#' @export
get_jaccard_sim <- function(set1, set2, simThr=0.6, IC_ll=NULL, annDB="org.Hs.eg.db", simMeasure="Rel"){
        if(is.null(set1) || is.null(set2))
        return(NULL)
        
	setU <- unique(c(set1[,"GOID"], set2[,"GOID"]))
	goVec <- setNames(rep(0, length(setU)), setU)
	
	GODB_vals = AnnotationDbi::select(GO.db, keys(GO.db, "GOID"), c("TERM", "ONTOLOGY"))
        GODB_vals_ll <- list()
        GODB_vals_ll[["BP"]] <- GODB_vals[GODB_vals$ONTOLOGY=="BP",c("GOID", "TERM")]
        GODB_vals_ll[["CC"]] <- GODB_vals[GODB_vals$ONTOLOGY=="CC",c("GOID", "TERM")]
        GODB_vals_ll[["MF"]] <- GODB_vals[GODB_vals$ONTOLOGY=="MF",c("GOID", "TERM")]

        simScore_list <- list()
        for(typ in c("BP", "CC", "MF")){
                #print(paste0("###### Computing For Ontology - ", typ, ' ######'))
                ont1 <- set1[set1$ONTOLOGY == typ,]
                ont2 <- set2[set2$ONTOLOGY == typ,]
                
                if(nrow(ont1)==0 || nrow(ont2)==0)
                next
                
		if(is.null(IC_ll) || length(which(names(IC_ll) %in% typ))==0){
			d <- GOSemSim::godata(annDB, ont=typ, computeIC=TRUE)
		}else{
			d <- IC_ll[[typ]]
		}
                
                ont_ids1 <- as.character(ont1[which(ont1$GOID %in% GODB_vals_ll[[typ]]$GOID), "GOID"])
                ont_ids2 <- as.character(ont2[which(ont2$GOID %in% GODB_vals_ll[[typ]]$GOID), "GOID"])

		scoreMat <- GOSemSim::mgoSim(ont_ids1, ont_ids2, semData=d, measure=simMeasure, combine=NULL)
		scoreDF <- as.data.frame(as.table(scoreMat), stringsAsFactors=FALSE)
		idx <- which(scoreDF[,3]>=simThr)
		if(length(idx)>0){
			simID <- unique(c(scoreDF[idx,1], scoreDF[idx,2]))
			goVec[simID] <- 1
		}
	}
	
	valI <- length(which(goVec==1))
	valU <- length(setU)
	val <- valI/valU
	return(val)
}

#' Get similarity between all pairs of elements between two lists of enriched annotations by using modified Jaccard with semsim.
#'
#' Get similarity between all pairs of elements between two lists of enriched annotations by using modified Jaccard with semsim.
#'
#' @importFrom GOSemSim godata mgoSim
#'
#' @param mod_ll1 list-of-list containing igraph objects, annotation enrichment and ranks for each module.
#' @param mod_ll2 list-of-list containing igraph objects, annotation enrichment and ranks for each module.
#' @param prefix1 prefix string for module names in mod_ll1; DEFAULT:NULL.
#' @param prefix2 prefix string for module names in mod_ll2; DEFAULT:NULL.
#' @param IC_ll list-of-list containing GOSemSimDATA objects for BP, MF and CC.
#' @return Numerical matrix representing the Jaccard similarity coefficient between pairs of mods from mod_ll1 vs mod_ll2.
#' @examples
#' \dontrun{
#' get_between_GO_sim(mod_ll1, mod_ll2, prefix1=NULL, prefix2=NULL, IC_ll=NULL)
#' }
#' @keywords internal
#' @export
get_between_GO_sim <- function(mod_ll1, mod_ll2, prefix1=NULL, prefix2=NULL, IC_ll=NULL){
	ae_list1 <- mod_ll1[["ae"]]
	ae_list2 <- mod_ll2[["ae"]]
	l1 <- length(ae_list1)
	l2 <- length(ae_list2)
	
	#if(is.null(prefix1))
	#prefix1 <- "A"
	#
	#if(is.null(prefix2))
	#prefix2 <- "B"

	names1 <- paste0(prefix1, names(ae_list1))
	#print(names1)
	names(ae_list1) <- names1
	names2 <- paste0(prefix2, names(ae_list2))
	#print(names2)
	names(ae_list2) <- names2
	
	res <- matrix(data=rep(0, l1*l2), nrow=l1, ncol=l2, dimnames=list(names1, names2))
	for(m1 in names(ae_list1)){
		if(is.null(ae_list1[[m1]]))
		next
		
		for(m2 in names(ae_list2)){
			if(is.null(ae_list2[[m2]]))
			next
		
			res[m1,m2] <- get_jaccard_sim(ae_list1[[m1]], ae_list2[[m2]], IC_ll=IC_ll)
			#print(res[m1,m2])
		}
	}
	return(res)
}

#' Get color gradients.
#'
#' Get color gradients by providing the reference colors for the lower and higher spectrum of the gradient and the size of the gradient.
#'
#' @importFrom grDevices colorRampPalette
#'
#' @param low Reference color for the lower end of the gradient spectrum; must be a valid argument to 'col2rgb()'.
#' @param high Reference color for the higher end of the gradient spectrum; must be a valid argument to 'col2rgb()'.
#' @param ncolors Number of colors in the gradient composition.
#' @return Gradient of colors as a vector.
#' @examples
#' \dontrun{
#' get_col_gradient(low="red", high="yellow", ncolors=123)
#' }
#' @keywords internal
#' @export
get_col_gradient <- function(low="red", high="yellow", ncolors=123) {
	rbPal <- grDevices::colorRampPalette(c(low,high))
	colGradient <- rbPal(ncolors)
	return(colGradient)
}

#' Plot heatmap of the similarity matrix.
#'
#' Plot heatmap of the similarity matrix.
#'
#' @importFrom stats as.dist hclust as.dendrogram
#' @importFrom gplots heatmap.2
#'
#' @param hSim Similarity matrix for plotting heatmap.
#' @param rSim Similarity matrix for clustering the items on the rows of hSim.
#' @param cSim Similarity matrix for clustering the items on the columns of hSim.
#' @return Heatmap plot with row and column clusters.
#' @examples
#' \dontrun{
#' plot_sim_heatmap(hSim, rSim=NULL, cSim=NULL)
#' }
#' @keywords internal
#' @export
plot_sim_heatmap <- function(hSim, rSim=NULL, cSim=NULL){
	rDend <- NULL
	if(!is.null(rSim)){
		rDat <- rSim
		rDat <- rDat[rowSums(rDat)>0,] 
		rDat <- rDat[,colSums(rDat)>0]
		rd <- stats::as.dist(1-rDat)
		rc <- stats::hclust(rd)
		rDend <- stats::as.dendrogram(rc)
	
	}

	cDend <- NULL
	if(!is.null(cSim)){
		cDat <- cSim
		cDat <- cDat[rowSums(cDat)>0,] 
		cDat <- cDat[,colSums(cDat)>0]
		cd <- stats::as.dist(1-t(cDat))
		cc <- stats::hclust(cd)
		cDend <- stats::as.dendrogram(cc)
	}

	hot_data <- hSim
	hot_data <- hot_data[rowSums(hot_data)>0,] 
	hot_data <- hot_data[,colSums(hot_data)>0]
	hot_data <- 1-hot_data

	#lowcol <- "yellow"
	#highcol <- "red"
	lowcol <- "red"
	highcol <- "yellow"
	plotcols <- get_col_gradient(lowcol, highcol)
	gplots::heatmap.2(hot_data,
		Rowv=rDend,
		Colv=cDend,
		dendrogram="both",
		scale="none",
		col=plotcols,
		cexCol=1,
		margins=c(8,8),
		zlim=c(0,1),
		trace="none"
	)
}

#' Update selected vertices with highlight color.
#'
#' Alter the color information in the selected nodes, when plotted these vertices will identified with the alternate color from the rest of the vertices.
#'
#' @importFrom igraph neighborhood V induced.subgraph
#'
#' @param iGraph igraph object representing the main graph from which the subgraph should be extracted.
#' @param selected_vertices Vector of genes to be used for identifying interconnected hub.
#' @return igraph object with highlight color for the selected vertices
#' @examples
#' \dontrun{
#' highlight_vertices(iGraph=main.graph, selected_vertices)
#' }
#' @keywords internal
#' @export
highlight_vertices <- function(iGraph, selected_vertices)
{
        vertices2Highlight <- which(get.vertex.attribute(iGraph,"name") %in% selected_vertices)
	igraph::V(iGraph)$color[vertices2Highlight] <- igraph::V(iGraph)$highlightcolor[vertices2Highlight]
	iGraph
}

#' Get enriched GO.
#'
#' Submit a list of genes and specify the ogranism specific annotation library to get enriched GO.
#'
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @importFrom AnnotationDbi select keys
#' @importFrom GSEABase getOBOCollection ids
#' @importFrom plyr ddply . summarise
#' @importFrom stats fisher.test
#'
#' @param genelist List of gene symbols for which you want enriched GO.
#' @param annType Type of GO annotation 'GO' or 'GO_SLIM' default:'GO'.
#' @param annDB Organism specific annotation library default:'org.Hs.eg.db'.
#' @return Data frame containing enriched GO with their corresponding pValue, EASE score and list of genes representing GO
#' @examples
#' \dontrun{
#' annotation_enrichment(genelist=list_of_genes.df, annType="GO", annDB="org.Hs.eg.db")
#' annotation_enrichment(genelist=list_of_genes.df)
#' }
#' @keywords internal
#' @export
annotation_enrichment <- function(genelist, annType="GO", annDB="org.Hs.eg.db"){
	keyType <- "GO"
	annDB <- get(annDB)
	genelist.len <- length(genelist)
	print("Creating Res DF...")
	selectDF <- AnnotationDbi::select(annDB, keys=genelist, columns=keyType, keytype="SYMBOL")
	selectDF <- selectDF[selectDF$EVIDENCE != "ND",]
	if(annType == "GO_SLIM"){
		print("Use GO SLIM Instead of ALL!")
		go_slim_generic_obo <- GSEABase::getOBOCollection("http://www.geneontology.org/ontology/subsets/goslim_generic.obo")
		selectDF <- selectDF[selectDF$GO %in% GSEABase::ids(go_slim_generic_obo),]
	}
	selectDF <- selectDF[,c(1,2)]
	if(length(which(is.na(selectDF[,2]))) > 0){
		print("Remove NA. Unmapped to annotation...")
		selectDF <- selectDF[-which(is.na(selectDF[,2])),]
	}
	colnames(selectDF) <- c("Symbol", "ID")
	selectDF <- unique(selectDF)
	#resDF <- ddply(selectDF,.(ID),nrow)
	#colnames(resDF) <- c("ID", "User_Genes")
	Symbol <- c()
	ID <- c()
	resDF <- plyr::ddply(selectDF, plyr::.(ID), plyr::summarise, paste(Symbol, collapse=","), length(Symbol))
	colnames(resDF) <- c("ID", "Symbols", "User_Genes")
	resDF$User_Genes_Not_In <- genelist.len-resDF$User_Genes
	print("Getting All Genes for mapped IDs...")
	selectAnnDF <- AnnotationDbi::select(annDB, keys=resDF$ID, columns="ENTREZID", keytype=keyType)
	colnames(selectAnnDF)[1] <- "ID"
	print("Getting Gene count for each ID...")
	resDF$ALL_Genes <- plyr::ddply(selectAnnDF,.(ID),nrow)[,2]
	print("Getting Total Gene Count...")
	genome.len <- length(unique(keys(annDB, keytype="ENTREZID")))
	resDF$ALL_Genes_Not_In <- genome.len-resDF$ALL_Genes
	print("Calculating PValue...")
	resDF$pValue <- apply(resDF, 1, function(x){
		#xx <- as.numeric(x[2:5])
		xx <- as.numeric(x[3:6])
		mat <- matrix(c(xx[1],xx[2],xx[3],xx[4]), nrow=2, dimnames=list(c("In", "Not_In"), c("User", "All")))
		testRes <- fisher.test(mat)
		testRes$p.value
	})
	print("Calculating EASE...")
	resDF$ease <- apply(resDF, 1, function(x){
		#xx <- as.numeric(x[2:5])
		xx <- as.numeric(x[3:6])
		genesIn <- xx[1]-1
		mat <- matrix(c(genesIn,xx[2],xx[3],xx[4]), nrow=2, dimnames=list(c("In", "Not_In"), c("User", "All")))
		testRes <- stats::fisher.test(mat)
		testRes$p.value
	})
	print("Completed. Returning...")
	resDF
}

#' Get progressive enrichment of GO annotation across 3 sets of different size.
#'
#' Get enriched annotations table from annotation_enrichment for three different size subgraphs and provide it to this function in the order of increasing size.
#'
#' @importFrom AnnotationDbi select
#'
#' @param genelists List of list containing the genes used to perform GO enrichment for the three sub-networks.
#' @param l1 Enriched GO annotation table for the small set.
#' @param l2 Enriched GO annotation table for the medium set.
#' @param l3 Enriched GO annotation table for the large set.
#' @return Data frame containing progressively enriched GO along with an enrichment score.
#' @examples
#' \dontrun{
#' progressive_enrichment(genelists, l1=enrichment_small, l2=enrichment_medium, l3=enrichment_large)
#' }
#' @keywords internal
#' @export
progressive_enrichment <- function(genelists, l1, l2, l3)
{
	print("Filtering the IDs, based on gene count...")
	l1_fltr <- l1[which(l1$User_Genes >=3), c(1,3,7)]
	l2_fltr <- l2[which(l2$User_Genes >=3), c(1,3,7)]
	l3_fltr <- l3[which(l3$User_Genes >=3), c(1,3,7)]
	
	print("~~ Creating Score=round(li_fltr$User_Genes/(length(genelists$li)/length(li_fltr$ID))) ~~")
	l1_fltr$Score <- round(l1_fltr$User_Genes/(length(genelists[[1]])/length(l1_fltr$ID)))
	l2_fltr$Score <- round(l2_fltr$User_Genes/(length(genelists[[2]])/length(l2_fltr$ID)))
	l3_fltr$Score <- round(l3_fltr$User_Genes/(length(genelists[[3]])/length(l3_fltr$ID)))

	l1_fltr <- l1_fltr[,c(1,2,4)]
	l2_fltr <- l2_fltr[,c(1,2,4)]
	l3_fltr <- l3_fltr[,c(1,2,4)]

	print("Set column names for mapping...")
	colnames(l1_fltr) <- c("ID", "L1_Count", "L1_Score")
	colnames(l2_fltr) <- c("ID", "L2_Count", "L2_Score")
	colnames(l3_fltr) <- c("ID", "L3_Count", "L3_Score")

	all_fltr_union <- Reduce(function(...) merge(..., all=TRUE), list(l1_fltr, l2_fltr, l3_fltr))
	all_fltr_union[is.na(all_fltr_union)] <- 0
	all_fltr_union <- all_fltr_union[which(all_fltr_union$L1_Score <= all_fltr_union$L2_Score & all_fltr_union$L2_Score <= all_fltr_union$L3_Score),]
	
	all_fltr_union_score <- cbind(all_fltr_union[1], round(rowSums(all_fltr_union[c("L1_Score","L2_Score","L3_Score")])))
	colnames(all_fltr_union_score) <- c("ID", "Score")
	all_fltr_union_score <- cbind(all_fltr_union_score, unname(sapply(all_fltr_union_score$ID, function(x) l3[which(l3$ID == x),2])))
	colnames(all_fltr_union_score) <- c("ID", "Score", "Genes")
	selectAnnDF <- AnnotationDbi::select(GO.db, keys=as.vector(all_fltr_union_score$ID), columns=c("TERM", "ONTOLOGY"), keytype="GOID")
	all_fltr_union_score <-  cbind(selectAnnDF, all_fltr_union_score[,c(2,3)])
	all_fltr_union_score
}

#' Summarize the GO information from the progressively enriched GO table.
#'
#' go_summarization takes as input the result from progressive_enrichment() function and segregrates the GO terms by type "BP", "CC" and "MF". 
#' For each type of GOs the semantic similarity is computed to summarize and cluster GO terms. The clustered GOs are reported as data frame, 
#' GO term with the highest score in each cluster is chosen as the representative of that cluster.
#'
#' @importFrom AnnotationDbi select keys
#' @importFrom GOSemSim mgoSim
#' @importFrom stats hclust as.dist cutree
#'
#' @param enriched_GO_DF Data frame containing progressively enriched GO along with an enrichment score.
#' @param score_col column from the enrichment data frame to use as score.
#' @param annDB Organism specific annotation library default:'org.Hs.eg.db'.
#' @param simMeasure GOSemSim measure to use for computing semantic similarity between GOs.
#' @param treeHeight Value used to cut the clustered tree and define clusters.
#' @param IC_ll list-of-list containing GOSemSimDATA objects for BP, MF and CC.
#' @param log_transform Log transform the score as -log10(score), for P.Values; DEFAULT:TRUE.
#' @return List of data frames containing clustered GO information, with major terms and individual terms. One data frame for each GO type BP, CC and MF.
#' @examples
#' \dontrun{
#' go_summarization(enriched_GO_DF,
#' score_col="EASE_Score",
#' annDB="org.Hs.eg.db",
#' simMeasure="Rel",
#' treeHeight=0.9,
#' IC_ll=NULL,
#' log_transform=TRUE
#' )
#' }
#' @keywords internal
#' @export
go_summarization <- function(enriched_GO_DF, score_col="EASE_Score", annDB="org.Hs.eg.db", simMeasure="Rel", treeHeight=0.9, IC_ll=NULL, log_transform=TRUE){
	GODB_vals = AnnotationDbi::select(GO.db, keys(GO.db, "GOID"), c("TERM", "ONTOLOGY"))
	GODB_vals_ll <- list()
	GODB_vals_ll[["BP"]] <- GODB_vals[GODB_vals$ONTOLOGY=="BP",c("GOID", "TERM")]
	GODB_vals_ll[["CC"]] <- GODB_vals[GODB_vals$ONTOLOGY=="CC",c("GOID", "TERM")]
	GODB_vals_ll[["MF"]] <- GODB_vals[GODB_vals$ONTOLOGY=="MF",c("GOID", "TERM")]

	simClust_DF <- list()
	for(typ in c("BP", "CC", "MF")){
		print(paste0("###### Clustering For Ontology - ", typ, ' ######'))
		ont <- enriched_GO_DF[enriched_GO_DF$ONTOLOGY == typ,]
		if(log_transform==TRUE){
			print("*** LOG TRANSFORMING ***")
			ont[,score_col] <- -log10(ont[,score_col])
		}
		ont_ids <- as.character(ont[which(ont$GOID %in% GODB_vals_ll[[typ]]$GOID), "GOID"])

		#simsem_matrix <- GOSemSim::mgoSim(ont_ids, ont_ids,  measure="Rel", ont=typ, organism="human", combine=NULL)

                #d <- godata(annDB, ont=typ, computeIC=TRUE)
		if(is.null(IC_ll) || length(which(names(IC_ll) %in% typ))==0){
			d <- GOSemSim::godata(annDB, ont=typ, computeIC=TRUE)
		}else{
			d <- IC_ll[[typ]]
		}

		simsem_matrix <- GOSemSim::mgoSim(ont_ids, ont_ids, semData=d, measure=simMeasure, combine=NULL)
                if(is.null(simsem_matrix)){
                        next
                }
		dim_vec <- dim(simsem_matrix)
                #print(dim_vec)
                if(is.null(dim_vec)){
                        next
                }
		dim_len <- dim_vec[1]
		#print(dim_len)
                valRowIdx <- which(rowSums(is.na(simsem_matrix))!=dim_len)
                valColIdx <- which(colSums(is.na(simsem_matrix))!=dim_len)
                if(length(valRowIdx)<2 || length(valColIdx)<2){
                        next
                }

		#simsem_matrix_clean <- simsem_matrix[rowSums(is.na(simsem_matrix))!=dim_len, colSums(is.na(simsem_matrix))!=dim_len]
		simsem_matrix_clean <- simsem_matrix[valRowIdx, valColIdx]
		simsem_matrix <- round(simsem_matrix_clean, 2)
                #print(str(simsem_matrix))
                #print(valRowIdx)
                #print(valColIdx)
		
		simClust <- stats::hclust(stats::as.dist(1-simsem_matrix))

		cut_height <- treeHeight
		simClust_cut <- stats::cutree(simClust, h=cut_height)
		simClust_DF[[typ]] <- data.frame(ID=character(), TERM=character(), Score=integer(), Representative=character(), stringsAsFactors=FALSE)
		
		for(i in unique(simClust_cut)){
			memberIDs <- names(simClust_cut[simClust_cut == i])
			#selRow <- ont[which(ont$GOID %in% memberIDs), c("TERM", "Score")]
                        #clustID <- selRow[which.max(selRow$Score), "TERM"]
			#localDF <- data.frame(ID=memberIDs, TERM=selRow$TERM, Score=selRow$Score, Representative=clustID, stringsAsFactors=F)
			selRow <- ont[which(ont$GOID %in% memberIDs), c("TERM", score_col)]
                        clustID <- selRow[which.max(selRow[,score_col]), "TERM"]
			localDF <- data.frame(ID=memberIDs, TERM=selRow$TERM, Score=selRow[,score_col], Representative=clustID, stringsAsFactors=FALSE)
			simClust_DF[[typ]] <- rbind(simClust_DF[[typ]], localDF)
		}
		print(paste0("###### Completed CLustering For Ontology - ", typ, ' ######'))
	}
	return(simClust_DF)
}

#' Plot the summarized progressively enriched GO table.
#'
#' Plot the summarized GO terms as tiles with area equivalent to the score, major GO term from each cluster is displayed as the representative.
#'
#' @importFrom treemap treemap
#'
#' @param input.DF Data frame containing the summarized progressively enriched GO of a specific type.
#' @param ont Type of GO term provided in the data frame input.DF, default: "BP".
#' @return treemap/tileplot representation of the summarized GO terms.
#' @examples
#' \dontrun{
#' plot_treemap(input.DF, ont="BP")
#' }
#' @keywords internal
#' @export
plot_treemap <- function(input.DF, ont="BP"){
	if(ont=="BP"){
		tmap_title <- "Gene Ontology: Biological Process"
	}else if(ont=="MF"){
		tmap_title <- "Gene Ontology: Molecular Function"
	}else if(ont=="CC"){
		tmap_title <- "Gene Ontology: Cellular Component"
	}

	treemap::treemap(
		input.DF,
		index = c("Representative","TERM"),
		vSize = "Score",
		type = "categorical",
		vColor = "Representative",
		title = tmap_title,
		inflate.labels = FALSE,
		lowerbound.cex.labels = 0,
		bg.labels = "#CCCCCCAA",
                #fontsize.labels = c(0,1), #Remove Category Name
                #inflate.labels = TRUE,
		position.legend = "none"
	)
}

#' Get visNetwork plot object representing the input igraph object.
#'
#' Convert the igraph into an interactive visNetwork object for plotting.
#'
#' @importFrom visNetwork visIgraph visIgraphLayout visOptions visEvents
#'
#' @param iG igraph object representing the inferred consensus matrix and annotated with the user specified aesthetics information.
#' @param plot_layout Graph layouts applicable on an igraph object, default: "nicely".
#' @param vBorderColor Valid R color for vertex borders, default: "black".
#' @param vShape Type of shape supported by visNetwork, default: "circle".
#' @param vFontColor Valid R color for vertex label, default: "#343434".
#' @param vSize Size of the vertext label font, default: 50.
#' @param eWidth Thickness of the edges connecting the vertices, default: 10.
#' @return visNetwork plot object created from the igraph information.
#' @examples
#' \dontrun{
#' get_visNetwork(iG,
#' plot_layout="nicely",
#' vBorderColor="black",
#' vShape="circle",
#' vFontColor="#343434",
#' vSize=50,
#' eWidth=10,
#' degDepth=2
#' )
#' }
#' @keywords internal
#' @export
#get_visNetwork <- function(iG, plot_layout, vBorderColor="black", vShape="circle", vFontColor="#343434", vSize=50, eWidth=10, vColTop="red", vColBottom="blue", vColType="score"){
get_visNetwork <- function(iG, plot_layout="nicely", vBorderColor="black", vShape="circle", vFontColor="#343434", vSize=50, eWidth=10, degDepth=2){
    visLayout <- paste0("layout_", plot_layout)
    vNet <- visNetwork::visIgraph(iG) %>%
        visNetwork::visIgraphLayout(layout = visLayout) %>%
        #visNetwork::visOptions(highlightNearest = list(enabled = FALSE, hover = FALSE)) %>% 
        visNetwork::visOptions(highlightNearest = list(enabled = TRUE, hover = TRUE, degree = degDepth)) %>% 
        visNetwork::visEvents(selectNode = "function(properties) {
            url='http://www.ncbi.nlm.nih.gov/gene?term=('+ properties.nodes +'%5BGene%20Name%5D)%20AND%20%22homo%20sapiens%22%5BOrganism%5D'
            window.open(url);}"
        )

    #if(vColType=="score"){
    vNet$x$nodes$color.background <- vNet$x$nodes[,"color"]
    #}else if(vColType=="rank"){
    #        rbPal <- colorRampPalette(c(vColTop,vColBottom))
    #        vNet$x$nodes$color.background <- rbPal(vcount(iG))
    #}
    vNet$x$nodes$color.border <- vBorderColor
    vNet$x$nodes[,c("color", "highlightcolor")] <- list(NULL)
    vNet$x$nodes$title <- paste0("<p><b>", vNet$x$nodes$id,"</b></p>")
    vNet$x$nodes$shape <- vShape
    vNet$x$nodes$value <- vSize
    vNet$x$nodes$font.size <- vSize
    if(vShape=="text"){
        vNet$x$nodes$font.color <- vNet$x$nodes$color.background
    }else{
        vNet$x$nodes$font.color <- vFontColor
    }
    #vNet$x$nodes$scaling.enabled <- FALSE
    vNet$x$options$nodes$shape <- NULL
    vNet$x$edges$width <- eWidth

    return(vNet)
}

#' Get master gene table for the inferred network with ranks, scores and module information.
#'
#' Get master gene table for the inferred network with ranks, scores and module information.
#'
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @importFrom AnnotationDbi select get
#'
#' @param ig igraph object representing the main graph from which the subgraph should be extracted.
#' @param rankedGenes vector of genes ordered by rank based on the selected attributes associated to each gene.
#' @param modLL list-of-list containing igraph objects, annotation enrichment and ranks for each module.
#' @param orgDB Organism specific annotation library default:'org.Hs.eg.db'.
#' @return Data Frame with gene rank, score and module memebership information.
#' @examples
#' \dontrun{
#' get_master_gene_table <- function(ig, rankedGenes, modLL, orgDB="org.Hs.eg.db")
#' }
#' @keywords internal
#' @export
get_master_gene_table <- function(ig, rankedGenes, modLL, orgDB="org.Hs.eg.db"){
	orgDB <- get(orgDB)
        mapped_df <- AnnotationDbi::select(orgDB, keys=rankedGenes, columns=c("ENSEMBL","ENTREZID","GENENAME"), keytype="SYMBOL")
        #print(dim(mapped_df))
        #print(head(mapped_df))
        tmpDF <- cbind(GENE_RANK=sapply(mapped_df$SYMBOL, function(x) which(rankedGenes %in% x)), mapped_df, SCORE=sapply(mapped_df$SYMBOL, function(x) V(ig)$score[which(V(ig)$name %in% x)]), P_VALUE=sapply(mapped_df$SYMBOL, function(x) V(ig)$pval[which(V(ig)$name %in% x)]), LFC=sapply(mapped_df$SYMBOL, function(x) V(ig)$lfc[which(V(ig)$name %in% x)]))
        mapped_df <- tmpDF[order(tmpDF$GENE_RANK),]
        print("Going to add member information...")

	mapped_df[,"MEMBERSHIP"] <- NA
	geneLL <- modLL[["names"]]
	for(i in 1:length(geneLL)){
		modName <- paste0("mod_", i)
		idx <- which(mapped_df$SYMBOL %in% geneLL[[i]])
		if(length(idx)>0){
			idx1 <- which(is.na(mapped_df[idx, "MEMBERSHIP"]))
			idx2 <- which(!is.na(mapped_df[idx, "MEMBERSHIP"]))
			
			if(length(idx1)>0){
				tmpIdx <- idx[idx1]
				mapped_df[tmpIdx, "MEMBERSHIP"] <- modName
			}
			
			if(length(idx2)>0){
				tmpIdx <- idx[idx2]
				mapped_df[tmpIdx, "MEMBERSHIP"] <- paste(mapped_df[tmpIdx, "MEMBERSHIP"], modName, sep=",")
			}
		}
	}
	idx <- which(is.na(mapped_df[, "MEMBERSHIP"]))
	if(length(idx)>0){
		mapped_df[idx, "MEMBERSHIP"] <- "orphan"
	}
        return(mapped_df)
}

#' Get module gene tables.
#'
#' Get gene tables for modules in the list-of-list object containing gene rank and scores.
#'
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @importFrom AnnotationDbi select get
#'
#' @param ig igraph object representing the main graph from which the subgraph should be extracted.
#' @param rankedGenes vector of genes ordered by rank based on the selected attributes associated to each gene.
#' @param modLL list-of-list containing igraph objects, annotation enrichment and ranks for each module.
#' @param orgDB Organism specific annotation library default:'org.Hs.eg.db'.
#' @return List of data frames per modules cotaining gene rank and scores.
#' @examples
#' \dontrun{
#' get_mod_gene_tables(ig, rankedGenes, modLL, orgDB="org.Hs.eg.db")
#' }
#' @keywords internal
#' @export
get_mod_gene_tables <- function(ig, rankedGenes, modLL, orgDB="org.Hs.eg.db"){
	orgDB <- get(orgDB)
	mapped_df_ll <- list()
	geneLL <- modLL[["names"]]
	modNames <- names(geneLL)
	for(i in 1:length(geneLL)){
		#modName <- paste0("mod_", i)
		modName <- modNames[i]
		print("modName:")
		print(modName)
		
		mapped_df <- AnnotationDbi::select(orgDB, keys=geneLL[[i]], columns=c("ENSEMBL","ENTREZID","GENENAME"), keytype="SYMBOL")
		tmpDF <- cbind(GENE_RANK=sapply(mapped_df$SYMBOL, function(x) which(rankedGenes %in% x)), mapped_df, SCORE=sapply(mapped_df$SYMBOL, function(x) V(ig)$score[which(V(ig)$name %in% x)]), PValue=sapply(mapped_df$SYMBOL, function(x) V(ig)$pval[which(V(ig)$name %in% x)]), LFC=sapply(mapped_df$SYMBOL, function(x) V(ig)$lfc[which(V(ig)$name %in% x)]))
		
		mapped_df <- tmpDF[order(tmpDF$GENE_RANK),]
		mapped_df_ll[[modName]] <- mapped_df
	}
        return(mapped_df_ll)
}

#' Resolve multiple attributes formed by creating union of igraph modules.
#'
#' Resolve multiple attributes formed by creating union of igraph modules.
#'
#' @importFrom igraph as_data_frame vertex_attr delete_vertex_attr edge_attr delete_edge_attr
#'
#' @param ig_union igraph object representing the main graph from which the subgraph should be extracted.
#' @param vertex_attr_names vector of names for vertex attributes from singular igraph to keep in union igraph.
#' @param edge_attr_names vector of names for edge attributes to keep.
#' @return List of data frames per modules cotaining gene rank and scores.
#' @examples
#' \dontrun{
#' resolve_ig_union_attrs(ig_union, vertex_attr_names, edge_attr_names)
#' }
#' @keywords internal
#' @export
resolve_ig_union_attrs <- function(ig_union, vertex_attr_names, edge_attr_names){
	#Resolve Vertex Attributes
	ig_union_df <- igraph::as_data_frame(ig_union, what="vertices")
	colIdxNames <- NULL
	
	#print(paste0("Vertex Attrs : ", vertex_attr_names))
	for(attr in vertex_attr_names){
		#print(paste0("Resolving Attr : ", attr))
		colIdx <- grep(paste0("^", attr, ".*"), colnames(ig_union_df))
		colIdxNames <- c(colIdxNames, colnames(ig_union_df)[colIdx])
		#print(paste0("Cols : ", colIdxNames))
		resAttrVal <- apply(ig_union_df[,colIdx], 1, function(x){idx<-which(!is.na(x));if(length(idx)>0){res<-x[idx[1]]}else{res<-NA};return(res)})
		igraph::vertex_attr(ig_union, name=attr) <- resAttrVal
	}
	#print(list.vertex.attributes(ig_union))
	
	#Remove redundant attributes
	#print("Will Remove Following:")
	#print(colIdxNames)
	
	for(nm in colIdxNames){
		#print(paste0("Removing Attr : " , nm))
		ig_union <- igraph::delete_vertex_attr(ig_union, nm)
	}
	
	#Resolve Edge Attributes
	ig_union_df <- igraph::as_data_frame(ig_union, what="edges")
	colIdxNames <- NULL
	
	for(attr in edge_attr_names){
		colIdx <- grep(paste0("^", attr, ".*"), colnames(ig_union_df))
		colIdxNames <- c(colIdxNames, colnames(ig_union_df)[colIdx])
		resAttrVal <- apply(ig_union_df[,colIdx], 1, function(x){idx<-which(!is.na(x));if(length(idx)>0){res<-x[idx[1]]}else{res<-NA};return(res)})
		igraph::edge_attr(ig_union, name=attr) <- resAttrVal
	}
	#print(list.edge.attributes(ig_union))
	
	#Remove redundant attributes
	for(nm in colIdxNames){
		print(paste0("Removing Attr : " , nm))
		ig_union <- igraph::delete_edge_attr(ig_union, nm)
	}
	return(ig_union)
}

#Set shiny reactive vars and environment to global variables
set_shiny_vars <- function(shinyVars, shinyEnv=NULL){
        shiny_app_vars <<- shinyVars
        if(!is.null(shinyEnv)){
                if(is.environment(shinyEnv)){
                        shiny_app_env <<- shinyEnv
                }
        }
}
