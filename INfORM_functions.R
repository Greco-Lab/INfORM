
#source("INfORM_class.R")
library(grid)
library(igraph)
library(TopKLists)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(GO.db)
library(plyr)
library(GSEABase)
library(GOSemSim)
library(treemap)
library(abind)
library(minet)
library(foreach)
library(parallel)
library(doParallel)
library(ggplot2)
library(visNetwork)

net_attr <- c("betweenness", "cc", "degree", "eccentricity", "closeness", "eigenvector")

methods <- c("clr","aracne","mrnet","mrnetb") # i
est.opt <- c("pearson","spearman","kendall","mi.empirical","mi.mm","mi.shrink","mi.sg")# j
disc.opt <- c("none","equalfreq","equalwidth","globalequalwidth") # k

combineList <- function(...){
	myParam <- list(...)
	#c(myParam)
	myParam
}

utils::globalVariables(names=c("GO.db"))

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
#' @param ncores Number of cores for running instances of MINET in parallel default:2.
#' @param debug_output Print help and status messages to help debug the running of the function default:FALSE.
#' @param updateProgress Shiny application can request for update of progress from this function default:NULL.
#' @return A binary symmetrix matrix representing the median of mutual information correlation computed across various MINET combinations
#' @examples
#' \dontrun{
#' calculate_correlation_matrix(gx_table=gene_expression.df, iMethods=c("clr","aracne","mrnet","mrnetb"), iEst=c("pearson","spearman","kendall","mi.empirical","mi.mm","mi.shrink","mi.sg"), iDisc=c("none","equalfreq","equalwidth","globalequalwidth"), ncores=12, debug_output=TRUE)
#' }
#' @keywords internal
#' @export
calculate_correlation_matrix <- function(gx_table, iMethods, iEst, iDisc, ncores=2, debug_output=FALSE, updateProgress=NULL){
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
	utils::capture.output(print(paste0("Starting Parallel Cluster For ", length(parList), " Combinations ", timestamp())), file="minet-log.txt", append=T)
	utils::capture.output(print(paste0("Starting Parallel Cluster For ", length(parList), " Combinations ", timestamp())), file="minet-to-run.txt", append=T)
	utils::capture.output(print(paste0("Starting Parallel Cluster For ", length(parList), " Combinations ", timestamp())), file="minet-completed.txt", append=T)
	#out.tmp.list <- foreach(i=1:length(parList)) %do% {
	out.tmp.list <- foreach::foreach(i=1:length(parList)) %dopar% {
		#capture.output(print("Start For Each"), file="minet-log.txt", append=T)
		mt <- parList[[i]]$mt
		est <- parList[[i]]$est
		disc <- parList[[i]]$disc
		#print(parList[[i]]$mt)

		if((est == "mi.empirical" | est == "mi.mm" | est == "mi.shrink" | est == "mi.sg") & disc == "none") {
			miMat <- -1
			miMatName <- "np"
		}
		else{
			if(debug_output==TRUE)
			utils::capture.output(print(paste("----",mt,est,disc, sep="__")), file="minet-log.txt", append=T)

			utils::capture.output(print(paste0("Iteration-", i, ": ", mt, "-", est, "-", disc)), file="minet-to-run.txt", append=T)
			#print("Before Updating text")
			#if (is.function(updateProgress)){
			#	text <- paste("MINET: ", mt, "-", est, "-", disc, "-", sep="")
			#	print("Updating text")
			#	updateProgress(detail = text)
			#}

			ptm <- proc.time()
			miMat <- minet::minet(stdGX, method=mt, estimator=est, disc=disc)
			utils::capture.output(print(paste0("Iteration-", i, ", ", mt, "-", est, "-", disc, ": ", "MINET Execution Time - ", round(proc.time() - ptm)[3], " sec")), file="minet-completed.txt", append=T)
			#capture.output(print(proc.time() - ptm), file="minet-log.txt", append=T)
			miMatName <- paste(mt,est,disc,sep="__")
		}
		out.list <- list("mat"=miMat, "name"=miMatName)
	}
	utils::capture.output(print("For Each Finished, Stopping Cluster..."), file="minet-log.txt", append=T)
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
            print("Making Median Matrix...")
            arr1<-abind::abind(llGX,along=3)
            matGX <- apply(arr1,c(1,2),median)
        }else{
            print("Only one matrix, not computing median!")
            matGX <- llGX[[1]]
        }

        print("Return Matrix")
        return(matGX)
}

#' Create binary inference matrix by combining information from different inference algorithms with the help of internal function calculate_correlation_matrix().
#'
#' get_ranked_consensus_binary_matrix uses the internal function calculate_correlation_matrix() to get a single consensus matrix per inference algorithm. User can specify
#' the inference algorithms, correlation calculation methods and discretization methods, a combination of parameters will be created per inference algorithm to run calculate_correlation_matrix(), 
#' this will generate a consensus matrix per inference algorithm. The consesus matrices from different inference algorithms are used to create a single binary matrix by rank based selection of edges.
#'
#' @importFrom TopKLists Borda
#'
#' @param gx_table Gene expression table as a data frame.
#' @param iMethods Vector of valid inference algorithms for MINET package.
#' @param iEst Vector of valid correlation methods for MINET package.
#' @param iDisc Vector of valid discretization methods for MINET package.
#' @param ncores Number of cores for running instances of MINET in parallel default:2.
#' @param debug_output Print help and status messages to help debug the running of the function default:FALSE.
#' @param updateProgress Shiny application can request for update of progress from this function default:NULL.
#' @return A binary symmetrix matrix representing the edge rank based consensus from different inference algorithms
#' @examples
#' \dontrun{
#' get_ranked_consensus_binary_matrix(gx_table=gene_expression.df, iMethods=c("clr","aracne","mrnet","mrnetb"), iEst=c("pearson","spearman","kendall","mi.empirical","mi.mm","mi.shrink","mi.sg"), iDisc=c("none","equalfreq","equalwidth","globalequalwidth"), ncores=12, debug_output=TRUE)
#' }
#' @keywords internal
#' @export
get_ranked_consensus_binary_matrix <- function(gx_table, iMethods, iEst, iDisc, ncores=2, debug_output=FALSE, updateProgress=NULL){
	mat_ll <- list()
    	ranked_edges_ll <- list()

	mthdCount <- 1
	totalMthds <- length(iMethods)+1
	for(mthd in iMethods){
		if (is.function(updateProgress)) {
			text <- paste0("'", mthd, "' Median")
			value <- mthdCount / totalMthds
			updateProgress(detail = text, value = value)
		}
		mthdCount <- mthdCount + 1

		print(paste0("Calculate correlation matrix for method : ", mthd))
		mat_ll[[mthd]] <- calculate_correlation_matrix(gx_table=gx_table, iMethods=mthd, iEst=iEst, iDisc=iDisc, ncores=ncores)

		print(paste0("Get ranked edges for method : ", mthd))
		mat_ll[[mthd]][lower.tri(mat_ll[[mthd]], diag=T)] <- NA
			
		edge_df <- as.data.frame(as.table(mat_ll[[mthd]]))
		edge_df <- edge_df[-which(is.na(edge_df$Freq)),]
		edge_df <- data.frame(edge=paste0(edge_df$Var1,";",edge_df$Var2), weight=edge_df$Freq, stringsAsFactors=F)

		ranked_edges_ll[[mthd]] <- edge_df[order(edge_df$weight, decreasing=T), "edge"]
	}

        if (is.function(updateProgress)) {
	    updateProgress(detail = "Consensus Binary", value = 1)
        }

	print("Perform Borda on list of list of ranked edges.")
	borda_res <- TopKLists::Borda(ranked_edges_ll)

	print("Get a consensus binary matrix by selecting the most significant ranked edges from median rank of Borda result.")
	bin_mat <- mat_ll[[1]]
	bin_mat[,] <- 0

	input_genes <- dim(gx_table)[1]
	genes <- NULL
	total_genes <- 0
	median_list <- borda_res$TopK$median
	for(i in c(1:length(median_list))){
		if(total_genes<input_genes){
			local_genes <- strsplit(median_list[i], ";")[[1]]
			bin_mat[local_genes[1],local_genes[2]] <- 1
			bin_mat[local_genes[2],local_genes[1]] <- 1
			genes[local_genes[1]] <- 1
			genes[local_genes[2]] <- 1
			total_genes <- length(genes)
		}else{
			break
		}
	}
	print("Binary matrix computed, returning!")
	return(bin_mat)
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
#' @importFrom igraph vertex_attr
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
#' set_vertex_color(iGraph=inferred.igraph, gx_data_table=gene_expression.df, dgx_table=differential_gene_expression.df, pos_cor_color=color.pos, pos_cor_highlight_color=color.pos.highlight, neg_cor_color=color.neg, neg_cor_highlight_color=color.neg.highlight, pos_perc=pos.threshold, neg_perc=neg.threshold)
#' }
#' @keywords internal
#' @export
set_vertex_color <- function(iGraph, gx_data_table, dgx_table, pos_cor_color="salmon", pos_cor_highlight_color="red", neg_cor_color="lightblue", neg_cor_highlight_color="royalblue", pos_perc=0.95, neg_perc=0.05){
        dgx_table <- as.matrix(dgx_table)
        posCor <- rownames(dgx_table)[which(dgx_table[,1]>= stats::quantile(as.vector(dgx_table), pos_perc))]
        negCor <- rownames(dgx_table)[which(dgx_table[,1]<= stats::quantile(as.vector(dgx_table), neg_perc))]

        col_length <- length(rownames(dgx_table))
        color_vector <- rep("lightgrey", col_length)
        highlight_color_vector <- rep("darkgrey", col_length)
        print("For Loop for populating vertex color vector...")
        for(i in c(which(is.element(rownames(gx_data_table),posCor)))){
                color_vector[i] <- pos_cor_color
                highlight_color_vector[i] <- pos_cor_highlight_color
        }

        for(i in c(which(is.element(rownames(gx_data_table),negCor)))){
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
#' set_edge_color(iGraph=inferred.igraph, gx_data_table=gene_expression.df, pos_cor_color=color.pos, pos_cor_highlight_color=color.pos.highlight, neg_cor_color=color.neg, neg_cor_highlight_color=color.neg.highlight)
#' }
#' @keywords internal
#' @export
set_edge_color <- function(iGraph, gx_data_table, pos_cor_color="salmon", pos_cor_highlight_color="red", neg_cor_color="lightblue", neg_cor_highlight_color="royalblue"){
	tGX <- t(gx_data_table)
	corGX <- cor(tGX)

	matColor <- ifelse(corGX<=0, neg_cor_color, pos_cor_color)
	matHighlightColor <- ifelse(corGX<=0, neg_cor_highlight_color, pos_cor_highlight_color)

	edgeIDs <- igraph::get.edges(iGraph, igraph::E(iGraph))

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
#' get_ranked_gene_list(iGraph=inferred.igraph, rank_list_attr=c("betweenness", "cc", "degree", "eccentricity", "closeness", "eigenvector", "score"), debug_output=FALSE)
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

#make_sub_ggplot <- function(candidates){
#        "Make ggplot object from sub iGraph."
#        local_igraph <- info_subGraph
#
#        print("Getting sub ggplot...")
#        print(igraph::list.vertex.attributes(local_igraph))
#        info_subGGplot <<- get_ggplot(local_igraph)
#}

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
	selectDF <- select(annDB, keys=genelist, columns=keyType, keytype="SYMBOL")
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
	selectAnnDF <- select(annDB, keys=resDF$ID, columns="ENTREZID", keytype=keyType)
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
	selectAnnDF <- select(GO.db, keys=as.vector(all_fltr_union_score$ID), columns=c("TERM", "ONTOLOGY"), keytype="GOID")
	all_fltr_union_score <-  cbind(selectAnnDF, all_fltr_union_score[,c(2,3)])
	all_fltr_union_score
}

#get_go_freq <- function(annDB="org.Hs.eg.db"){
#	annDB <- get(annDB)
#	selectAnnDF <- select(annDB, keys=keys(annDB, keytype="GO"), columns="ENTREZID", keytype="GO")
#	colnames(selectAnnDF)[1] <- "ID"
#	go_map_DF <- plyr::ddply(selectAnnDF,.(ID),nrow)
#	colnames(go_map_DF) <- c("ID", "ALL_Genes")
#	genome.len <- length(unique(keys(annDB, keytype="ENTREZID")))
#	go_map_DF$Freq <- round((go_map_DF$ALL_Genes/genome.len)*100,2)
#	return(go_map_DF)
#}

#get_semsim <- function(input_GO_DF){
#	out_semsim_table <- list()
#
#	BP_semsim <- mgoSim(input_GO_DF[input_GO_DF$ONTOLOGY == "BP", 1], input_GO_DF[input_GO_DF$ONTOLOGY == "BP", 1],  measure="Rel", ont="BP", combine=NULL)
#	BP_semsim <- round(BP_semsim, 2)
#	out_semsim_table[["BP"]] <- as.data.frame(as.table(BP_semsim))
#	out_semsim_table[["BP"]] <- out_semsim_table[["BP"]][which(upper.tri(BP_semsim, diag=FALSE)),]
#	colnames(out_semsim_table[["BP"]]) <- c("ID1", "ID2", "SIM")
#
#	CC_semsim <- mgoSim(input_GO_DF[input_GO_DF$ONTOLOGY == "CC", 1], input_GO_DF[input_GO_DF$ONTOLOGY == "CC", 1],  measure="Rel", ont="CC", combine=NULL)
#	CC_semsim <- round(CC_semsim, 2)
#	out_semsim_table[["CC"]] <- as.data.frame(as.table(CC_semsim))
#	out_semsim_table[["CC"]] <- out_semsim_table[["CC"]][which(upper.tri(CC_semsim, diag=FALSE)),]
#	colnames(out_semsim_table[["CC"]]) <- c("ID1", "ID2", "SIM")
#
#	MF_semsim <- mgoSim(input_GO_DF[input_GO_DF$ONTOLOGY == "MF", 1], input_GO_DF[input_GO_DF$ONTOLOGY == "MF", 1],  measure="Rel", ont="MF", combine=NULL)
#	MF_semsim <- round(MF_semsim, 2)
#	out_semsim_table[["MF"]] <- as.data.frame(as.table(MF_semsim))
#	out_semsim_table[["MF"]] <- out_semsim_table[["MF"]][which(upper.tri(MF_semsim, diag=FALSE)),]
#	colnames(out_semsim_table[["MF"]]) <- c("ID1", "ID2", "SIM")
#
#	return(out_semsim_table)
#}

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
#' @return List of data frames containing clustered GO information, with major terms and individual terms. One data frame for each GO type BP, CC and MF.
#' @examples
#' \dontrun{
#' go_summarization(enriched_GO_DF, annDB="org.Hs.eg.db", simMeasure="Rel", treeHeight=9)
#' }
#' @keywords internal
#' @export
go_summarization <- function(enriched_GO_DF, annDB="org.Hs.eg.db", simMeasure="Rel", treeHeight=9){
	GODB_vals = select(GO.db, keys(GO.db, "GOID"), c("TERM", "ONTOLOGY"))
	GODB_vals_ll <- list()
	GODB_vals_ll[["BP"]] <- GODB_vals[GODB_vals$ONTOLOGY=="BP",c("GOID", "TERM")]
	GODB_vals_ll[["CC"]] <- GODB_vals[GODB_vals$ONTOLOGY=="CC",c("GOID", "TERM")]
	GODB_vals_ll[["MF"]] <- GODB_vals[GODB_vals$ONTOLOGY=="MF",c("GOID", "TERM")]

	simClust_DF <- list()
	for(typ in c("BP", "CC", "MF")){
		print(paste0("###### Clustering For Ontology - ", typ, ' ######'))
		ont <- enriched_GO_DF[enriched_GO_DF$ONTOLOGY == typ,]
		ont_ids <- as.character(ont[which(ont$GOID %in% GODB_vals_ll[[typ]]$GOID), "GOID"])

		#simsem_matrix <- GOSemSim::mgoSim(ont_ids, ont_ids,  measure="Rel", ont=typ, organism="human", combine=NULL)
                d <- godata(annDB, ont=typ, computeIC=TRUE)
		simsem_matrix <- GOSemSim::mgoSim(ont_ids, ont_ids, semData=d, measure="Rel", combine=NULL)
		dim_len <- dim(simsem_matrix)[1]
		print(dim_len)
		simsem_matrix_clean <- simsem_matrix[rowSums(is.na(simsem_matrix))!=dim_len, colSums(is.na(simsem_matrix))!=dim_len]
		simsem_matrix <- round(simsem_matrix_clean, 2)
		
		simClust <- stats::hclust(stats::as.dist(1-simsem_matrix))

		cut_height <- 0.9
		simClust_cut <- stats::cutree(simClust, h = cut_height)
		simClust_DF[[typ]] <- data.frame(ID=character(), TERM=character(), Score=integer(), Representative=character(), stringsAsFactors=F)
		
		for(i in unique(simClust_cut)){
			memberIDs <- names(simClust_cut[simClust_cut == i])
			#clustID <- paste0(typ, "_", i)
			selRow <- ont[which(ont$GOID %in% memberIDs), c("TERM", "Score")]
                        clustID <- selRow[which.max(selRow$Score), "TERM"]
			#clustID <- selRow[which(selRow$Score==max(selRow$Score)), "TERM"]
			localDF <- data.frame(ID=memberIDs, TERM=selRow$TERM, Score=selRow$Score, Representative=clustID, stringsAsFactors=F)
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

#get_ggplot = function(iG, plot_layout="with_fr"){
#        "Get ggplot object representing the input iGraph."
#        local_igraph <- iG
#
#        #iG.fr <- layout.fruchterman.reingold(local_igraph)
#        #layout_format = "with_sphere"
#        layout_format <- plot_layout
#        iG.fr <- layout_(local_igraph, layout=get(layout_format)())
#        iG.fr.df <- as.data.frame(iG.fr)
#        iG.fr.df$ver <- igraph::V(local_igraph)$name
#        iG.fr.df$verColor <- igraph::V(local_igraph)$color
#        gFrame <- get.data.frame(local_igraph)
#        gFrame$from.x <- iG.fr.df$V1[match(gFrame$from, iG.fr.df$ver)]
#        gFrame$from.y <- iG.fr.df$V2[match(gFrame$from, iG.fr.df$ver)]
#        gFrame$to.x <- iG.fr.df$V1[match(gFrame$to, iG.fr.df$ver)]
#        gFrame$to.y <- iG.fr.df$V2[match(gFrame$to, iG.fr.df$ver)]
#        #cPalette <- unique(sort(V(local_igraph)$color))
#        cPalette <- unique(sort(c(igraph::V(local_igraph)$color, igraph::V(local_igraph)$highlightcolor, igraph::E(local_igraph)$color, igraph::E(local_igraph)$highlightcolor)))
#        names(cPalette) <- cPalette
#        #ePalette <- unique(sort(c(E(local_igraph)$color, E(local_igraph)$highlightcolor)))
#
#        local_ggplot <- ggplot2::ggplot() +
#        geom_segment(data=gFrame,aes(x=from.x,xend=to.x,y=from.y,yend=to.y,colour=color)) +
#        #scale_colour_manual(values=cPalette) +
#        scale_colour_manual(values=ePalette) +
#        geom_point(data=iG.fr.df,aes(x=V1,y=V2),size=21,colour="black") +
#        geom_point(data=iG.fr.df,aes(x=V1,y=V2,colour=verColor),size=20) +
#        #scale_colour_manual(values=cPalette) +
#        scale_colour_manual(name="", values=cPalette) +
#        geom_text(data=iG.fr.df,aes(x=V1,y=V2,label=ver),size=3) +
#        scale_x_continuous(expand=c(0,1))+
#        scale_y_continuous(expand=c(0,1))+
#        theme_bw()+
#        theme(
#                axis.text.x=element_blank(),
#                axis.text.y=element_blank(),
#                axis.ticks=element_blank(),
#                axis.title.x=element_blank(),
#                axis.title.y=element_blank(),
#                legend.position="none",
#                panel.grid.major=element_blank(),
#                panel.grid.minor=element_blank(),
#                panel.background=element_blank()
#        )
#
#        write.table(iG.fr.df, file="iG.fr.df.txt", sep="\t", quote=FALSE)	
#        local_ggplot
#}

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
#' get_visNetwork(iG, plot_layout="nicely", vBorderColor="black", vShape="circle", vFontColor="#343434", vSize=50, eWidth=10, degDepth=2)
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

