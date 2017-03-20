library(shiny)
library(shinyjs)
library(shinyBS)
library(shinydashboard)
library(DT)
library(R.utils)
library(ggplot2)
library(GSEABase)
library(visNetwork)
source("~/git/INfORM/INfORM_functions.R")

options(shiny.maxRequestSize=500*1024^2)

net_attr <- c("betweenness", "cc", "degree", "closeness", "eigenvector", "score")
rank_attr <- c("betweenness", "cc", "degree", "closeness", "eigenvector", "score")

methods <- c("clr","aracne","mrnet","mrnetb") # i
est.opt <- c("pearson","spearman","kendall","mi.empirical","mi.mm","mi.shrink","mi.sg")# j
disc.opt <- c("none","equalfreq","equalwidth","globalequalwidth") # k

gxTable <- NULL
gxCorMat <- NULL
tmpMethod <- NULL

shinyServer(
	function(input, output, session){

		#myValues <- shiny::reactiveValues(gxTable=NULL, dgxTable=NULL, method=methods, est=est.opt, disc=disc.opt, cores=2, rankAttr=rank_attr, gxCorMat=NULL, updateTabCorMat=0, cList=NULL, mList=NULL, GO_clust_summ_list=NULL)
		myValues <- shiny::reactiveValues(gxTable=NULL, dgxTable=NULL, gxCorMat=NULL, iGraph=NULL, updateTabCorMat=0, cList=NULL, mList=NULL, GO_clust_summ_list=NULL)

		myValues$gxTable <- shiny::reactive({
			gxFile <- input$gx
			if (is.null(gxFile) || is.null(myValues$dgxTable()))
			return(NULL)

			print("Reading loaded gx File...")
			gx <- read.csv(gxFile$datapath, row.names=1, header=TRUE, sep="\t", check.names=FALSE)
                        #gx <- gx[order(rownames(gx)),]
                        gx <- gx[rownames(myValues$dgxTable()),]
                        gx
		})

		myValues$dgxTable <- shiny::reactive({
			print("Checking dgx File...")
			dgxFile <- input$dgx
			if (is.null(dgxFile))
			return(NULL)

			print("Reading loaded dgx File...")
			dgx <- read.csv(dgxFile$datapath, row.names=1, header=FALSE, sep="\t", check.names=FALSE)
                        dgx <- dgx[order(rownames(dgx)),,drop=F]
                        print(str(dgx))
                        dgx
		})

                myValues$method <- shiny::reactive({
                        if(is.null(input$method))
                        return(methods)

                        input$method
                })

                myValues$est <- shiny::reactive({
                        if(is.null(input$est))
                        return(est.opt)

                        input$est
                })

                myValues$disc <- shiny::reactive({
                        if(is.null(input$disc))
                        return(disc.opt)

                        input$disc
                })

                myValues$cores <- shiny::reactive({
                        if(is.null(input$cores))
                        return(cores<-2)

                        input$cores
                })

                myValues$rankAttr <- shiny::reactive({
                        if(is.null(input$rankAttr))
                        return(rank_attr)

                        input$rankAttr
                })

		observeEvent(input$runINfORM, {
			print("Run INfORM.....")

                        #shiny::updateTabsetPanel(session, "display", selected="net_display")

                        shiny::validate(
				need(!is.null(myValues$gxTable()), "No Gene Expression Table Provided!")
			)
                        shiny::validate(
				need(!is.null(myValues$dgxTable()), "No Differential Expression Table Provided!")
			)
			shiny::validate(
				need(!is.null(myValues$method()), "No Inference Algorithm Selected!")
			)
			shiny::validate(
				need(!is.null(myValues$est()), "No Correlation Selected!")
			)
			shiny::validate(
				need(!is.null(myValues$disc()), "No Discretization Method Selected!")
			)
			print("Passed Validations.....")

                        progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}
                        on.exit(progress$close())

			localGxTable <- myValues$gxTable()
			localDgxTable <- myValues$dgxTable()

                        corMat <- NULL
                        matFile <- input$mat
			if(is.null(matFile)){
                                myMethod <- myValues$method()
                                myEst <- myValues$est()
                                myDisc <- myValues$disc()
                                myCores <- as.numeric(myValues$cores())

                                #print(paste0("methods : ", myMethod))
                                #print(paste0("ESTs : ", myEst))
                                #print(paste0("Discs : ", myDisc))

                                progress$set(message="Inference", value=0)

                                corMat <- get_ranked_consensus_binary_matrix(localGxTable, iMethods=myMethod, iEst=myEst, iDisc=myDisc, ncores=myCores, debug_output=FALSE, updateProgress=updateProgress)
                        }else{
                                print("Reading loaded correlation matrix...")
                                myValues$updateTabCorMat <- 1
                                corMat <- as.matrix(read.csv(matFile$datapath, row.names=1, header=TRUE, sep="\t", check.names=FALSE))
                        }

			localGxCorMat <- corMat

                        myValues$totalGenes <- nrow(localGxCorMat)

			##Clean Correlation Matrix
			#if(input$checkUnconnected)
			#{
			#	progress$set(message="Cleaning Adjacency Matrix", value=0)

			#	print("Clean Correlation Matrix...")
			#	connectedGenes <- unique(sort(unlist(apply(localGxCorMat, 1, function(x){rNames<-names(x); idxs<-which(x==1); rNames[idxs]}))))
			#	print(paste0(length(connectedGenes), " Connected Genes"))
			#	localGxCorMat <- localGxCorMat[connectedGenes,connectedGenes]

			#	print("Clean GX Table...")
			#	print(dim(localGxTable))
			#	whichVec <- which(rownames(localGxTable) %in% connectedGenes)
			#	print(length(whichVec))
			#	localGxTable <- localGxTable[whichVec,]
			#	print("GX Table Cleaned...")
			#	print(dim(localGxTable))

			#	print("Clean DGX Table...")
			#	print(dim(localDgxTable))
			#	whichVec <- which(rownames(localDgxTable) %in% connectedGenes)
			#	print(length(whichVec))
			#	localDgxTable <- localDgxTable[whichVec,,drop=F]
			#	print("DGX Table Cleaned...")
			#	print(dim(localDgxTable))
			#}

			print("Get IGraph Object...")
                        progress$set(message="Getting iGraph Object", value=0)
                        localIGraph <- get_iGraph(localGxCorMat)

                        print("Setting DGX values as 'score'")
                        vertex_attr(localIGraph, name="score") <- localDgxTable[V(localIGraph)$name,]

                        myValues$iGraph <- localIGraph
                        rankAttr <- myValues$rankAttr()

                        progress$set(message="Getting Gene Ranks", value=0)
                        localRankedGenes <- get_ranked_gene_list(localIGraph, rank_list_attr=rankAttr, debug_output=FALSE)
                        myValues$rankedGenes <- localRankedGenes

			progress$set(message="Getting Candidate Genes", value=0)
                        if(input$candid_type == "count"){
                                candidates <- localRankedGenes[c(1:as.integer(input$candid_val))]
                        }
                        else{
                                countForPerc <- round(as.integer(input$candid_val)*length(localRankedGenes)/100)
                                candidates <- localRankedGenes[c(1:countForPerc)]
                        }

			print("Selected Candidates...")
			print(candidates)
			myValues$candidates <- candidates

                        progress$set(message="Getting Candidates...", value=0)
			updateProgress(detail=paste0("Top: ",input$l1_val), value=1/3)
			l1_candidates <- localRankedGenes[c(1:as.integer(input$l1_val))]

			updateProgress(detail=paste0("Top: ",input$l2_val), value=2/3)
			l2_candidates <- localRankedGenes[c(1:as.integer(input$l2_val))]

			updateProgress(detail=paste0("Top: ",input$l3_val), value=3/3)
			l3_candidates <- localRankedGenes[c(1:as.integer(input$l3_val))]

			print("Candidates1:")
			print(l1_candidates)
			print("Candidates2:")
			print(l2_candidates)
			print("Candidates3:")
			print(l3_candidates)

			myValues$cList <- list(l1_candidates, l2_candidates, l3_candidates)

			progress$set(message="Getting Sub-Graphs...", value=0)
			updateProgress(detail=paste0("Top: ",input$l1_val), value=1/3)
			l1_subGraph <- get_sub_graph(localIGraph, l1_candidates, vertex_level=1, shortest_paths=TRUE)
			myValues$smallSubPlot <- l1_subGraph

			updateProgress(detail=paste0("Top: ",input$l2_val), value=2/3)
			l2_subGraph <- get_sub_graph(localIGraph, l2_candidates, vertex_level=1, shortest_paths=TRUE)
			myValues$mediumSubPlot <- l2_subGraph

			updateProgress(detail=paste0("Top: ",input$l3_val), value=3/3)
			l3_subGraph <- get_sub_graph(localIGraph, l3_candidates, vertex_level=1, shortest_paths=TRUE)
			myValues$largeSubPlot <- l3_subGraph

			l1_genelist <- get.vertex.attribute(l1_subGraph, "name")
			l2_genelist <- get.vertex.attribute(l2_subGraph, "name")
			l3_genelist <- get.vertex.attribute(l3_subGraph, "name")

			print("subgraph1:")
			print(l1_genelist)
			print("subgraph2:")
			print(l2_genelist)
			print("subgraph3:")
			print(l3_genelist)

			localSubGraphGeneLists <- list(l1_genelist, l2_genelist, l3_genelist)
			myValues$mList <- localSubGraphGeneLists

			#rOrgDB <- unlist(strsplit(input$organism, ";"))[2]
			rOrgDB <- unlist(strsplit(input$organism, ";"))[1]
			progress$set(message="Getting Enriched Annotation...", value=0)
			updateProgress(detail=paste0("Top: ",input$l1_val), value=1/3)
			l1_enrichmentDF <- annotation_enrichment(l1_genelist, annDB=rOrgDB)
			updateProgress(detail=paste0("Top: ",input$l2_val), value=2/3)
			l2_enrichmentDF <- annotation_enrichment(l2_genelist, annDB=rOrgDB)
			updateProgress(detail=paste0("Top: ",input$l3_val), value=3/3)
			l3_enrichmentDF <- annotation_enrichment(l3_genelist, annDB=rOrgDB)

			progress$set(message="Getting Progressive Enrichment...", value=0)
			updateProgress(detail="Processing: ", value=0.5)
			localEnrichmentDF <- progressive_enrichment(localSubGraphGeneLists, l1_enrichmentDF, l2_enrichmentDF, l3_enrichmentDF)
			myValues$enrichmentDF <- localEnrichmentDF
			updateProgress(detail="Obtained: ", value=1)

			progress$set(message="Creating Tile Plots...", value=0)
			sim_cutoff=0.4

			updateProgress(detail="Clustering GO Terms", value=1/3)
                        GO_clust_summ_list <- go_summarization(localEnrichmentDF)
                        myValues$GO_clust_summ_list <- GO_clust_summ_list
			updateProgress(detail="Finishing: ", value=3/3)

			myValues$gxCorMat <- localGxCorMat
                })

		shiny::observeEvent(input$calcMat, {
			cat("INSIDE OBSERVE EVENT.....")
			shiny::updateTabsetPanel(session, "display", selected="corMat_display")
		})

		output$totalGeneBox <- shinydashboard::renderValueBox({
			if(is.null(myValues$gxCorMat)){
				geneCount <- "NA"
			}
			else{
			    	localGxCorMat <- myValues$gxCorMat
				geneCount <- nrow(localGxCorMat)
			}
			
			shinydashboard::valueBox(
				paste0(geneCount), "Total Genes", icon=icon("list")
			)
		})

		output$connGeneBox <- shinydashboard::renderValueBox({
			if(is.null(myValues$gxCorMat)){
				connectedGeneCount <- "NA"
			}
			else{
			    	localGxCorMat <- myValues$gxCorMat
				connectedGeneCount <- length(unique(sort(unlist(apply(localGxCorMat, 1, function(x){rNames<-names(x); idxs<-which(x==1); rNames[idxs]})))))
			}
                        
			shinydashboard::valueBox(
				paste0(connectedGeneCount), "Connected Genes", icon=icon("thumbs-o-up"), color="green"
			)
		})

		output$unGeneBox <- shinydashboard::renderValueBox({
			if(is.null(myValues$gxCorMat)){
				connectedGeneCount <- "NA"
				unConnectedGeneCount <- "NA"
			}
			else{
			    	localGxCorMat <- myValues$gxCorMat
				connectedGeneCount <- length(unique(sort(unlist(apply(localGxCorMat, 1, function(x){rNames<-names(x); idxs<-which(x==1); rNames[idxs]})))))
				unConnectedGeneCount <- dim(localGxCorMat)[1] - connectedGeneCount
			}
                        
			shinydashboard::valueBox(
				paste0(unConnectedGeneCount), "Un-connected Genes", icon=icon("thumbs-o-down"), color="red"
			)
		})

		output$unGeneList <- shiny::renderTable({
			shiny::validate(
                                need(!is.null(myValues$gxCorMat), "")
                        )

                        localgxcormat <- myValues$gxCorMat

			connectedGeneList <- unique(sort(unlist(apply(localgxcormat, 1, function(x){rNames<-names(x); idxs<-which(x==1); rNames[idxs]}))))
			unConnectedGeneList <- rownames(localgxcormat)[-which(rownames(localgxcormat) %in% connectedGeneList)]
			#print(unConnectedGeneList)
			tmpDF <- as.data.frame(unConnectedGeneList)
			#print(class(tmpDF))
			colnames(tmpDF) <- "Un-connected Genes"
			#print(tmpDF)
			tmpDF
		})

		output$downloadUnGenes <- shiny::downloadHandler(
			filename = function(){
				paste("Un-Connected_Genes_", Sys.Date(), '.txt', sep='')
			},
			content = function(con){
				connectedGeneList <- unique(sort(unlist(apply(myValues$loadedCorMat(), 1, function(x){rNames<-names(x); idxs<-which(x==1); rNames[idxs]}))))
	                        unConnectedGeneList <- rownames(myValues$loadedCorMat())[-which(rownames(myValues$loadedCorMat()) %in% connectedGeneList)]
				tmpDF <- as.data.frame(unConnectedGeneList)
				colnames(tmpDF) <- "Un-connected Genes"
				write.table(tmpDF, con, row.names=FALSE, quote=FALSE, sep="\t")
			}
		)

		output$downloadMat <- shiny::downloadHandler(
			filename = function(){
                                paste("Adjacency_Matrix_", Sys.Date(), '.txt', sep='')
                        },
			content = function(con){
				write.table(myValues$gxCorMat, con, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
			}
		)

		shiny::observe({
			if(is.null(input$mat)){
				shinyjs::enable("calcMat")
				shinyjs::enable("method")
				shinyjs::enable("est")
				shinyjs::enable("disc")
				shinyjs::enable("cores")
			}
			else{
				shinyjs::disable("calcMat")
				shinyjs::disable("method")
				shinyjs::disable("est")
				shinyjs::disable("disc")
				shinyjs::disable("cores")
			}

                        if(is.null(input$gx) || is.null(input$dgx)){
				shinyjs::disable("runINfORM")
                        }else{
                                shinyjs::enable("runINfORM")
                        }

			if(is.null(myValues$gxCorMat)){
				shinyjs::disable("downloadMat")
				shinyjs::disable("downloadUnGenes")
			}
			else{
				shinyjs::enable("downloadMat")
				
			    	localGxCorMat <- myValues$gxCorMat
				connectedGeneCount <- length(unique(sort(unlist(apply(localGxCorMat, 1, function(x){rNames<-names(x); idxs<-which(x==1); rNames[idxs]})))))
				unConnectedGeneCount <- dim(localGxCorMat)[1] - connectedGeneCount
				if(unConnectedGeneCount>0){
				    shinyjs::enable("downloadUnGenes")
				}
			}
		})

                myValues$iGraph_dynamic <- shiny::reactive({
                        shiny::validate(
				need(!is.null(myValues$iGraph), "No Graph to Plot!")
			)
                        #if (is.null(myValues$iGraph))
			#return(NULL)

                        progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}
                        on.exit(progress$close())

			progress$set(message="Formatting Graph", value=0)

			localGxTable <- myValues$gxTable()
			localDgxTable <- myValues$dgxTable()
                        localIGraph <- myValues$iGraph
                        localRankedGenes <- myValues$rankedGenes

                        #if(input$checkUnconnected)
			#{
			#	print("Clean Correlation Matrix...")
			#	connectedGenes <- V(localIGraph)$name
			#	print(paste0(length(connectedGenes), " Connected Genes"))

			#	print("Clean GX Table...")
			#	print(dim(localGxTable))
			#	whichVec <- which(rownames(localGxTable) %in% connectedGenes)
			#	print(length(whichVec))
			#	localGxTable <- localGxTable[whichVec,]
			#	print("GX Table Cleaned...")
			#	print(dim(localGxTable))

			#	print("Clean DGX Table...")
			#	print(dim(localDgxTable))
			#	whichVec <- which(rownames(localDgxTable) %in% connectedGenes)
			#	print(length(whichVec))
			#	localDgxTable <- localDgxTable[whichVec,,drop=F]
			#	print("DGX Table Cleaned...")
			#	print(dim(localDgxTable))
			#}

                        if(input$vColType=="score"){
                                localIGraph <- set_vertex_color(localIGraph, localGxTable, localDgxTable, pos_cor_color=input$vColP, pos_cor_highlight_color=input$vhColP, neg_cor_color=input$vColN, neg_cor_highlight_color=input$vhColN, pos_perc=as.numeric(input$posPerc), neg_perc=as.numeric(input$negPerc))
                        }else if(input$vColType=="rank"){
                                rbPal <- colorRampPalette(c(input$vColT,input$vColB))
                                rankGradient <- rbPal(vcount(localIGraph))
                                names(rankGradient) <- localRankedGenes
                                igraph::vertex_attr(localIGraph, name="color") <- rankGradient[V(localIGraph)$name]
                        }
                        localIGraph <- set_edge_color(localIGraph, localGxTable, pos_cor_color=input$eColP, pos_cor_highlight_color=input$ehColP, neg_cor_color=input$eColN, neg_cor_highlight_color=input$ehColN)
                        return(localIGraph)
                })

                output$IGraphVis <- visNetwork::renderVisNetwork({
                        shiny::validate(
                                #need(!is.null(myValues$iGraph), "No graph to plot!")
                                need(!is.null(myValues$iGraph_dynamic), "Waiting for updated graph...")
                        )

			progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}
                        on.exit(progress$close())

			progress$set(message="Plotting", value=0)

                        localIGraph <- myValues$iGraph_dynamic()

                        #mainVisNetwork <-  get_visNetwork(localIGraph, plot_layout=input$graphLayout, vBorderColor=input$vBrdCol, vShape=input$vShape, vFontColor=input$vLblCol, vSize=input$vSize, eWidth=input$eWidth, vColTop=input$vColT, vColBottom=input$vColB, vColType=input$vColType)
                        mainVisNetwork <-  get_visNetwork(localIGraph, plot_layout=input$graphLayout, vBorderColor=input$vBrdCol, vShape=input$vShape, vFontColor=input$vLblCol, vSize=input$vSize, eWidth=input$eWidth, degDepth=input$dDepth)
                        mainVisNetwork
                })

		output$downloadGraph <- shiny::downloadHandler(
			filename = function(){
				paste("network_", Sys.Date(), '.', input$graphExportFormat, sep='')
			},
			content = function(con){
				igraph::write_graph(myValues$iGraph, con, format=input$graphExportFormat)
			}
		)

		output$downloadGraphGenes <- shiny::downloadHandler(
			filename = function(){
				paste("network_genes_", Sys.Date(), '.txt', sep='')
			},
			content = function(con){
				write.table(V(myValues$iGraph)$name, con, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
			}
		)

		output$downloadGO <- shiny::downloadHandler(
			filename = function(){
				paste("enriched_GO_", Sys.Date(), '.txt', sep='')
			},
			content = function(con){
				write.table(myValues$enrichmentDF, con, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
			}
		)

		output$enrichmentDT <- DT::renderDataTable({
			shiny::validate(
                                need(!is.null(myValues$enrichmentDF), "Waiting for annotation enrichment results...")
                        )

			print("In Enrichment Data Frame, renderDataTable!")
			DT::datatable(myValues$enrichmentDF, filter = list(position='top', clear=FALSE), options = list(search = list(regex=TRUE, caseInsensitive=FALSE), scrollX=TRUE))
		})

		output$tilePlotBP <- shiny::renderPlot({
			shiny::validate(
                                need(!is.null(myValues$enrichmentDF), "Waiting for annotation enrichment results...")
                        )
			print("INSIDE TILE PLOT!")

			plot_treemap(myValues$GO_clust_summ_list$BP, ont="BP")
		})

		output$tilePlotCC <- shiny::renderPlot({
			shiny::validate(
                                need(!is.null(myValues$enrichmentDF), "Waiting for annotation enrichment results...")
                        )
			print("INSIDE TILE PLOT!")

			plot_treemap(myValues$GO_clust_summ_list$CC, ont="CC")
		})

		output$tilePlotMF <- shiny::renderPlot({
			shiny::validate(
                                need(!is.null(myValues$enrichmentDF), "Waiting for annotation enrichment results...")
                        )
			print("INSIDE TILE PLOT!")

			plot_treemap(myValues$GO_clust_summ_list$MF, ont="MF")
		})

		###Sub-Networks
                output$smallSubNetVis <- visNetwork::renderVisNetwork({
                        shiny::validate(
                                #need(!is.null(myValues$enrichmentDF), "Waiting for annotation enrichment results...")
                                need(!is.null(myValues$iGraph_dynamic), "Waiting for formatted graph...")
                        )

                        progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}
                        on.exit(progress$close())

			progress$set(message="Extracting Sub-Graph", value=0)

                        localIGraph <- myValues$iGraph_dynamic()

                        localRankedGenes <- myValues$rankedGenes
			l1_candidates <- localRankedGenes[c(1:as.integer(input$l1_val))]
			l1_subGraph <- get_sub_graph(localIGraph, l1_candidates, vertex_level=1, shortest_paths=TRUE)
                        if(input$vColType=="score"){
                                l1_subGraph <- highlight_vertices(l1_subGraph, l1_candidates)
                        }

			progress$set(message="Plotting", value=0)

                        smallVisNetwork <-  get_visNetwork(l1_subGraph, plot_layout=input$graphLayout, vBorderColor=input$vBrdCol, vShape=input$vShape, vFontColor=input$vLblCol, vSize=input$vSize, eWidth=input$eWidth, degDepth=input$dDepth)
                        smallVisNetwork
                })

                output$downloadSmallSubNet <- shiny::downloadHandler(
                        filename = function(){
                                paste("top", input$l1_val, "_subNet_", Sys.Date(), '.', input$smallSubNetExportFormat, sep='')
                        },
                        content = function(con){
                                write.graph(myValues$smallSubPlot, con, format=input$smallSubNetExportFormat)
                        }
                )

                output$downloadSmallSubNetGenes <- shiny::downloadHandler(
                        filename = function(){
                                paste("top", input$l1_val, "_subNet_genes_", Sys.Date(), '.txt', sep='')
                        },
                        content = function(con){
                                write.table(V(myValues$smallSubPlot)$name, con, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
                        }
                )

                output$mediumSubNetVis <- visNetwork::renderVisNetwork({
                        shiny::validate(
                                #need(!is.null(myValues$enrichmentDF), "Waiting for annotation enrichment results...")
                                need(!is.null(myValues$iGraph_dynamic), "Waiting for formatted graph...")
                        )

                        progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}
                        on.exit(progress$close())

			progress$set(message="Extracting Sub-Graph", value=0)

                        localIGraph <- myValues$iGraph_dynamic()

                        localRankedGenes <- myValues$rankedGenes
			l2_candidates <- localRankedGenes[c(1:as.integer(input$l2_val))]
			l2_subGraph <- get_sub_graph(localIGraph, l2_candidates, vertex_level=1, shortest_paths=TRUE)
                        if(input$vColType=="score"){
                                l2_subGraph <- highlight_vertices(l2_subGraph, l2_candidates)
                        }

			progress$set(message="Plotting", value=0)

                        mediumVisNetwork <-  get_visNetwork(l2_subGraph, plot_layout=input$graphLayout, vBorderColor=input$vBrdCol, vShape=input$vShape, vFontColor=input$vLblCol, vSize=input$vSize, eWidth=input$eWidth, degDepth=input$dDepth)
                        mediumVisNetwork
                })

                output$downloadMediumSubNet <- shiny::downloadHandler(
                        filename = function(){
                                paste("top", input$l2_val, "_subNet_", Sys.Date(), '.', input$mediumSubNetExportFormat, sep='')
                        },
                        content = function(con){
                                write.graph(myValues$mediumSubPlot, con, format=input$mediumSubNetExportFormat)
                        }
                )

                output$downloadMediumSubNetGenes <- shiny::downloadHandler(
                        filename = function(){
                                paste("top", input$l2_val, "_subNet_genes_", Sys.Date(), '.txt', sep='')
                        },
                        content = function(con){
                                write.table(V(myValues$mediumSubPlot)$name, con, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
                        }
                )

                output$largeSubNetVis <- visNetwork::renderVisNetwork({
                        shiny::validate(
                                #need(!is.null(myValues$enrichmentDF), "Waiting for annotation enrichment results...")
                                need(!is.null(myValues$iGraph_dynamic), "Waiting for formatted graph...")
                        ) 

                        progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}
                        on.exit(progress$close())

			progress$set(message="Extracting Sub-Graph", value=0)

                        localIGraph <- myValues$iGraph_dynamic()

                        localRankedGenes <- myValues$rankedGenes
			l3_candidates <- localRankedGenes[c(1:as.integer(input$l3_val))]
			l3_subGraph <- get_sub_graph(localIGraph, l3_candidates, vertex_level=1, shortest_paths=TRUE)
                        if(input$vColType=="score"){
                                l3_subGraph <- highlight_vertices(l3_subGraph, l3_candidates)
                        }

			progress$set(message="Plotting", value=0)

                        largeVisNetwork <-  get_visNetwork(l3_subGraph, plot_layout=input$graphLayout, vBorderColor=input$vBrdCol, vShape=input$vShape, vFontColor=input$vLblCol, vSize=input$vSize, eWidth=input$eWidth, degDepth=input$dDepth)
                        largeVisNetwork
                })

                output$downloadLargeSubNet <- shiny::downloadHandler(
                        filename = function(){
                                paste("top", input$l3_val, "_subNet_", Sys.Date(), '.', input$largeSubNetExportFormat, sep='')
                        },
                        content = function(con){
                                write.graph(myValues$largeSubPlot, con, format=input$largeSubNetExportFormat)
                        }
                )

                output$downloadLargeSubNetGenes <- shiny::downloadHandler(
                        filename = function(){
                                paste("top", input$l3_val, "_subNet_genes_", Sys.Date(), '.txt', sep='')
                        },
                        content = function(con){
                                write.table(V(myValues$largeSubPlot)$name, con, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
                        }
                )

		myValues$rankDF <- shiny::reactive({
			shiny::validate(
                                need(!is.null(myValues$enrichmentDF), "Waiting for enrichment!")
                        )

			progress <- shiny::Progress$new()
			progress$set(message="Getting Rank Table...", value=0)
			on.exit(progress$close())

			#orgDB <- unlist(strsplit(input$organism, ";"))[2]
			orgDB <- unlist(strsplit(input$organism, ";"))[1]
			orgDB <- get(orgDB)
			mapped_df <- select(orgDB, keys=myValues$rankedGenes, columns=c("ENSEMBL","ENTREZID","GENENAME"), keytype="SYMBOL")
			print(dim(mapped_df))
			print(head(mapped_df))
			#tmpDF <- cbind(gene_rank=sapply(mapped_df$SYMBOL, function(x) which(myINfORM$info_rank_list %in% x)), mapped_df)
			#tmpDF <- cbind(GENE_RANK=sapply(mapped_df$SYMBOL, function(x) which(myINfORM$info_rank_list %in% x)), mapped_df)
			tmpDF <- cbind(GENE_RANK=sapply(mapped_df$SYMBOL, function(x) which(myValues$rankedGenes %in% x)), mapped_df, SCORE=sapply(mapped_df$SYMBOL, function(x) V(myValues$iGraph)$score[which(V(myValues$iGraph)$name %in% x)]))
			#print(head(as.vector(sapply(mapped_df$SYMBOL, function(x) V(myINfORM$info_igraph)$corscore[which(V(myINfORM$info_igraph)$name %in% x)]))))
			mapped_df <- tmpDF[order(tmpDF$GENE_RANK),]
			print("Going to add member information...")

			print(mapped_df[which(mapped_df$SYMBOL %in% myValues$cList[[1]]),])
			mapped_df$CANDIDATE_GROUP <- "NA"
			mapped_df$CANDIDATE_GROUP[which(mapped_df$SYMBOL %in% myValues$cList[[3]])] <- input$l3_val
			mapped_df$CANDIDATE_GROUP[which(mapped_df$SYMBOL %in% myValues$cList[[2]])] <- input$l2_val
			mapped_df$CANDIDATE_GROUP[which(mapped_df$SYMBOL %in% myValues$cList[[1]])] <- input$l1_val

			mapped_df$MODULE_GROUP <- "NA"
			mapped_df$MODULE_GROUP[which(mapped_df$SYMBOL %in% myValues$mList[[3]])] <- input$l3_val
			mapped_df$MODULE_GROUP[which(mapped_df$SYMBOL %in% myValues$mList[[2]])] <- input$l2_val
			mapped_df$MODULE_GROUP[which(mapped_df$SYMBOL %in% myValues$mList[[1]])] <- input$l1_val

			mapped_df
		})
		output$rankDT <- DT::renderDataTable({
			shiny::validate(
                                need(!is.null(myValues$rankDF()), "Waiting for updated Gene List...")
                        )

			print("In Rank Data Frame, renderDataTable!")
			DT::datatable(myValues$rankDF(), filter = list(position='top', clear=FALSE), options = list(search = list(regex=TRUE, caseInsensitive=FALSE), scrollX=TRUE))
		})

		output$downloadRankDT <- shiny::downloadHandler(
			filename = function(){
				paste("Enriched_Gene_Table_", Sys.Date(), '.txt', sep='')
			},
			content = function(con){
				write.table(myValues$rankDF(), con, quote=FALSE, row.names=FALSE,  sep="\t")
			}
		)

		output$selOrganism <- shiny::renderUI({
			org_ann_libs <- rownames(installed.packages()[grep(".*\\.eg\\.db", rownames(installed.packages())),])
			lapply(org_ann_libs, require, character.only = TRUE)
			org_choices <- sapply(org_ann_libs, function(x){y=unname(packageDescription(x, fields="organism")); val=paste0(x,";",y); names(val) <- y; return(val)}, USE.NAMES=FALSE)
			shiny::selectInput("organism", "Choose Organism", choices=org_choices, multiple=FALSE, selected="org.Hs.eg.db;Homo sapiens")
		})

		output$selMethod <- shiny::renderUI({
			shiny::selectInput("method", "Select Inference Algorithm", choices=methods, multiple=TRUE, selected=methods)
		})
		
		output$selEst <- shiny::renderUI({
			shiny::selectInput("est", label="Select Correlation", choices=est.opt, multiple=TRUE, selected=est.opt)
		})
		
		output$selDisc <- shiny::renderUI({
			shiny::selectInput("disc", "Select Discretization Method", choices=disc.opt, multiple=TRUE, selected=disc.opt)
		})

		output$selCores <- shiny::renderUI({
			shiny::selectInput("cores", "Select Number of Cores", choices=(1:detectCores()), multiple=FALSE, selected=2)
		})

		output$rankAttr <- shiny::renderUI({
			#rank_attr <- c(net_attr[!net_attr %in% "eccentricity"], "corscore")
			shiny::selectInput("rankAttr", "Select Attributes for Ranking Candidates", choices=rank_attr, multiple=TRUE, selected=rank_attr)
		})
		
		shinyjs::onclick("toggleColScheme", shinyjs::toggle(id="colScheme", anim=TRUE))
		shinyjs::onclick("showAdvancedOptions", shinyjs::toggle(id="advancedOptions", anim=TRUE))
		shinyjs::onclick("showColorOptionsMain", shinyjs::toggle(id="colorOptionsMain", anim=TRUE))

                observeEvent(input$vColType, {
                        print("Observing Color Type Event!")
                        if(input$vColType=="score"){
                                shinyjs::hide(id="colByRank", anim=TRUE)
                                shinyjs::show(id="colByScore", anim=TRUE)
                        }else if(input$vColType=="rank"){
                                shinyjs::hide(id="colByScore", anim=TRUE)
                                shinyjs::show(id="colByRank", anim=TRUE)
                        }
                })

                observeEvent(input$bgCol, {
                        print("Observing Background Color Event!")
                        #cssStr <- paste0("background: ", input$bgCol, ";")
                        #inlineCSS(list(
                        #        "#mPlot" = cssStr
                        #))
                        js$plot_bg_col(input$bgCol)
                })

		shinyBS::addPopover(session, "IGraphVis_info", "Main Graph Info", content=paste0("<p>Graphical representation of the inferred network, nodes represent genes and edges represent the connection between the genes.</p> ", 
				"<p>Node color represents the positive and negative trend of gene by differential expression. ",
				"Edge color represents the positive or negative correlation between the genes.</p> ",
				"<p>Scroll up/down within the plot area to zoom in/zoom out, click and drag within the plot area to move the plot. ",
				"Hover over a node to see the node name, clicking a node opens NCBI ENTREZ GENE database record of the gene represented by that node.</p> ",
				"<p>Plot aesthetics and representation can be changed from the 'Show/Hide Aesthetics Options' section, click on the label to see the panel of options. ",
				"This panel provides option to change graph layout, specify percentile threshold to denote positively and negatively associated nodes. ",
				"Aesthetics can be altered to show gene rank as node color, increase/decrease node label size, alter the shape of the node and ",
				"specify edge width. Cutomization options for node color, edge color, background color, node border color are also provided.</p>"
			),
			placement="right",
			trigger="focus",
			#options=list(container="body", width="70%", "max-width"="70%")
			options=list(container="#main_info")
		)

		shinyBS::addPopover(session, "smallSubNetVis_info", "Small Set Graph Info", content=paste0("<p>Graphical representation of the sub-network extracted from the main network ", 
				"by using the top ranked genes specified in 'Small Set' option (default: 5) as seed genes. </p>",
				"<p>The size of top ranked genes is set from the 'Show/hide Advanced Options' panel. ",
				"Seed genes are highlighted as per the color scheme in 'Show/Hide Aesthetics Options' panel. ",
				"Aesthetics customization is same as main graph.</p>"
			),
			placement="right",
			trigger="focus",
			options=list(container="#small_info")
		)

		shinyBS::addPopover(session, "mediumSubNetVis_info", "Medium Set Graph Info", content=paste0("<p>Graphical representation of the sub-network extracted from the main network. ", 
				"by using the top ranked genes specified in 'Medium Set' option (default: 10) as seed genes. </p>",
				"<p>The size of top ranked genes is set from the 'Show/Hide Advanced Options' panel. ",
				"Seed genes are highlighted as per the color scheme in 'Show/hide Aesthetics Options' panel. ",
				"Aesthetics customization is same as main graph.</p>"
			),
			placement="right",
			trigger="focus",
			options=list(container="#medium_info")
		)

		shinyBS::addPopover(session, "largeSubNetVis_info", "Large Set Graph Info", content=paste0("<p>Graphical representation of the sub-network extracted from the main network. ", 
				"by using the top ranked genes specified in 'Large Set' option (default: 20) as seed genes. </p>",
				"<p>The size of top ranked genes is set from the 'Show/Hide Advanced Options' panel. ",
				"Seed genes are highlighted as per the color scheme in 'Show/hide Aesthetics Options' panel. ",
				"Aesthetics customization is same as main graph.</p>"
			),
			placement="right",
			trigger="focus",
			options=list(container="#large_info")
		)

		shinyBS::addPopover(session, "enrichmentDT_info", "Annotation Enrichment Info", content=paste0("<p>Annotation enrichment is performed by Fisher's Exact Test on the set of genes from each extracted sub-networks. </p>",
				"<p>GO term are filtered as enriched in small, medium and large gene set and the number of genes representing the GO term should be greater than or equal to the enrichment result from smaller gene set.</p>",
				"<p>Score associated with each GO term is the sum of normalized gene counts across the three enrichment results. </p>",
				"<p>List of genes from the enrichment result of the largest gene set is also reported. </p>"
			),
			placement="right",
			trigger="focus",
			options=list(container="#et_info")
		)

		shinyBS::addPopover(session, "tile_plot_panel_info", "Tile Plot Info", content=paste0("<p>The enriched GO terms are summarized based on their semantic similarity. </p>",
				"<p>This is represented as tile plot with GO terms clustered on the basis of their semantic similarity, major term from each cluster is chosen as the cluster representative. </p>",
				"<p>The area of tiles is representative of the associated score for each GO and aggregrated scores for the clusters. </p>",
				"<p>There is a tile plots for each GO type, Biological Process, Cellular Component and Molecular Function. </p>"
			),
			placement="right",
			trigger="focus",
			options=list(container="#tile_info")
		)

		shinyBS::addPopover(session, "rankDT_info", "Rank Table Info", content=paste0("<p>Table containing the final rank of the gene from the whole network inferred by INfORM, ",
				"HGNC Gene Symbol of the gene, unique identifier from the NCBI ENTREZ Gene database, descriptive name of the gene, user provided differential expression score, ",
				"membership of gene in top5, top10 and top20 candidate genes, membership of gene in response modules extracted by using the top5, top10 and top20 candidate genes. </p>",
				"<p>Membership information is cumulative, genes of compact membership also belong to the broader membership, eg: gene membership 5+10+20 = total genes in 20. </p>"
			),
			placement="right",
			trigger="focus",
			options=list(container="#rt_info")
		)
	}
)
