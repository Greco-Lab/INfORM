suppressMessages(library(shiny))
suppressMessages(library(shinyjs))
suppressMessages(library(shinyBS))
suppressMessages(library(shinydashboard))
suppressMessages(library(DT))
suppressMessages(library(R.utils))
suppressMessages(library(GSEABase))
suppressMessages(library(visNetwork))
suppressMessages(library(randomcoloR))
suppressMessages(library(radarchart))
#suppressMessages(library(Rserve))
suppressMessages(library(WriteXLS))
suppressMessages(library(igraph))
suppressMessages(library(plyr))
suppressMessages(library(GOSemSim))
print("Print Source Directory")
print(dirname(getSrcDirectory(function(x){x})))
functions_R <- file.path(dirname(getSrcDirectory(function(x){x})), "INfORM_functions.R")
source(functions_R)

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
		myValues <- shiny::reactiveValues(gxTable=NULL, dgxTable=NULL, gxCorMat=NULL, iGraph=NULL, updateTabCorMat=0, cList=NULL, mList=NULL, GO_clust_summ_list=NULL)

                ##Close R session on closing the graphical window
                #session$onSessionEnded(function(){
                #        stopApp()
                #        q("no")
                #})

		myValues$sepChoices <- c("TAB", ",", ";", "SPACE", "OTHER")
		myValues$quoteChoices <- c(NA, "SINGLE", "DOUBLE")
                myValues$envir <- environment()

		myValues$inputGx <- observeEvent(input$upload_gx_submit, {
                        gxFile <- input$gx
			##if (is.null(gxFile) || is.null(myValues$dgxTable()))
			#if (is.null(gxFile) || is.null(myValues$dgxTable))
			if (is.null(gxFile))
			return(NULL)

                        sepS <- input$gxSepS
			sepT <- input$gxSepT
			sepChar=NULL
			if(sepS=="OTHER"){
				sepChar <- sepT
			}else{
				if(sepS=="TAB"){
					sepChar="\t"
				}else if(sepS=="SPACE"){
					sepChar=" "
				}else{
					sepChar=sepS
				}
			}
			print(sepChar)

			quote <- input$gxQuote
			print(quote)
			print(is.na(quote))
			if(is.na(quote) || quote=="NA"){
				quote <- ""
			}else if(quote=="SINGLE"){
				quote <- "'"
			}else if(quote=="DOUBLE"){
				quote <- '"'
			}

			print("Reading loaded gx File...")
			#gx <- read.csv(gxFile$datapath, row.names=1, header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE, quote="")
			gx <- read.csv(gxFile$datapath, row.names=1, header=TRUE, sep=sepChar, stringsAsFactors=FALSE, quote=quote, as.is=TRUE, strip.white=TRUE, check.names=FALSE)
                        myValues$inputGx <- gx

			updateTextInput(session, "gxUploadDisp", value=gxFile$name)
			shinyBS::toggleModal(session, "importGxModal", toggle="close")
                        shinyBS::updateButton(session, "import_gx_submit", style="success", icon=icon("check-circle"))
		})

		myValues$gxTable <- shiny::reactive({
			#gxFile <- input$gx
			##if (is.null(gxFile) || is.null(myValues$dgxTable()))
			#if (is.null(gxFile) || is.null(myValues$dgxTable))
			#return(NULL)

			if (is.null(myValues$inputGx) || is.null(myValues$dgxTable))
			return(NULL)

			gx <- myValues$inputGx
                        gx <- gx[rownames(myValues$dgxTable),]
                        return(gx)
		})

		#myValues$dgxTable <- shiny::reactive({
		#	print("Checking dgx File...")
		#	dgxFile <- input$dgx
		#	if (is.null(dgxFile))
		#	return(NULL)

		#	print("Reading loaded dgx File...")
		#	dgx <- read.csv(dgxFile$datapath, row.names=1, header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE, quote="")
                #        dgx <- dgx[order(rownames(dgx)),,drop=F]
                #        print(str(dgx))
                #        dgx
		#})

		myValues$inputDgx <- eventReactive(input$load_dgx_submit, {
			if(is.null(input$dgx))
			return(NULL)

			dgxFile <- input$dgx
			sepS <- input$sepS
			sepT <- input$sepT
			sepChar=NULL
			if(sepS=="OTHER"){
				sepChar <- sepT
			}else{
				if(sepS=="TAB"){
					sepChar="\t"
				}else if(sepS=="SPACE"){
					sepChar=" "
				}else{
					sepChar=sepS
				}
			}
			print(sepChar)

			quote <- input$quote
			print(quote)
			print(is.na(quote))
			if(is.na(quote) || quote=="NA"){
				quote <- ""
			}else if(quote=="SINGLE"){
				quote <- "'"
			}else if(quote=="DOUBLE"){
				quote <- '"'
			}

			rowNames <- NULL
			con <- file(dgxFile$datapath, "r", blocking=FALSE)
			fileLines <- readLines(con)
			fileLines <- gsub("\t\t", "\tNA\t", fileLines)
			fileLines <- gsub("\t$", "\tNA", fileLines)
			close(con)
			colNumByRowDist <- table(sapply(fileLines, function(x) {length(strsplit(x, sepChar)[[1]])}, USE.NAMES=FALSE))
			if(length(colNumByRowDist) > 1){
				fileLines2 <- fileLines[-1]
				colNumByRowDist2 <- table(sapply(fileLines2, function(x) {length(strsplit(x, sepChar)[[1]])}, USE.NAMES=FALSE))
				if(length(colNumByRowDist2) > 1){
					shinyjs::info(paste0("Separating character '", sepChar, "', results in inconsistent number of columns!\n\nPlease check the input file format and select the correct separating character!"))
					myValues$dgxLoaded <- NULL
					return(NULL)
				}
				rowNames <- 1
			}

			myValues$dgxLoaded <- 1
			myValues$dgxFileName <- dgxFile$name
			dgxTable <- read.csv(dgxFile$datapath, row.names=1, header=TRUE, sep=sepChar, stringsAsFactors=FALSE, quote=quote, as.is=TRUE, strip.white=TRUE, check.names=FALSE)

                        print("str(dgxTable) -- after:")
                        print(str(dgxTable))
			return(dgxTable)
		})

		myValues$dgxColChoices <- reactive({
			if(is.null(myValues$dgxLoaded))
			return(c("NA"))

			choicesVec <- seq(1,ncol(myValues$inputDgx()))
			choicesNames <- paste0("Column ", choicesVec)
			names(choicesVec) <- choicesNames
			return(choicesVec)
		})

		output$dgxDT <- DT::renderDataTable({
			shiny::validate(
				need(!is.null(myValues$inputDgx()), "No Differential Expression file!")
			)

                        dgxTable <- myValues$inputDgx()
                        colnames(dgxTable) <- paste0(colnames(dgxTable), " [", c(1:ncol(dgxTable)), "]")
			DT::datatable(dgxTable, filter="none", 
				options = list(
					ch = list(regex=TRUE, caseInsensitive=FALSE), 
					scrollX=TRUE, 
					pageLength=2,
					lengthMenu=c(1,2,3),
					ordering=FALSE
				)
			)
		},server=TRUE)

		observeEvent(input$upload_dgx_submit, {
			shiny::validate(
				need(!is.null(myValues$inputDgx()), "No Differential Expression File Provided!")
			)
			
                        dgxTable <- myValues$inputDgx()
			pvColID <- as.integer(input$pvCol)
			lfcColID <- as.integer(input$lfcCol)

                        pvals <- dgxTable[,pvColID]
                        lfcs <- dgxTable[,lfcColID]
                        if(any(pvals=="") || any(pvals==" ") || any(is.na(pvals))){
                                shinyjs::info(paste0("The specified P.Value column contains BLANK and/or NA values!\n\nPlease check the differential expression data and ensure that the selected P.Value column has complete information!"))
                                return(NULL)
                        }else if(any(pvals<0) || any(pvals>1)){
                                shinyjs::info(paste0("The specified P.Value column contains values out of 0-1 range!\n\nPlease check the differential expression data and ensure that the selected P.Value column has complete information!"))
                                return(NULL)
                        }else if(any(lfcs=="") || any(lfcs==" ") || any(is.na(lfcs))){
                                shinyjs::info(paste0("The specified LogFC column contains BLANK and/or NA values!\n\nPlease check the differential expression data and ensure that the selected LogFC column has complete information!"))
                                return(NULL)
                        }

			colIdx <- which(colnames(dgxTable) %in% c(pvColID,lfcColID))
			if(length(colIdx)>0){
				dgxTable <- dgxTable[,colIdx]
			}
                        myValues$dgxTable <- dgxTable
			myValues$pvColID <- pvColID
			myValues$lfcColID <- lfcColID

			updateTextInput(session, "dgxUploadDisp", value=myValues$dgxFileName)
			shinyBS::toggleModal(session, "importDgxModal", toggle="close")
                        shinyBS::updateButton(session, "import_dgx_submit", style="success", icon=icon("check-circle"))
                })

                myValues$method <- shiny::reactive({
                        if(is.null(input$method))
                        #return(methods)
                        return(NULL)

                        input$method
                })

                myValues$est <- shiny::reactive({
                        if(is.null(input$est))
                        #return(est.opt)
                        return(NULL)

                        input$est
                })

                myValues$disc <- shiny::reactive({
                        if(is.null(input$disc))
                        #return(disc.opt)
                        return(NULL)

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
                       
                        if(input$ensembleStrat!="user"){
                                if(is.null(myValues$method())){
                                        shinyjs::info(paste0("No Inference Algorithm Selected!"))
                                }

                                if(is.null(myValues$est())){
                                        shinyjs::info(paste0("No Correlation Selected!"))
                                }

                                if(is.null(myValues$disc())){
                                        shinyjs::info(paste0("No Discretization Method Selected!"))
                                }

                                if(any(
                                        is.null(myValues$gxTable()), 
                                        is.null(myValues$dgxLoaded), 
                                        is.null(myValues$method()), 
                                        is.null(myValues$est()), 
                                        is.null(myValues$disc())
                                ))
                                return(NULL)
                        }

			localDgxTable <- myValues$dgxTable
			rOrgDB <- unlist(strsplit(input$organism, ";"))[1]

                        ##Check "SYMBOL" mapping to annotation library
                        annDB <- get(rOrgDB)
                        allSymbols <- keys(annDB, keytype="SYMBOL")

                        #Check differential expression file
                        rowNamesDgx <- rownames(localDgxTable)
                        missingSymbolsIdx <- which(!rowNamesDgx %in% allSymbols)
                        if(length(missingSymbolsIdx)>0){
                                shinyjs::info(paste0("Unable to map ", length(missingSymbolsIdx), " of ", length(rowNamesDgx), 
                                "Gene Symbols from Differential Expression Table to the selected annotation library ", rOrgDB, "!\n\n",
                                "Missing IDs: \n",
                                paste0(rowNamesDgx[missingSymbolsIdx], collapse=", ")
                                ))
                                if((length(missingSymbolsIdx)/length(rowNamesDgx))*100>20){
                                        shinyjs::info("Missing too many genes from Differential Expression Table.\nEXITING!!!")
                                        return(NULL)
                                }
                                #return(NULL)
                                localDgxTable <- localDgxTable[-missingSymbolsIdx,]
                                rowNamesDgx <- rownames(localDgxTable)
                        }

                        localGxTable <- myValues$gxTable()
                        localGxTable <- localGxTable[rowNamesDgx,]

                        #Check differential expression vs expression matrix
                        rowNamesDgx <- rownames(localDgxTable)
                        rowNamesGx <- rownames(localGxTable)
                        missingSymbolsIdx <- which(!rowNamesDgx %in% rowNamesGx)
                        if(length(missingSymbolsIdx)>0){
                                shinyjs::info(paste0("Unable to map ", length(missingSymbolsIdx), " of ", length(rowNamesDgx), 
                                "Gene Symbols from Differential Expression Table to the Gene Expression Matrix!\n\n",
                                "Missing IDs: \n",
                                paste0(rowNamesDgx[missingSymbolsIdx], collapse=", ")
                                ))

                                if((length(missingSymbolsIdx)/length(rowNamesDgx))*100>20){
                                        shinyjs::info("Missing too many genes from Differential Expression Table.\nEXITING!!!")
                                        return(NULL)
                                }
                                localDgxTable <- localDgxTable[-missingSymbolsIdx,]
                                rowNamesDgx <- rownames(localDgxTable)
                                localGxTable <- localGxTable[rowNamesDgx,]
                        }

                        #Check gene expression matrix
                        rowNamesGx <- rownames(localGxTable)
                        missingSymbolsIdx <- which(!rowNamesGx %in% allSymbols)
                        if(length(missingSymbolsIdx)>0){
                                shinyjs::info(paste0("Unable to map ", length(missingSymbolsIdx), " of ", length(rowNamesGx), 
                                "Gene Symbols from Gene Expression Table to the selected annotation library ", rOrgDB, "!\n\n",
                                "Missing IDs: \n",
                                paste0(rowNamesDgx[missingSymbolsIdx], collapse=", ")
                                ))
                                localGxTable <- localGxTable[-missingSymbolsIdx,]
                                rowNamesGx <- rownames(localGxTable)
                                localDgxTable <- localDgxTable[rowNamesGx,]
                                #return(NULL)
                        }

                        runMINET <- FALSE
                        if(length(grep("minet", input$ensembleStrat))>0){
                                runMINET <- TRUE
                        }
                        
                        if(length(grep("user", input$ensembleStrat))>0){
                                matDF <- input$mat
                                userMatList <- list()
                                fileNames <- matDF$name
                                filePaths <- matDF$datapath
                                for(i in c(1:length(filePaths))){
                                        userMat <- as.matrix(read.csv(filePaths[i], row.names=1, header=TRUE, sep="\t", check.names=FALSE))
                                        if(any(is.na(userMat))){
                                                shinyjs::info(paste0("Uploaded matrix contains NA values!\n\nPlease ensure that the input matrix contains SCORES or RANKS as weights!"))
                                                return(NULL)
                                        }
                                        if(all(userMat %in% c(0,1))){
                                                shinyjs::info(paste0("Uploaded matrix is a binary matrix with values 0 and 1!\n\nPlease ensure that the input matrix contains SCORES or RANKS as weights!"))
                                                return(NULL)
                                        }

                                        #Check integrity of user uploaded adjacency matrix
                                        if(nrow(userMat)!=ncol(userMat)){
                                                shinyjs::info(paste0("File:", fileNames[i], " has differring number of rows and columns, ROWS:", nrow(userMat), " & COLS:", ncol(userMat), "\n\n Please check for consistency and rerun analysis! EXITING!!!" 
                                                ))
                                                return(NULL)
                                        }

                                        rowNamesMat <- rownames(userMat)
                                        colNamesMat <- colnames(userMat)
                                        errorIdx <- which(!rowNamesMat %in% colNamesMat)
                                        errorIdx2 <- which(!colNamesMat %in% rowNamesMat)
                                        if(any(length(errorIdx)>0, length(errorIdx2)>0)){
                                                shinyjs::info(paste0("File:", fileNames[i], " has mismatches betweem column names and row names.\n\nPlease check for consistency and rerun analysis! EXITING!!!" 
                                                ))
                                                return(NULL)
                                        }

                                        userMat <- userMat[,rowNamesMat]
                                        #Check differential expression vs expression matrix
                                        rowNamesDgx <- rownames(localDgxTable)
                                        missingSymbolsIdx <- which(!rowNamesDgx %in% rowNamesMat)
                                        if(length(missingSymbolsIdx)>0){
                                                shinyjs::info(paste0("Unable to map ", length(missingSymbolsIdx), " of ", length(rowNamesDgx), 
                                                "Gene Symbols from Differential Expression Table to File:", fileNames[i], "!\n\n",
                                                "Missing IDs: \n",
                                                paste0(rowNamesDgx[missingSymbolsIdx], collapse=", ")
                                                ))

                                                if((length(missingSymbolsIdx)/length(rowNamesDgx))*100>20){
                                                        shinyjs::info("Missing too many genes from Differential Expression Table.\nEXITING!!!")
                                                        return(NULL)
                                                }
                                                localDgxTable <- localDgxTable[-missingSymbolsIdx,]
                                                rowNamesDgx <- rownames(localDgxTable)
                                                localGxTable <- localGxTable[rowNamesDgx,]
                                                #userMat <- userMat[rowNamesDgx,rowNamesDgx]
                                        }
                                        userMatList[[i]] <- userMat
                                }

                                if(length(userMatList)>0){
                                        rowNamesDgx <- rownames(localDgxTable)
                                        for(i in c(1:length(userMatList))){
                                                userMatList[[i]] <- userMatList[[i]][rowNamesDgx,rowNamesDgx]
                                        }
                                }
                                #return(NULL)
                        }

                        if(input$ensembleStrat=="minet+user"){
                                rowNamesDgx <- rownames(localDgxTable)
                                localGxTable <- localGxTable[rowNamesDgx,]
                        }

			print("Passed Validations.....")

                        progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}
                        on.exit({
				progress$close()
				shiny::updateTabsetPanel(session, "display", selected="corMat_display")
			})

                        #print("localGxTable : ")
                        #print(dim(localGxTable))
                        #print(str(localGxTable))
                        #print("localDgxTable : ")
                        #print(dim(localDgxTable))
                        #print(str(localDgxTable))
                        #return(NULL)

                        corMat <- NULL
                        #matFile <- input$mat
			if(runMINET || length(userMatList)>1){
                                myMethod <- myValues$method()
                                myEst <- myValues$est()
                                myDisc <- myValues$disc()
                                myCores <- as.numeric(myValues$cores())

                                #print(paste0("methods : ", myMethod))
                                #print(paste0("ESTs : ", myEst))
                                #print(paste0("Discs : ", myDisc))

                                progress$set(message="Inference", value=0)

                                rankMat <- get_ranked_consensus_matrix(gx_table=localGxTable, iMethods=myMethod, iEst=myEst, iDisc=myDisc, summ_by=input$summType, score_type=input$bordaScore, ncores=myCores, matList=userMatList, mat_weights=input$matWeights, ensemble_strategy=input$ensembleStrat, debug_output=FALSE, updateProgress=updateProgress)
                                resL <- parse_edge_rank_matrix(rankMat, mat_weights=input$matWeights, edge_selection_strategy=input$selEdge, topN=input$topCutOff)
                                corMat <- resL[["bin_mat"]]
                                edgeRank <- resL[["edge_rank"]]
                        }else{
                                rankMat <- userMatList[[1]]
                                resL <- parse_edge_rank_matrix(rankMat, mat_weights=input$matWeights, edge_selection_strategy=input$selEdge, topN=input$topCutOff)
                                corMat <- resL[["bin_mat"]]
                                edgeRank <- resL[["edge_rank"]]

                                print("Bin Mat Str : ")
                                print(str(corMat))
                                print("Edge Rank : ")
                                print(length(edgeRank))
                                print(max(edgeRank))
                                myValues$updateTabCorMat <- 1
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
                        print("Counts localIGraph : ")
                        print(vcount(localIGraph))
                        print(ecount(localIGraph))

                        #Setting edge attribute rank
                        edgeIDs <- igraph::ends(localIGraph, igraph::E(localIGraph), names=TRUE)
			igraph::edge_attr(localIGraph, name="rank") <- apply(edgeIDs, 1, function(x) {
				rankMat[x[1],x[2]]
			})

                        print("Setting DGX values as 'score'")
                        ###vertex_attr(localIGraph, name="score") <- localDgxTable[V(localIGraph)$name,]
			##pval <- localDgxTable[V(localIGraph)$name,"P.Value"]
			##lfc <- localDgxTable[V(localIGraph)$name,"logFC"]
			#pval <- localDgxTable[V(localIGraph)$name, 2]
			#lfc <- localDgxTable[V(localIGraph)$name, 1]
			pvColID <- myValues$pvColID
			lfcColID <- myValues$lfcColID
			pval <- localDgxTable[V(localIGraph)$name, pvColID]
			lfc <- localDgxTable[V(localIGraph)$name, lfcColID]
			score <- lfc * -log10(pval)
			print("COL IDs:")
			print("pvColID:")
			print(pvColID)
			print("lfcColID:")
			print(lfcColID)
			print("DGX TABLE INFO:")
			print(str(localDgxTable))
			print(head(localDgxTable))
			print("EXTRACTED PVAL INFO:")
			print(head(pval))
			print("EXTRACTED LFC INFO:")
			print(head(lfc))
                        igraph::vertex_attr(localIGraph, name="score") <- score
			igraph::vertex_attr(localIGraph, name="pval") <- pval
			igraph::vertex_attr(localIGraph, name="lfc") <- lfc

                        progress$set(message="Ranking Genes", value=0)
                        rankAttr <- myValues$rankAttr()
			rankAttr.c <- rank_attr[-which(rank_attr=="score")]
                        localRankedGenes <- get_ranked_gene_list(localIGraph, rank_list_attr=rankAttr, debug_output=FALSE)
                        localRankedGenes.c <- get_ranked_gene_list(localIGraph, rank_list_attr=rankAttr.c, debug_output=FALSE)
                        localRankedGenes.pv <- V(localIGraph)$name[order(V(localIGraph)$pval, decreasing=FALSE)]
                        localRankedGenes.lfc <- V(localIGraph)$name[order(abs(V(localIGraph)$lfc), decreasing=TRUE)]

                        myValues$rankedGenes <- localRankedGenes

			progress$set(message="Module Detection", value=0)
			updateProgress(detail="Getting Modules...", value=1/4)

                        modMethod <- input$modMethod
			minModSize <- input$minModSize

                        modules <- get_modules(iGraph=localIGraph, method=modMethod)
			modulesCount <- sum(sizes(modules)>=minModSize)
			print("modulesCount")
			print(modulesCount)
			if(modulesCount==0){
                                        shinyjs::info(paste0("No modules detected by ", modMethod, "!\n\n",
					"Please adjust MINET setup and module detection parameters and rerun!"
					))
                                        return(NULL)
			}
			modPalette <- setNames(randomcoloR::distinctColorPalette(modulesCount), paste0(modMethod, "_", c(1:modulesCount)))
			members <- membership(modules)

			#modMembers <- NULL
			idx <- which(table(members)<minModSize)
			if(length(idx)>0){
				orphanIdx <- names(idx)
				j <- 1
				for(i in 1:c(length(modules))){
					#i <- as.character(i)
					if(i %in% orphanIdx){
						#tmpVecNames <- names(which(members %in% i))
						#tmpVec <- setNames("orphan", tmpVecNames)
						#modMembers <- c(modMembers, tmpVec)
						members[which(members %in% i)] <- "orphan"
					}else{
						#tmpVecNames <- names(which(members %in% i))
						#tmpVec <- setNames(paste0(modMethod, "_", j), tmpVecNames)
						#modMembers <- c(modMembers, tmpVec)
						members[which(members %in% i)] <- paste0(modMethod, "_", j)
						j <- j+1
					}
				}
			}else{
				members <- setNames(paste0(modMethod, "_", members), names(members))
			}
			print("modMems:")
			print(members)
			igraph::vertex_attr(localIGraph, name="group") <- members[V(localIGraph)$name]

			updateProgress(detail="Annotating Modules...", value=2/4)
			modules_ll <- annotate_modules(iGraph=localIGraph, modules=modules, rl=localRankedGenes, rl.c=localRankedGenes.c, rl.pv=localRankedGenes.pv, rl.lfc=localRankedGenes.lfc, rl.edge=edgeRank, annDB=rOrgDB, p_value=input$sigPval, min_mod_size=minModSize, prefix=modMethod)

			updateProgress(detail="Generating Module Score Table...", value=3/4)
			#Create IC list for semsim
			IC_ll <- sapply(c("BP", "CC", "MF"), function(typ) godata(rOrgDB, ont=typ, computeIC=TRUE))
			#modules_rank_table <- as.rank.table(iGraph=localIGraph, annotated_modules_ll=modules_ll, IC_ll=IC_ll)
			modules_rank_table <- as.rank.table(iGraph=localIGraph, annotated_modules_ll=modules_ll) ##Not computing GO sim

			updateProgress(detail="Computing Module Similarity by GO...", value=4/4)
			within_go_sim <- get_between_GO_sim(modules_ll, modules_ll, IC_ll=IC_ll)
			#simHeatMap <- plot_sim_heatmap(hSim=within_go_sim, rSim=within_go_sim, cSim=within_go_sim)

			updateProgress(detail="Finishing: ", value=3/3)

                        myValues$iGraph <- localIGraph
			myValues$gxCorMat <- localGxCorMat
			myValues$rankMat <- rankMat
			myValues$rl <- localRankedGenes
			myValues$rl.c <- localRankedGenes.c
			myValues$rl.pv <- localRankedGenes.pv
			myValues$rl.lfc <- localRankedGenes.lfc
			myValues$edgeRank <- edgeRank
			myValues$modules <- modules
			myValues$modules_ll <- modules_ll
			myValues$modules_rank_table <- modules_rank_table
			myValues$modPalette <- modPalette
			myValues$IC_ll <- IC_ll
			#myValues$simHeatMap <- simHeatMap
			myValues$within_go_sim <- within_go_sim
                        myValues$modMethodVal <- modMethod
			myValues$minModSizeVal <- minModSize

                        #Response chart variables
                        myValues$GO_clust_summ_list <- NULL
                        myValues$optAnnotation <- NULL
                        myValues$optAE <- NULL
                        myValues$morphed_gene_ll <- NULL
                        myValues$morphed_ae_ll <- NULL
                        myValues$optModVal <- NULL
                        myValues$conModVal <- NULL
                        myValues$indModVal <- NULL
                })

		observeEvent(input$runDetect, {
			print("Run Module Detection.....")

                        if(any(
                                is.null(myValues$gxTable()), 
                                is.null(myValues$dgxLoaded), 
                                is.null(myValues$method()), 
                                is.null(myValues$est()), 
                                is.null(myValues$disc())
                        ))
                        return(NULL)

			print("Passed Validations.....")

                        progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}
                        on.exit({
				progress$close()
				shiny::updateTabsetPanel(session, "display", selected="corMat_display")
			})

			progress$set(message="Module Detection", value=0)
			updateProgress(detail="Getting Modules...", value=1/4)

                        localIGraph <- myValues$iGraph
			localRankedGenes <- myValues$rl
			localRankedGenes.c <- myValues$rl.c
			localRankedGenes.pv <- myValues$rl.pv
			localRankedGenes.lfc <- myValues$rl.lfc
			edgeRank <- myValues$edgeRank
			rOrgDB <- unlist(strsplit(input$organism, ";"))[1]
                        modMethod <- input$modMethod
			minModSize <- input$minModSize

                        modules <- get_modules(iGraph=localIGraph, method=modMethod)
			modulesCount <- sum(sizes(modules)>=minModSize)
			modPalette <- setNames(randomcoloR::distinctColorPalette(modulesCount), paste0(modMethod, "_", c(1:modulesCount)))
			members <- membership(modules)

			idx <- which(table(members)<minModSize)
			if(length(idx)>0){
				orphanIdx <- names(idx)
				j <- 1
				for(i in 1:c(length(modules))){
					if(i %in% orphanIdx){
						members[which(members %in% i)] <- "orphan"
					}else{
						members[which(members %in% i)] <- paste0(modMethod, "_", j)
						j <- j+1
					}
				}
			}else{
				members <- setNames(paste0(modMethod, "_", members), names(members))
			}
			print("modMems:")
			print(members)
			igraph::vertex_attr(localIGraph, name="group") <- members[V(localIGraph)$name]

			updateProgress(detail="Annotating Modules...", value=2/4)
			modules_ll <- annotate_modules(iGraph=localIGraph, modules=modules, rl=localRankedGenes, rl.c=localRankedGenes.c, rl.pv=localRankedGenes.pv, rl.lfc=localRankedGenes.lfc, rl.edge=edgeRank, annDB=rOrgDB, p_value=input$sigPval, min_mod_size=minModSize, prefix=modMethod)

			updateProgress(detail="Generating Module Score Table...", value=3/4)
			#Create IC list for semsim
			IC_ll <- sapply(c("BP", "CC", "MF"), function(typ) godata(rOrgDB, ont=typ, computeIC=TRUE))
			#modules_rank_table <- as.rank.table(iGraph=localIGraph, annotated_modules_ll=modules_ll, IC_ll=IC_ll)
			modules_rank_table <- as.rank.table(iGraph=localIGraph, annotated_modules_ll=modules_ll) ##Not computing GO sim

			updateProgress(detail="Computing Module Similarity by GO...", value=4/4)
			within_go_sim <- get_between_GO_sim(modules_ll, modules_ll, IC_ll=IC_ll)
			#simHeatMap <- plot_sim_heatmap(hSim=within_go_sim, rSim=within_go_sim, cSim=within_go_sim)

			updateProgress(detail="Finishing: ", value=3/3)

                        print("Updating Select Input Varaibles#########")
                        tmpChoices <- names(modules_ll[["ig"]])
                        print(tmpChoices)
                        updateSelectInput(session, "optMod", choices=tmpChoices, selected="")
                        updateSelectInput(session, "conMod", choices=tmpChoices, selected="")
                        updateSelectInput(session, "indMod", choices=tmpChoices, selected="")
                        modsOptCon <- unique(c(input$optMod, input$conMod))
                        modsOptInd <- unique(c(input$optMod, input$indMod))
                        modsConInd <- unique(c(input$conMod, input$indMod))
                        print(modsOptCon)
                        print(modsOptInd)
                        print(modsConInd)
			myValues$modules <- modules
			myValues$modules_ll <- modules_ll
			myValues$modules_rank_table <- modules_rank_table
			myValues$modPalette <- modPalette
			myValues$IC_ll <- IC_ll
			myValues$within_go_sim <- within_go_sim
                        myValues$modMethodVal <- modMethod
			myValues$minModSizeVal <- minModSize

                        #Response chart variables
                        myValues$GO_clust_summ_list <- NULL
                        myValues$optAnnotation <- NULL
                        myValues$optAE <- NULL
                        myValues$morphed_gene_ll <- NULL
                        myValues$morphed_ae_ll <- NULL
                        myValues$optModVal <- NULL
                        myValues$conModVal <- NULL
                        myValues$indModVal <- NULL
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
				paste0(geneCount), "Total Genes", icon=icon("list"), color="light-blue", width=3
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
				paste0(connectedGeneCount), "Connected Genes", icon=icon("thumbs-o-up"), color="olive", width=3
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
				paste0(unConnectedGeneCount), "Un-connected Genes", icon=icon("thumbs-o-down"), color="red", width=3
			)
		})

		output$netDensityBox <- shinydashboard::renderValueBox({
			if(is.null(myValues$iGraph)){
				netDensity <- "NA"
			}
			else{
			    	netDensity <- paste0(round(igraph::graph.density(myValues$iGraph)*100, 2), "%")
			}
                        
			shinydashboard::valueBox(
				netDensity, "Network Density", icon=icon("area-chart"), color="purple", width=3
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
                                #paste("Adjacency_Matrix_", Sys.Date(), '.txt', sep='')
                                paste("Edge_Rank_Matrix_", Sys.Date(), '.txt', sep='')
                        },
			content = function(con){
				#write.table(myValues$gxCorMat, con, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
				write.table(myValues$rankMat, con, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
			}
		)

		shiny::observe({
			if(input$ensembleStrat=="user"){
				shinyjs::disable("calcMat")
				shinyjs::disable("method")
				shinyjs::disable("est")
				shinyjs::disable("disc")
				shinyjs::disable("cores")
				#shinyjs::disable("selEdge")
			}
			else{
				shinyjs::enable("calcMat")
				shinyjs::enable("method")
				shinyjs::enable("est")
				shinyjs::enable("disc")
				shinyjs::enable("cores")
				#shinyjs::enable("selEdge")
			}

                        if(input$ensembleStrat=="minet"){
                                if(is.null(input$gx) || is.null(input$dgx)){
                                        shinyjs::disable("runINfORM")
                                }else{
                                        shinyjs::enable("runINfORM")
                                }
                        }else if(input$ensembleStrat=="user"){
                                if(is.null(input$dgx) || is.null(input$mat)){
                                        shinyjs::disable("runINfORM")
                                }else{
                                        shinyjs::enable("runINfORM")
                                }
                        }else if(input$ensembleStrat=="minet+user"){
                                if(is.null(input$gx) || is.null(input$dgx) || is.null(input$mat)){
                                        shinyjs::disable("runINfORM")
                                }else{
                                        shinyjs::enable("runINfORM")
                                }
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

			if(!is.null(input$sepS) && input$sepS=="OTHER"){
				shinyjs::enable("sepT")
			}else{
				shinyjs::disable("sepT")
			}

                        if(!is.null(input$gxSepS) && input$gxSepS=="OTHER"){
				shinyjs::enable("gxSepT")
			}else{
				shinyjs::disable("gxSepT")
			}

			if(input$selEdge=="default"){
				shinyjs::disable("topCutOff")
                                shinyjs::hide(id="topCutOff", anim=TRUE)
			}else{
				shinyjs::enable("topCutOff")
                                shinyjs::show(id="topCutOff", anim=TRUE)
			}

                        if(is.null(myValues$modules)){
				shinyjs::disable("runDetect")
                        }else{
				shinyjs::enable("runDetect")
                        }

                        if(is.null(input$optMod) || length(input$optMod)<2 || input$optMod==""){
				shinyjs::disable("rChartButton")
				shinyBS::addTooltip(session, id="rChartButton", title="Please provide at least two module to merge!", placement="bottom")
                        }else{
				shinyjs::enable("rChartButton")
				shinyBS::removeTooltip(session, id="rChartButton")
                        }

                        if(input$ensembleStrat=="minet"){
				shinyjs::hide(id="matOptions", anim=TRUE)
                        }else{
				shinyjs::show(id="matOptions", anim=TRUE)
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
			#localDgxTable <- myValues$dgxTable()
			localDgxTable <- myValues$dgxTable
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

                        if(input$chkGrpCol==TRUE){
                                modPalette <- myValues$modPalette
                                modVec <- igraph::vertex_attr(localIGraph, name="group")
				print("modVec")
				print(modVec)
                                mods <- unique(modVec)
                                for(mod in mods){
					if(mod == "orphan")
					next

                                        idx <- which(modVec %in% mod)
					if(length(idx)>0){
						igraph::vertex_attr(localIGraph, name="color")[idx] <- modPalette[mod]
					}
                                }
                        }else if(input$vColType=="score"){
                        #if(input$vColType=="score"){
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
                                need(!is.null(myValues$iGraph_dynamic()), "Waiting for updated graph...")
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

                        if(input$chkGrpCol==TRUE){
                                modPalette <- myValues$modPalette
                                print("str(modPalette):")
                                print(str(modPalette))
                                #modVec <- igraph::vertex_attr(localIGraph, name="group")
                                #mods <- unique(modVec)
                                addNodes <- data.frame(label=names(modPalette), shape="circle", icon.color=modPalette)
                                print("addNodes:")
                                print(addNodes)
                                mainVisNetwork <- mainVisNetwork %>%
                                visLegend(addNodes=addNodes, useGroups=FALSE, main="Modules")
                        }

			#Add export control
			networkImageExportFileName <- paste0("network_", Sys.Date())
			imageType <- "png"
                        bgColor <- input$bgCol
			mainVisNetwork <- mainVisNetwork %>%
			visExport(type=imageType, name=networkImageExportFileName, label=paste0("Export as ", imageType), background=bgColor, float="right", style=NULL, loadDependencies=TRUE)
			
                        mainVisNetwork
                })

                #observeEvent(input$chkGrpCol, {
		#	#Color Groups 
		#	chkGrpCol <- input$chkGrpCol
		#	isolate({
		#		if(chkGrpCol==TRUE){
		#			modPalette <- myValues$modPalette
                #                        modVec <- igraph::vertex_attr(localIGraph, name="group")
                #                        mods <- unique(modVec)
                #                        addNodes <- data.frame(label=mods, shape="circle", icon.color=modPalette[mods])
                #                        for(mod in mods){
                #                                visNetworkProxy("IGraphVis") %>%
                #                                visGroups(groupname=mod, color=modPalette[mod])
                #                        }
		#		}
                #                visNetworkProxy("IGraphVis") %>%
                #                visLegend()
		#	})
		#})

		observeEvent(input$focusMod, {
			#Focus on 
			gr <- input$focusMod
			isolate({
				if(gr != "ALL"){
					ig <- myValues$iGraph_dynamic()
					nodeDF <- as_data_frame(ig, what="vertices")
					names <- nodeDF$name[nodeDF$group %in% gr]
				}else{
					names <- NULL
				}
				visNetworkProxy("IGraphVis") %>%
				visFit(nodes=names, animation=list(duration=200, easingFunction="easeInOutQuad"))
			})
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
				paste("enriched_GO_", Sys.Date(), '.xlsx', sep='')
			},
			content = function(con){
                                modules_ll <- myValues$modules_ll
                                ae_list <- modules_ll[["ae"]]
                                idx <- which(unlist(lapply(ae_list, is.null), use.names=FALSE))
                                if(length(idx)>0){
                                            ae_list <- ae_list[-idx]
                                }
                                #names(ae_list) <- paste0("mod_", names(ae_list))
                                WriteXLS(ae_list, ExcelFileName=con, col.names=TRUE, AdjWidth=TRUE, BoldHeaderRow=TRUE)
			}
		)

		output$enrichmentDT <- DT::renderDataTable({
			shiny::validate(
                                need(!is.null(myValues$modules_ll), "Waiting for module detection...")
                        )
                        shiny::validate(
                                need(!is.null(input$modEnrich), "Waiting for module detection...")
                        )
                        shiny::validate(
                                need(!is.na(input$modEnrich), "Waiting for module detection...")
                        )

                        mod <- input$modEnrich
			modules_ll <- myValues$modules_ll
                        enrichmentDF <- modules_ll[["ae"]][[mod]]

			DT::datatable(enrichmentDF, filter = list(position='top', clear=FALSE), options = list(search = list(regex=TRUE, caseInsensitive=FALSE), scrollX=TRUE))
		})

		###Module Network
                myValues$mod_dynamic <- shiny::reactive({
                        #shiny::validate(
			#	need(!is.null(myValues$iGraph), "No Graph to Plot!")
			#)
			if(is.null(myValues$iGraph))
			return(NULL)

                        progress <- shiny::Progress$new()
                        on.exit(progress$close())

			progress$set(message="Formatting Modules Graph...", value=0)
			mods <- input$mod
			localGxTable <- myValues$gxTable()

			if(is.null(mods))
			return(NULL)

                        if(length(mods)>1){
                                localIGraph <- igraph:::union.igraph(myValues$modules_ll[["ig"]][mods])
                                vAttrNames <- list.vertex.attributes(myValues$modules_ll[["ig"]][[1]])[-1]
                                eAttrNames <- list.edge.attributes(myValues$modules_ll[["ig"]][[1]])
                                #print("vAttrNames")
                                #print(vAttrNames)
                                #print("eAttrNames")
                                #print(eAttrNames)
                                localIGraph <- resolve_ig_union_attrs(localIGraph, vAttrNames, eAttrNames)
                        }else{
                                localIGraph <- myValues$modules_ll[["ig"]][[mods[1]]]
                        }
                        localIGraph <- set_edge_color(localIGraph, localGxTable, pos_cor_color=input$eColP, pos_cor_highlight_color=input$ehColP, neg_cor_color=input$eColN, neg_cor_highlight_color=input$ehColN)
                        modPalette <- myValues$modPalette
                        modVec <- igraph::vertex_attr(localIGraph, name="group")
                        #mods <- unique(modVec)
                        for(mod in mods){
                                idx <- which(modVec %in% mod)
                                igraph::vertex_attr(localIGraph, name="color")[idx] <- modPalette[mod]
                        }
                        localIGraph
                })

                output$modVis <- visNetwork::renderVisNetwork({
                        shiny::validate(
                                #need(!is.null(myValues$iGraph_dynamic()), "Waiting for formatted graph...")
                                need(!is.null(myValues$mod_dynamic()), "Waiting for formatted graph...")
                        )

                        progress <- shiny::Progress$new()
                        on.exit(progress$close())
			
                        #progress$set(message="Extracting Mod", value=0)
			#mods <- input$mod
                        #print("class(mods) : ")
                        #print(class(mods))
                        ##localIGraph <- myValues$modules_ll[["ig"]][input$mod]
                        #if(length(mods)>1){
                        #        localIGraph <- igraph:::union.igraph(myValues$modules_ll[["ig"]][mods])
                        #        vAttrNames <- list.vertex.attributes(myValues$modules_ll[["ig"]][[1]])[-1]
                        #        eAttrNames <- list.edge.attributes(myValues$modules_ll[["ig"]][[1]])
                        #        print("localIGraph : ")
                        #        print(vcount(localIGraph))
                        #        print(ecount(localIGraph))
                        #        print(list.vertex.attributes(localIGraph))
                        #        print(list.edge.attributes(localIGraph))
                        #        localIGraph <- resolve_ig_union_attrs(localIGraph, vAttrNames, eAttrNames)
                        #}else{
                        #        print(V(myValues$modules_ll[["ig"]][[mods[1]]]))
                        #        localIGraph <- myValues$modules_ll[["ig"]][[mods[1]]]
                        #}

                        #modPalette <- myValues$modPalette
                        #print("str(modPalette):")
                        #print(str(modPalette))
                        #modVec <- igraph::vertex_attr(localIGraph, name="group")
                        #mods <- unique(modVec)
                        #for(mod in mods){
                        #        idx <- which(modVec %in% mod)
                        #        igraph::vertex_attr(localIGraph, name="color")[idx] <- modPalette[mod]
                        #}

			progress$set(message="Plotting", value=0)

                        localIGraph <- myValues$mod_dynamic()
                        modVisNetwork <- get_visNetwork(localIGraph, plot_layout=input$graphLayout, vBorderColor=input$vBrdCol, vShape=input$vShape, vFontColor=input$vLblCol, vSize=input$vSize, eWidth=input$eWidth, degDepth=input$dDepth)

			networkImageExportFileName <- paste0(paste(input$mod, collapse="-"), "_", Sys.Date())
			imageType <- "png"
                        bgColor <- input$bgCol
                        modPalette <- myValues$modPalette
			mods <- input$mod
                        addNodes <- data.frame(label=mods, shape="circle", icon.color=modPalette[mods])
                        print("addNodes:")
                        print(addNodes)
			modVisNetwork <- modVisNetwork %>%
			visLegend(addNodes=addNodes, useGroups=FALSE, main="Modules") %>%
			#Add export control
			visExport(type=imageType, name=networkImageExportFileName, label=paste0("Export as ", imageType), background=bgColor, float="right", style=NULL, loadDependencies=TRUE)

                        modVisNetwork
                })

                output$downloadMod <- shiny::downloadHandler(
                        filename = function(){
                                #paste0(input$mod, "_", Sys.Date(), '.', input$modExportFormat)
                                paste0(paste(input$mod, collapse="-"), "_", Sys.Date(), '.', input$modExportFormat)
                        },
                        content = function(con){
                                #write.graph(myValues$modules_ll[["ig"]][input$mod], con, format=input$modExportFormat)
                                write.graph(myValues$mod_dynamic(), con, format=input$modExportFormat)
                        }
                )

                output$downloadModGenes <- shiny::downloadHandler(
                        filename = function(){
                                #paste0(input$mod, "_genes_", Sys.Date(), '.txt')
                                paste0(paste(input$mod, collapse="-"), "_genes_", Sys.Date(), '.txt')
                        },
                        content = function(con){
                                #write.table(V(myValues$modules_ll[["ig"]][input$mod])$name, con, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
                                write.table(V(myValues$mod_dynamic())$name, con, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
                        }
                )

		##Radar Chart
		output$radarChart <- radarchart::renderChartJSRadar({
                        shiny::validate(
                                need(!is.null(myValues$modules_rank_table), "Waiting for module detection...")
                        )
			rank_table <- myValues$modules_rank_table
			if(nrow(rank_table)<3){
				rank_table <- as.data.frame(t(rank_table[,-1]))
			}else{
				rank_table <- rank_table[,c(-1)]
			}
			radarchart::chartJSRadar(scores=rank_table, labs=rownames(rank_table), labelSize=36, polyAlpha=0.3, lineAlpha=0.9, gridLines.lineWidth=10, pointLabels.fontSize=20)
		})

                ##Module Optimization Radar Chart
                observeEvent(input$rChartButton, {
                        if(is.null(myValues$modules_ll))
                        return(NULL)

			progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}
                        on.exit(progress$close())

			progress$set(message="Response Module Optimization", value=0)
			updateProgress(detail="Reorganizing Modules...", value=1/4)

                        optMods <- input$optMod
                        conMods <- input$conMod
                        indMods <- input$indMod

                        localIGraph <- myValues$iGraph
			gene_ll <- myValues$modules_ll[["names"]]
			ae_ll <- myValues$modules_ll[["ae"]]
			rl <- myValues$rl
			rl.c <- myValues$rl.c
			rl.pv <- myValues$rl.pv
			rl.lfc <- myValues$rl.lfc
			edgeRank <- myValues$edgeRank
			rOrgDB <- unlist(strsplit(input$organism, ";"))[1]
			IC_ll <- myValues$IC_ll
                        
                        if(is.null(optMods) && is.null(conMods) && is.null(indMods))
                        return(NULL)

                        morphed_gene_ll <- list()
                        morphed_ae_ll <- list()
                        if(!is.null(optMods)){
                                idx <- which(names(gene_ll) %in% optMods)
                                morphed_gene_ll[["chosen"]] <- unique(as.character(unlist(gene_ll[idx])))
                                optAeDF <- ldply(ae_ll[idx], data.frame)
                                #idx <- which(duplicated(optAeDF$GOID))
                                #if(length(idx)>1){
                                #        optAeDF <- optAeDF[-idx,]
                                #}
                                optAeDF <- optAeDF[,-1]
                                optAeDF <- transform(optAeDF, Score=(-log10(optAeDF$EASE_Score))*optAeDF$Gene_Count)
                                resDF <- plyr::ddply(optAeDF, plyr::.(GOID), plyr::summarize, sum(Score))
                                optAeDF <- plyr::join(resDF, optAeDF, by="GOID", type="left", match="first")
                                morphed_ae_ll[["chosen"]] <- optAeDF
                        }
                        if(!is.null(conMods)){
                                idx <- which(names(gene_ll) %in% conMods)
                                morphed_gene_ll[["contrast"]] <- unique(as.character(unlist(gene_ll[idx])))

                                conAeDF <- ldply(ae_ll[idx], data.frame)
                                #idx <- which(duplicated(conAeDF$GOID))
                                #if(length(idx)>1){
                                #        conAeDF <- conAeDF[-idx,]
                                #}
                                conAeDF <- conAeDF[,-1]
                                conAeDF <- transform(conAeDF, Score=(-log10(conAeDF$EASE_Score))*conAeDF$Gene_Count)
                                resDF <- plyr::ddply(conAeDF, plyr::.(GOID), plyr::summarize, sum(Score))
                                conAeDF <- plyr::join(resDF, conAeDF, by="GOID", type="left", match="first")
                                morphed_ae_ll[["contrast"]] <- conAeDF
                        }
                        if(!is.null(indMods)){
                                idx <- which(names(gene_ll) %in% indMods)
                                if(length(idx)>0){
                                        morphed_gene_ll <- c(morphed_gene_ll, gene_ll[idx])
                                        morphed_ae_ll <- c(morphed_ae_ll, ae_ll[idx])
                                }
                        }

			updateProgress(detail="Getting New Ranks...", value=2/4)
                        #rank_table <- get_rank_table(iGraph=localIGraph, gene_ll=morphed_gene_ll, ae_ll=morphed_ae_ll,  rl=rl, rl.c=rl.c, rl.pv=rl.pv, rl.lfc=rl.lfc, rl.edge=edgeRank, annDB=rOrgDB, IC_ll=IC_ll)
                        rank_table <- get_rank_table(iGraph=localIGraph, gene_ll=morphed_gene_ll, rl=rl, rl.c=rl.c, rl.pv=rl.pv, rl.lfc=rl.lfc, rl.edge=edgeRank) ##Not computing GO sim
                        rank_table <- as.data.frame(t(rank_table[,-1]))

			updateProgress(detail="Summarizing Response GOs...", value=3/4)
			sim_cutoff=0.4
                        optAE <- morphed_ae_ll[["chosen"]]
			#GO_clust_summ_list <- go_summarization(optAE, score_col="EASE_Score", annDB=rOrgDB, IC_ll=IC_ll)
			GO_clust_summ_list <- go_summarization(optAE, score_col="Score", annDB=rOrgDB, IC_ll=IC_ll, log_transform=FALSE, simMeasure=input$selSemSim, treeHeight=input$treeHeight)

			updateProgress(detail="Completed", value=4/4)

                        shinyBS::updateButton(session, "rChartButton", style="success", icon=icon("check-circle"))
			myValues$GO_clust_summ_list <- GO_clust_summ_list
                        myValues$optAnnotation <- rank_table
                        myValues$optAE <- optAE
                        myValues$morphed_gene_ll <- morphed_gene_ll
                        myValues$morphed_ae_ll <- morphed_ae_ll
                        myValues$optModVal <- sort(optMods)
                        myValues$conModVal <- sort(conMods)
                        myValues$indModVal <- sort(indMods)
                })

		output$optRadarChart <- radarchart::renderChartJSRadar({
                        #shiny::validate(
                        #        need(!is.null(myValues$optAE), "Waiting for module optimization results...")
                        #)
                        rank_table <- myValues$optAnnotation
			radarchart::chartJSRadar(scores=rank_table, labs=rownames(rank_table), labelSize=36, polyAlpha=0.3, lineAlpha=0.9, gridLines.lineWidth=5, pointLabels.fontSize=20)
		})

		output$downloadOptRankDT <- shiny::downloadHandler(
			filename = function(){
                                paste("Gene_Tables_Optimized_", Sys.Date(), '.xlsx', sep='')
			},
			content = function(con){
                                localIGraph <- myValues$iGraph
                                rl <- myValues$rl
                                modules_ll <- list() 
                                modules_ll[["names"]] <- myValues$morphed_gene_ll
                                orgDB <- unlist(strsplit(input$organism, ";"))[1]
                                masterDF <- myValues$rankDF()
                                modDF_LL <- get_mod_gene_tables(ig=localIGraph, rankedGenes=rl, modLL=modules_ll, orgDB=orgDB)
                                outLL <- list(master=masterDF)
                                outLL <- c(outLL, modDF_LL)
                                WriteXLS(outLL, ExcelFileName=con, col.names=TRUE, AdjWidth=TRUE, BoldHeaderRow=TRUE)
			}
		)

		##Module GO Similarity Heatmap
		output$simHeatmap <- renderPlot({
                        shiny::validate(
                                need(!is.null(myValues$within_go_sim), "Waiting for module detection...")
                        )
                        shiny::validate(
                                need(nrow(myValues$within_go_sim)>1, "No Heatmap to plot for single module!")
                        )
			#modMethod <- input$modMethod
			#modules_ll <- myValues$modules_ll
			#IC_ll <- myValues$IC_ll
			#within_go_sim <- get_between_GO_sim(modules_ll, modules_ll, prefix1=modMethod, prefix2=modMethod, IC_ll=IC_ll)
			within_go_sim <- myValues$within_go_sim
			simHeatMap <- plot_sim_heatmap(hSim=within_go_sim, rSim=within_go_sim, cSim=within_go_sim)
			simHeatMap <- myValues$simHeatMap
			simHeatMap
		})

		output$downloadSimMat <- shiny::downloadHandler(
			filename = function(){
                                paste0("Modules_Similarity_by_GO_", Sys.Date(), '.txt')
			},
			content = function(con){
				within_go_sim <- myValues$within_go_sim
				write.table(within_go_sim, con, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
			}
		)

		output$tilePlotBP <- shiny::renderPlot({
			shiny::validate(
				need(!is.null(myValues$optAE), "Waiting for module optimization results...")
			)
			shiny::validate(
				need(!is.null(myValues$GO_clust_summ_list$BP), "No results to plot for BP!")
			)
			print("INSIDE TILE PLOT!")
			plot_treemap(myValues$GO_clust_summ_list$BP, ont="BP")
		})

		output$tilePlotCC <- shiny::renderPlot({
			shiny::validate(
				need(!is.null(myValues$optAE), "Waiting for module optimization results...")
			)
			shiny::validate(
				need(!is.null(myValues$GO_clust_summ_list$CC), "No results to plot for CC!")
			)
			print("INSIDE TILE PLOT!")
			plot_treemap(myValues$GO_clust_summ_list$CC, ont="CC")
		})

		output$tilePlotMF <- shiny::renderPlot({
			shiny::validate(
				need(!is.null(myValues$optAE), "Waiting for module optimization results...")
			)
			shiny::validate(
				need(!is.null(myValues$GO_clust_summ_list$MF), "No results to plot for MF!")
			)
			print("INSIDE TILE PLOT!")
			plot_treemap(myValues$GO_clust_summ_list$MF, ont="MF")
		})

                output$downloadTilePlot <- shiny::downloadHandler(
			filename = function(){
				paste("Tileplot_", Sys.Date(), '.PDF', sep='')
			},
			content = function(con){
                                #Disable Warning
                                oldw <- getOption("warn")
                                options(warn = -1)

				tempPDF <- file.path(tempdir(), "tileplot.Rmd")
				file.copy("tileplot.Rmd", tempPDF, overwrite=TRUE)

				params <- list(BP=myValues$GO_clust_summ_list$BP, CC=myValues$GO_clust_summ_list$CC, MF=myValues$GO_clust_summ_list$MF)
				rmarkdown::render(tempPDF, output_file=con,
					params=params,
					envir=new.env(parent=globalenv())
				)

                                #Enable Warning
                                options(warn = oldw)
			}
		)

		myValues$rankDF <- shiny::reactive({
			shiny::validate(
                                need(!is.null(myValues$modules_ll), "Waiting for module detection...")
                        )

			progress <- shiny::Progress$new()
			progress$set(message="Getting Rank Table...", value=0)
			on.exit(progress$close())

                        localIGraph <- myValues$iGraph
			rl <- myValues$rl
			modules_ll <- myValues$modules_ll
			#orgDB <- unlist(strsplit(input$organism, ";"))[2]
			orgDB <- unlist(strsplit(input$organism, ";"))[1]
                        mapped_df <- get_master_gene_table(ig=localIGraph, rankedGenes=rl, modLL=modules_ll, orgDB=orgDB)

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
                                paste("Gene_Tables_", Sys.Date(), '.xlsx', sep='')
			},
			content = function(con){
                                localIGraph <- myValues$iGraph
                                rl <- myValues$rl
                                modules_ll <- myValues$modules_ll
                                orgDB <- unlist(strsplit(input$organism, ";"))[1]
                                masterDF <- myValues$rankDF()
                                modDF_LL <- get_mod_gene_tables(ig=localIGraph, rankedGenes=rl, modLL=modules_ll, orgDB=orgDB)
                                outLL <- list(master=masterDF)
                                outLL <- c(outLL, modDF_LL)
                                WriteXLS(outLL, ExcelFileName=con, col.names=TRUE, AdjWidth=TRUE, BoldHeaderRow=TRUE)
			}
		)

                ###RSERVE SESSION
                #observeEvent(input$rserve_submit, {
                #        if(is.null(myValues$modules_ll)){
                #                shinyjs::info("RServe socket daemon not started. \n\nFailed to find detected modules!!!")
                #                return(NULL)
                #        }

                #        shinyjs::info(paste0("Rserve socket daemon started!!!\n\n",
                #        "NOTE: YOU HAVE STARTED A SOCKET SERVER SESSION. THIS SHINY APP WILL BE UNRESPONSIVE UNTIL THE SERVER IS SHUTDOWN FROM A REMOTE SESSION.\n\n",
                #        "KEEP CALM AND USE RSclient IN ANOTHER R SESSION TO ACCESS THE DAEMON.\n\n",
                #        " --- R CODE ---\n",
                #        "##Connect to INfORM\n",
                #        "c <- RSclient::RSconnect(host = \"localhost\", port = 6311)\n\n",
                #        "##List names of the variables from INfORM\n",
                #        "RSclient::RSeval(c, \"names(shiny_app_vars)\")\n\n",
                #        "##Store the shiny varialbes locally\n",
                #        "shiny_app_vars <- RSclient::RSeval(c, \"shiny_app_vars\")\n\n",
                #        "##Shutdown server and close connection\n",
                #        "RSclient::RSshutdown(c)\n",
                #        "RSclient::RSclose(c)\n",
                #        " --- R CODE ---"
                #        ))
                #        #print("Print on deamon start button:")
                #        #print(ls())
                #        #print(environment())
                #        #print("Global environment:")
                #        #print(myValues$envir)
                #        envir <- myValues$envir
                #        set_shiny_vars(shinyVars=myValues, shinyEnv=envir)
                #        #print(ls(envir))
                #        Rserve::run.Rserve()
                #})

		##DYNAMIC UI WIDGETS
		output$selOrganism <- shiny::renderUI({
			org_ann_libs <- rownames(installed.packages()[grep(".*\\.eg\\.db", rownames(installed.packages())),])
			lapply(org_ann_libs, require, character.only = TRUE)
			org_choices <- sapply(org_ann_libs, function(x){y=unname(packageDescription(x, fields="organism")); val=paste0(x,";",y); names(val) <- y; return(val)}, USE.NAMES=FALSE)
			shiny::selectInput("organism", "Choose Organism", choices=org_choices, multiple=FALSE, selected="org.Hs.eg.db;Homo sapiens")
		})

		output$selPvCol <- renderUI({
			selectInput("pvCol", "P.Value Column", choices=myValues$dgxColChoices())
		})

		output$selLfcCol <- renderUI({
			selectInput("lfcCol", "LogFC Column", choices=myValues$dgxColChoices())
		})

		output$selSep <- renderUI({
			selectInput("sepS", "Column Seperator", choices=myValues$sepChoices, selected=myValues$sepChoices[1])
		})

		output$selQuote <- renderUI({
			selectInput("quote", "Quotes", choices=myValues$quoteChoices, selected=myValues$quoteChoices[1])
		})

                output$selGxSep <- renderUI({
			selectInput("gxSepS", "Column Seperator", choices=myValues$sepChoices, selected=myValues$sepChoices[1])
		})

		output$selGxQuote <- renderUI({
			selectInput("gxQuote", "Quotes", choices=myValues$quoteChoices, selected=myValues$quoteChoices[1])
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

		output$selMod <- shiny::renderUI({
			if(is.null(myValues$modules_ll)){
				tmpChoices <- NA
				sel <- NA
			}else{
				tmpChoices <- names(myValues$modules_ll[["ig"]])
				sel <- tmpChoices[1]
			}
			shiny::selectInput("mod", "Select Module", choices=tmpChoices, multiple=TRUE, selected=sel)
		})

                output$selModEnrich <- shiny::renderUI({
			if(is.null(myValues$modules_ll)){
				tmpChoices <- NA
				sel <- NA
			}else{
				tmpChoices <- names(myValues$modules_ll[["ig"]])
				sel <- tmpChoices[1]
			}
			shiny::selectInput("modEnrich", "Select Module", choices=tmpChoices, multiple=FALSE, selected=sel)
		})

		output$selFocusMod <- shiny::renderUI({
			if(is.null(myValues$modules_ll)){
				tmpChoices <- NA
				sel <- NA
			}else{
				tmpChoices <- c("ALL", names(myValues$modules_ll[["ig"]]))
				sel <- tmpChoices[1]
			}
			shiny::selectInput("focusMod", "Focus Module", choices=tmpChoices, multiple=FALSE, selected=sel)
		})

                output$selOptMod <- shiny::renderUI({
			if(is.null(myValues$modules_ll)){
				tmpChoices <- NA
				#sel <- NA
			}else{
				#tmpChoices <- c("none", names(myValues$modules_ll[["ig"]]))
				#sel <- tmpChoices[1]
				tmpChoices <- names(myValues$modules_ll[["ig"]])
			}
			shiny::selectInput("optMod", "Select Modules of Interest", choices=tmpChoices, multiple=TRUE, selected=NULL)
		})

                output$selConMod <- shiny::renderUI({
			if(is.null(myValues$modules_ll)){
				tmpChoices <- NA
				#sel <- NA
			}else{
				#tmpChoices <- c("none", names(myValues$modules_ll[["ig"]]))
				#sel <- tmpChoices[1]
				tmpChoices <- names(myValues$modules_ll[["ig"]])
			}
			shiny::selectInput("conMod", "Select Modules for Contrast", choices=tmpChoices, multiple=TRUE, selected=NULL)
		})

                output$selIndMod <- shiny::renderUI({
			if(is.null(myValues$modules_ll)){
				tmpChoices <- NA
				#sel <- NA
			}else{
				#tmpChoices <- c("none", names(myValues$modules_ll[["ig"]]))
				#sel <- tmpChoices[1]
				tmpChoices <- names(myValues$modules_ll[["ig"]])
			}
			shiny::selectInput("indMod", "Select Modules as Independents", choices=tmpChoices, multiple=TRUE, selected=NULL)
		})

                shiny::observe({
                        #Update detect module button if options are change sin UI
                        if(isTRUE(input$modMethod!=myValues$modMethodVal) || isTRUE(input$minModSize!=myValues$minModSizeVali)){
                                toDetect <- TRUE
                        }else{
                                toDetect <- FALSE
                        }
                        if(toDetect==TRUE){
                                shinyBS::updateButton(session, "runDetect", style="danger", icon=icon("exclamation-circle"))
                        }else{
                                shinyBS::updateButton(session, "runDetect", style="success", icon=icon("check-circle"))
                        }

			#Logic for selectInputs in module optimization
			if(!is.null(myValues$modules_ll)){
				print("%%% IN OBSERVE GLOBAL %%%") 
				print("%%% OBSERVE MODS %%%") 
				tmpChoices <- names(myValues$modules_ll[["ig"]])
				modsAll <- unique(c(input$optMod, input$conMod, input$indMod))
                                modChk <- which(modsAll %in% tmpChoices)
                                #if(length(modChk)>0){
                                if(length(modChk)==length(modsAll)){
                                        print(tmpChoices)
                                        modsOptCon <- unique(c(input$optMod, input$conMod))
                                        modsOptInd <- unique(c(input$optMod, input$indMod))
                                        modsConInd <- unique(c(input$conMod, input$indMod))
                                        print(modsOptCon)
                                        print(modsOptInd)
                                        print(modsConInd)
                                        if(is.null(modsOptCon)){
                                                tmpChoicesInd <- tmpChoices
                                        }else{
                                                idx <- which(tmpChoices %in% modsOptCon)
                                                if(length(idx)>0){
                                                        tmpChoicesInd <- tmpChoices[-idx]
                                                }
                                        }
                                        if(is.null(modsOptInd)){
                                                tmpChoicesCon <- tmpChoices
                                        }else{
                                                idx <- which(tmpChoices %in% modsOptInd)
                                                if(length(idx)>0){
                                                        tmpChoicesCon <- tmpChoices[-idx]
                                                }
                                        }
                                        if(is.null(modsConInd)){
                                                tmpChoicesOpt <- tmpChoices
                                        }else{
                                                idx <- which(tmpChoices %in% modsConInd)
                                                if(length(idx)>0){
                                                        tmpChoicesOpt <- tmpChoices[-idx]
                                                }
                                        }

                                        if(!(length(myValues$optChoicesVal)==length(tmpChoicesOpt) && all(myValues$optChoicesVal %in% sort(tmpChoicesOpt)))){
                                                updateSelectInput(session, "optMod", choices=tmpChoicesOpt, selected=input$optMod)
                                        }
                                        if(!(length(myValues$conChoicesVal)==length(tmpChoicesCon) && all(myValues$conChoicesVal %in% sort(tmpChoicesCon)))){
                                                updateSelectInput(session, "conMod", choices=tmpChoicesCon, selected=input$conMod)
                                        }
                                        if(!(length(myValues$indChoicesVal)==length(tmpChoicesInd) && all(myValues$indChoicesVal %in% sort(tmpChoicesInd)))){
                                                updateSelectInput(session, "indMod", choices=tmpChoicesInd, selected=input$indMod)
                                        }
                                        myValues$optChoicesVal <- sort(tmpChoicesOpt)
                                        myValues$conChoicesVal <- sort(tmpChoicesCon)
                                        myValues$indChoicesVal <- sort(tmpChoicesInd)

                                        toPlot <- FALSE
                                        if(!(length(myValues$optModVal)==length(input$optMod) && all(myValues$optModVal %in% sort(input$optMod)))){
                                                toPlot <- TRUE
                                        }
                                        if(!(length(myValues$conModVal)==length(input$conMod) && all(myValues$conModVal %in% sort(input$conMod)))){
                                                toPlot <- TRUE
                                        }
                                        if(!(length(myValues$indModVal)==length(input$indMod) && all(myValues$indModVal %in% sort(input$indMod)))){
                                                toPlot <- TRUE
                                        }
                                        if(length(input$optMod)==0 && length(input$conMod)==0 && length(input$indMod)==0){
                                                toPlot <- TRUE
                                        }

                                        if(toPlot==TRUE){
                                                shinyBS::updateButton(session, "rChartButton", style="danger", icon=icon("exclamation-circle"))
                                        }else{
                                                shinyBS::updateButton(session, "rChartButton", style="success", icon=icon("check-circle"))
                                        }
                                }
			}
		})

                #observeEvent(input$optMod, {
		#	if(is.null(myValues$modules_ll))
                #        return(NULL)
                #     
		#	print("%%% IN OBSERVE OPT MOD INPUT %%%") 
                #        modsOptCon <- unique(c(input$optMod, input$conMod))
                #        modsOptInd <- unique(c(input$optMod, input$indMod))
                #        tmpChoices <- names(myValues$modules_ll[["ig"]])
		#	if(is.null(modsOptCon)){
		#		tmpChoicesInd <- tmpChoices
		#	}else{
		#		idx <- which(tmpChoices %in% modsOptCon)
		#		tmpChoicesInd <- tmpChoices[-idx]
		#	}
		#	if(is.null(modsOptInd)){
		#		tmpChoicesCon <- tmpChoices
		#	}else{
		#		idx <- which(tmpChoices %in% modsOptInd)
		#		tmpChoicesCon <- tmpChoices[-idx]
		#	}
                #        updateSelectInput(session, "conMod", choices=tmpChoicesCon, selected=input$conMod)
                #        updateSelectInput(session, "indMod", choices=tmpChoicesInd, selected=input$indMod)
                #})

                #observeEvent(input$conMod, {
		#	if(is.null(myValues$modules_ll))
                #        return(NULL)
                #      
		#	print("%%% IN OBSERVE CON MOD INPUT %%%") 
                #        modsOptCon <- unique(c(input$optMod, input$conMod))
                #        modsConInd <- unique(c(input$conMod, input$indMod))
                #        tmpChoices <- names(myValues$modules_ll[["ig"]])
		#	if(is.null(modsOptCon)){
		#		tmpChoicesInd <- tmpChoices
		#	}else{
		#		idx <- which(tmpChoices %in% modsOptCon)
		#		tmpChoicesInd <- tmpChoices[-idx]
		#	}
		#	if(is.null(modsConInd)){
		#		tmpChoicesOpt <- tmpChoices
		#	}else{
		#		idx <- which(tmpChoices %in% modsConInd)
		#		tmpChoicesOpt <- tmpChoices[-idx]
		#	}
                #        updateSelectInput(session, "optMod", choices=tmpChoicesOpt, selected=input$optMod)
                #        updateSelectInput(session, "indMod", choices=tmpChoicesInd, selected=input$indMod)
                #})

                #observeEvent(input$indMod, {
		#	if(is.null(myValues$modules_ll))
                #        return(NULL)
                #      
		#	print("%%% IN OBSERVE IND MOD INPUT %%%") 
                #        modsOptInd <- unique(c(input$optMod, input$indMod))
                #        modsConInd <- unique(c(input$conMod, input$indMod))
                #        tmpChoices <- names(myValues$modules_ll[["ig"]])
		#	if(is.null(modsOptInd)){
		#		tmpChoicesCon <- tmpChoices
		#	}else{
		#		idx <- which(tmpChoices %in% modsOptInd)
		#		tmpChoicesCon <- tmpChoices[-idx]
		#	}
		#	if(is.null(modsConInd)){
		#		tmpChoicesOpt <- tmpChoices
		#	}else{
		#		idx <- which(tmpChoices %in% modsConInd)
		#		tmpChoicesOpt <- tmpChoices[-idx]
		#	}
                #        updateSelectInput(session, "optMod", choices=tmpChoicesOpt, selected=input$optMod)
                #        updateSelectInput(session, "conMod", choices=tmpChoicesCon, selected=input$conMod)
                #})

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
				"<p>Default node color represents the positive and negative trend of gene by differential expression. ",
				"Edge color represents the positive or negative correlation between the genes.</p> ",
				"<p>Detected modules can be highlighted by checking box 'Show Modules', node color is updated to represent the module. ",
				"Color association of modules is presented in the legend, graph view can be focused to a specific module by selection the module of interest from 'Focus Module'.</p> ",
				"<p>Plot aesthetics and representation can be changed from the 'Show/Hide Aesthetics Options' section, click on the label to see the panel of options. ",
				"This panel provides option to change graph layout, specify percentile threshold to denote positively and negatively associated nodes. ",
				"Aesthetics can be altered to show gene rank as node color, increase/decrease node label size, alter the shape of the node and ",
				"specify edge width. Cutomization options for node color, edge color, background color, node border color are also provided.</p>",
				"<p>Scroll up/down within the plot area to zoom in/zoom out, click and drag within the plot area to move the plot. ",
				"Hover over a node to see the node name, clicking a node opens NCBI ENTREZ GENE database record of the gene represented by that node.</p> "
			),
			placement="right",
			trigger="focus",
			#options=list(container="body", width="70%", "max-width"="70%")
			options=list(container="#main_info")
		)

		shinyBS::addPopover(session, "modVis_info", "Module Info", content=paste0("<p>Graphical representation of the modules detected from the main network ", 
				"by using one the community detection methods (DEFAULT:Walktrap) and minimum module size (DEFAULT:10). </p>",
				"<p>These parameters can be set from the 'Show/hide Advanced Options' panel. ",
				"Selected modules are highlighted with specific colors, module color association is presented as a legend. ",
				"General aesthetics customization and layout is same as main graph.</p>"
			),
			placement="right",
			trigger="focus",
			options=list(container="#mod_info")
		)

		shinyBS::addPopover(session, "enrichmentDT_info", "Annotation Enrichment Info", content=paste0("<p>Annotation enrichment is performed by Fisher's Exact Test on the set of genes from each detected module. </p>",
				"<p>Enrichment result is filtered by EASE score of '0.05'.</p>",
				"<p>Enrichment table for a module can be displayed by selecting the module of interest from the selection box 'Select Module'. </p>",
				"<p>Enrichment table for all the modules can be exported as a single spreasheet from the 'Download & Export' section. </p>"
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

		shinyBS::addPopover(session, "simHeatmap_info", "GO based Similarity Heatmap Info", content=paste0("<p>Similarity between detected modules is computed as the Jaccard Index for set of significant Gene Ontology terms enriched by the genes constituting the module. </p>",
				"<p>The Jaccard Index is assisted with Semantic Similarity between GO terms, so each individual match is not exact but is signigficantly close related terms based on Sematic Similarity. </p>",
				"<p>The Heatmap is generated by converting the similarity matrix into a distance matrix (1-sim), the computed distance matrix is also used to draw the cluster dendogram on rows and columns. </p>",
				"<p>The color scale ranges from red to yellow, meaning the most similar are red and most distant are yellow. </p>"
			),
			placement="right",
			trigger="focus",
			options=list(container="#heatmap_info")
		)

		shinyBS::addPopover(session, "radarChart_info", "Radar Chart Info", content=paste0("<p>Radar chart represents the multivariate information in 2-dimensional plane. </p>",
				"<p>This radar chart plots an axis for each detected module and values for each property related to the module are plotted on their axis, points for each property have a specific color and are joined to points on other axis creating a color polygon for that property. </p>",
				"<p>In this radar chart 'rank_centrality' represent the median rank of the genes by their centrality properties in the whole network. 'rank_pv' represents the median rank of the genes by differential analysis P.Value (EASE_Score) in the whole gene list. 'rank_lfc' represents the median rank of the genes by differential analysis log Fold Change in the whole gene list. 'rank_edge' represents the median edge rank of the edges in the modules in the overall edge ranks of the while network. 'size' represents the normalized number of nodes in the module vs the total number of nodes in the whole network. 'GO' represents the similairty within the set of significantly enriched GO terms in the module as semantic similarity assisted Jaccrd Index. </p>",
				"<p>User can hide/display the points for a property by clicking on the property name in the legend. </p>"
			),
			placement="right",
			trigger="focus",
			options=list(container="#radar_info")
		)

		shinyBS::addPopover(session, "optRadarChart_info", "Optimization Radar Chart Info", content=paste0("<p>Radar chart represents the multivariate information in 2-dimensional plane. </p>",
				"<p>This radar chart plots an axis for each property and values representing each of the 'chosen', 'contrast' and individual modules are plotted point on their axis, points for each module have a specific color and are joined to points on other axis creating a color polygon for that module. </p>",
				"<p>In this radar chart axis 'rank_centrality' represent the median rank of the genes by their centrality properties in the whole network. 'rank_pv' represents the median rank of the genes by differential analysis P.Value (EASE_Score) in the whole gene list. 'rank_lfc' represents the median rank of the genes by differential analysis log Fold Change in the whole gene list. 'rank_edge' represents the median edge rank of the edges in the modules in the overall edge ranks of the while network. 'size' represents the normalized number of nodes in the module vs the total number of nodes in the whole network. 'GO' represents the similairty within the set of significantly enriched GO terms in the module as semantic similarity assisted Jaccrd Index. </p>",
				"<p>The module 'chosen' represents the grouping of selected modules, module 'contrast' represents the grouping of modules selected as contrast and selected individual modules are plotted with their original labels. </p>",
				"<p>User can hide/display the points for a module by clicking on the property name in the legend. </p>"
			),
			placement="right",
			trigger="focus",
			options=list(container="#opt_radar_info")
		)

		shinyBS::addPopover(session, "rankDT_info", "Rank Table Info", content=paste0("<p>Table containing the final rank of the gene from the whole network inferred by INfORM, ",
				"HGNC Gene Symbol of the gene, unique identifier from the NCBI ENTREZ Gene database, descriptive name of the gene, user provided differential expression score, ",
				"P.Value and log fold change of the differentially expressed genes provided by the user, ",
				"membership of gene in modules extracted by using community detection method. </p>",
				"<p>Genes not belonging to any module over the minimum size threshold are marked as 'orphan' members. </p>",
				"<p>Overall rank table along with module specific rank tables can be exported as a single spreadsheet. </p>"
			),
			placement="right",
			trigger="focus",
			options=list(container="#rt_info")
		)
	}
)
