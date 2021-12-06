suppressMessages(library(shiny))
suppressMessages(library(shinyjs))
suppressMessages(library(colourpicker))
suppressMessages(library(shinyBS))
suppressMessages(library(shinydashboard))
suppressMessages(library(DT))
suppressMessages(library(visNetwork))
suppressMessages(library(radarchart))

jsCode <- 'shinyjs.plot_bg_col = function(col){
        var el = $("." + "vPlot");
        el.css("background-color", col);
}'

appCSS <- "
	.nav-tabs-custom .nav-tabs { border-top-left-radius: 3px; border-top-right-radius: 3px; }
	.nav-tabs-custom .nav-tabs li.active { border-left: 1px #3c8dbc solid; border-right: 1px #3c8dbc solid; }
	.nav-tabs-custom .nav-tabs li a {
		border-left: 1px #3c8dbc solid;
		border-right: 1px #3c8dbc solid;
		border-top: 1px #3c8dbc solid;
		border-radius: 4px 4px 0 0;
	}
	.nav-tabs-custom .nav-tabs li { margin-right: 3px }
	.nav-tabs-custom .tab-content { border: 1px #3c8dbc solid; }

	//For plot height
	//.shiny-plot-output { height:100vh !important; }

	//For plot loading
	//.plot-container { position: relative; }
	//#loading-spinner {
	//	position: absolute;
	//	left: 50%;
	//	top: 50%;
	//	z-index: -1;
	//	margin-top: -33px;  /* half of the spinner's height */
	//	margin-left: -33px; /* half of the spinner's width */
	//}
	//#hcPlot.recalculating { z-index: -2; }
	//#.shiny-plot-output .recalculating { z-index: -2; }
"

dashboardPage(
	dashboardHeader(title="INfORM - Inference of NetwOrk Response Module", titleWidth="25%"),
	dashboardSidebar(disable=TRUE),
	dashboardBody(
		tags$head(tags$style(HTML('
		   div.wrap .popover {
			width: 700px;
			max-width: 100%;
		    }
		'))),
		useShinyjs(),
		extendShinyjs(text = jsCode, functions = "plot_bg_col"),
		inlineCSS(appCSS),
                shinyBS::bsModal("importGxModal", "Import Gene Expression Table", "import_gx_submit", size="large",
			fluidRow(
				column(3,
					fileInput("gx", label="File")
				),column(3,
					uiOutput("selGxSep")
				),column(3,
					textInput("gxSepT", "Other Seperator", value=":")
				),column(3,
					uiOutput("selGxQuote")
				)
			),fluidRow(
				column(12, align="right",
					shinyBS::bsButton("upload_gx_submit", label="Import", style="info", icon=icon("hand-o-right"))
				)
			)
		),
		shinyBS::bsModal("importDgxModal", "Import Differential Gene Expression Table", "import_dgx_submit", size="large",
			fluidRow(
				column(3,
					fileInput("dgx", label="File")
				),column(3,
					uiOutput("selSep")
				),column(3,
					textInput("sepT", "Other Seperator", value=":")
				),column(3,
					uiOutput("selQuote")
				)
			),fluidRow(
				column(1,
					actionButton("load_dgx_submit", "Load")
				)
			),hr(),
			fluidRow(
				column(4,
					uiOutput("selPvCol")
				),column(4,
					uiOutput("selLfcCol")
				)
			),fluidRow(
				column(12, align="right",
					shinyBS::bsButton("upload_dgx_submit", label="Import", style="info", icon=icon("hand-o-right"))
				)
			),hr(),
			fluidRow(
				column(12,
					DT::dataTableOutput("dgxDT")
				)
			)
		),
		fluidRow(
			box(
				title="Upload", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, width=12,
				fluidRow(
			#		column(12,
			#		    h4("Upload Gene Expression Data")
			#		)
			#	),fluidRow(
					column(12,
						h4("Upload Data"),
						wellPanel(
							fluidRow(
								column(4,
									#fileInput(inputId="gx", label="Gene Expression Table")
                                                                        div(id="outerGxDiv", class="form-group shiny-input-container",
										tags$label("Gene Expression Table"),
										div(id="innerGxDiv", class="input-group",
											div(id="gxBtnDiv", class="input-group-btn",
												shinyBS::bsButton("import_gx_submit", label="Upload", style="danger", icon=icon("exclamation-circle"))
											),
											tags$input(id="gxUploadDisp", class="form-control", placeholder="No file uploaded", readonly="readonly", type="text")
										)
									),
									shinyBS::bsTooltip("import_gx_submit", "Launch a graphical window, to configure import of gene expression data from a file!", placement="bottom")
								),column(4,
									#fileInput(inputId="dgx", label="Differential Gene Expression Table"),
								#),column(2,
									#shinyBS::bsButton("import_dgx_submit", label="Import Differential Gene Expression Table", style="danger", icon=icon("exclamation-circle")),
									div(id="outerDgxDiv", class="form-group shiny-input-container",
										tags$label("Differential Gene Expression Table"),
										div(id="innerDgxDiv", class="input-group",
											div(id="dgxBtnDiv", class="input-group-btn",
												shinyBS::bsButton("import_dgx_submit", label="Upload", style="danger", icon=icon("exclamation-circle"))
											),
											tags$input(id="dgxUploadDisp", class="form-control", placeholder="No file uploaded", readonly="readonly", type="text")
										)
									),
									shinyBS::bsTooltip("import_dgx_submit", "Launch a graphical window, to configure import of differential expression data from a file!", placement="bottom")
								),column(4,
                                                                        selectInput("ensembleStrat", "Choose Ensemble Strategy",
                                                                                choices=c("MINET Inference"="minet",
                                                                                        "User Matrices"="user",
                                                                                        "MINET Inference + User Matrices"="minet+user"
                                                                                ),
                                                                                multiple=FALSE,
                                                                                selected="minet"
                                                                        )
								)
							),fluidRow(
								column(12,hidden(div(id="matOptions",
                                                                        fluidRow(
                                                                                column(4,
                                                                                        fileInput(inputId="mat", label="User Provided Matrices", multiple=TRUE)
                                                                                ),column(4,
                                                                                        selectInput("matWeights", "Edge Weight Property",
                                                                                                choices=c("Ranks"="rank",
                                                                                                        "Scores"="score"
                                                                                                ),
                                                                                                multiple=FALSE,
                                                                                                selected="rank"
                                                                                        )
                                                                                )
                                                                        )
                                                                )))
							)
						)
					)
				),fluidRow(
					column(6,
						#checkboxInput("checkUnconnected", "Remove Unconnected Genes", value=FALSE),
						uiOutput("selOrganism"),
						hr(),
						actionButton("runINfORM", "Run INfORM"),
						bsTooltip("runINfORM", "Run INfORM with provided Gene Expression Table and Differential Expression Table!", placement="right")
					),column(6,
						uiOutput("selCores")
					)
				),fluidRow(
					column(12,
						a(id = "showAdvancedOptions", "Show/Hide Advanced Options"),
						hidden(div(id="advancedOptions",
                                                        fluidRow(column(6,
                                                                h4("Select Parameters to Run MINET (Mutual Information NETworks)"),
                                                                wellPanel(fluidRow(column(12,
                                                                        #uiOutput("selMethod"),
                                                                        #uiOutput("selEst"),
                                                                        #uiOutput("selDisc")
                                                                        selectInput("method", "Select Inference Algorithm",
                                                                                choices=c("clr",
                                                                                        "aracne",
                                                                                        "mrnet",
                                                                                        "mrnetb"
                                                                                ),
                                                                                multiple=TRUE,
                                                                                selected=c("clr",
                                                                                        "aracne",
                                                                                        "mrnet",
                                                                                        "mrnetb"
                                                                                )
                                                                        ),
                                                                        selectInput("est", "Select Correlation",
                                                                                choices=c("pearson",
                                                                                        "spearman",
                                                                                        "kendall",
                                                                                        "mi.empirical",
                                                                                        "mi.mm",
                                                                                        "mi.shrink",
                                                                                        "mi.sg"
                                                                                ),
                                                                                multiple=TRUE,
                                                                                selected=c("pearson",
                                                                                        "spearman",
                                                                                        "kendall",
                                                                                        "mi.empirical",
                                                                                        "mi.mm",
                                                                                        "mi.shrink",
                                                                                        "mi.sg"
                                                                                )
                                                                        ),
                                                                        selectInput("disc", "Select Discretization Method",
                                                                                choices=c("none",
                                                                                        "equalfreq",
                                                                                        "equalwidth",
                                                                                        "globalequalwidth"
                                                                                ),
                                                                                multiple=TRUE,
                                                                                selected=c("none",
                                                                                        "equalfreq",
                                                                                        "equalwidth",
                                                                                        "globalequalwidth"
                                                                                )
                                                                        )
                                                                ))),
                                                                h4("Edge Selection For Final Ensemble"),
                                                                wellPanel(fluidRow(column(6,
                                                                        selectInput("selEdge", "Edge Selection Strategy",
                                                                                choices=c("default", "top"),
                                                                                multiple=FALSE,
                                                                                selected="default",
                                                                        )
                                                                ),column(6,
                                                                        numericInput("topCutOff", "Top 'n' Percent Edges",
                                                                                value=10,
                                                                                min=1,
                                                                                max=100,
                                                                                step=1
                                                                        )
                                                                ))),
                                                                h4("Module Detection"),
                                                                wellPanel(fluidRow(column(12,
                                                                        fluidRow(
                                                                                column(4,
                                                                                        selectInput("modMethod", "Method",
                                                                                                choices=c(
                                                                                                        "Walktrap"="walktrap",
                                                                                                        "Spinglass"="spinglass",
                                                                                                        "Louvain"="louvain",
                                                                                                        "Greedy"="greedy"
                                                                                                ),
                                                                                                multiple=FALSE,
                                                                                                selected="walktrap"
                                                                                        )
                                                                                ),column(4,
                                                                                        numericInput("minModSize", "Minimum Module Size",
                                                                                                value=10,
                                                                                                min=1,
                                                                                                max=NA,
                                                                                                step=1
                                                                                        )
                                                                                ),column(4,
                                                                                        div(id="outerDgxDiv", class="form-group shiny-input-container",
                                                                                                tags$label(" "),
                                                                                                div(id="innerDgxDiv", class="input-group",
                                                                                                        shinyBS::bsButton("runDetect", label="Detect Modules", style="danger", icon=icon("exclamation-circle"))
                                                                                                )
                                                                                        )
                                                                                )
                                                                        )
                                                                )))
                                                        ),column(6,
                                                                h4("Intra Inference Algorithm Consensus"),
                                                                wellPanel(fluidRow(column(12,
                                                                        selectInput("summType", "Summarization Score",
                                                                                choices=c("median",
                                                                                        "mean",
                                                                                        "max"
                                                                                ),
                                                                                multiple=FALSE,
                                                                                selected="median"
                                                                        )
                                                                ))),
                                                                h4("Inter Inference Algorithm Consensus"),
                                                                wellPanel(fluidRow(column(12,
                                                                        selectInput("bordaScore", "Rank Aggregation Score",
                                                                                choices=c("median",
                                                                                        "mean",
                                                                                        "geo.mean",
                                                                                        "l2norm"
                                                                                ),
                                                                                multiple=FALSE,
                                                                                selected="median"
                                                                        )
                                                                ))),
                                                                h4("Configure Gene Ranking"),
                                                                wellPanel(fluidRow(column(12,
                                                                        uiOutput("rankAttr")
                                                                ))),
                                                                h4("Gene Ontology Annotation Options"),
                                                                wellPanel(fluidRow(column(4,
                                                                        numericInput("sigPval", "Significanct P.Value",
                                                                                value=0.05,
                                                                                min=0.0000001,
                                                                                max=0.1,
                                                                                step=0.000001
                                                                        )
                                                                ),column(4,
                                                                        selectInput("selSemSim", "Semantic Similarity Method",
                                                                                choices=c("Rel", "Resnik", "Lin", "Jiang", "Wang"),
                                                                                multiple=FALSE,
                                                                                selected="Rel",
                                                                        )
                                                                ),column(4,
                                                                        numericInput("treeHeight", "Cluster Tree Cut Height",
                                                                                value=0.9,
                                                                                min=0.1,
                                                                                max=1,
                                                                                step=0.05
                                                                        )
                                                                )))
                                                        ))
						))
					)
				)
			)
		),fluidRow(column(12,
			tabBox(
				id="display", title="Display Area", width=12,
				tabPanel(value="corMat_display", title="Adjacency Matrix",
					fluidRow(
						valueBoxOutput('totalGeneBox', width=3),
						valueBoxOutput('connGeneBox', width=3),
						valueBoxOutput('unGeneBox', width=3),
						valueBoxOutput('netDensityBox', width=3)
					),fluidRow(
						column(12,
							fluidRow(
								column(4,
									downloadButton('downloadMat', label='Download Adjacency Matrix'),
									hr()
								)
							),fluidRow(
								column(4,
									downloadButton('downloadUnGenes', label='Download Un-Connected Genes List')
								)
							),fluidRow(
								column(4,
									tableOutput('unGeneList')
								)
							)
						)
					)
				),tabPanel(value="net_display", title="Network",
					fluidRow(column(12,
					a(id = "showColorOptionsMain", "Show/Hide Aesthetics Options"),
					hidden(div(id="colorOptionsMain",
					       fluidRow(
							column(6,
								h4("Node Association Colors"),
								radioButtons("vColType", NULL, choices=c("Color by Score"="score", "Color Gradient by Rank"="rank"), selected="score", inline=TRUE),
								wellPanel(
									div(id="colByScore",
										h4("Color Nodes by Score"),
										fluidRow(
											column(4,colourpicker::colourInput("vColP", "Node Color +ve", value="#FA8072")),
											column(4,colourpicker::colourInput("vColN", "Node Color -ve", value="#ADD8E6")),
											column(4,colourpicker::colourInput("vColD", "Node Color Default", value="#D3D3D3"))
										),fluidRow(
											column(4,colourpicker::colourInput("vhColP", "Node Highlight Color +ve", value="#FF0000")),
											column(4,colourpicker::colourInput("vhColN", "Node Highlight Color -ve", value="#4169E1")),
											column(4,colourpicker::colourInput("vhColD", "Node Highlight Color Default", value="#A9A9A9"))
										),
										h4("Score For Node Association"),
										fluidRow(
											column(6, selectInput("posPerc", "Precentile Threshold: +ve Association", choices=seq(from=0.01, to=0.99, by=0.01), multiple=FALSE, selected="0.5")),
											column(6, selectInput("negPerc", "Percentile Threshold: -ve Association", choices=seq(from=0.01, to=0.99, by=0.01), multiple=FALSE, selected="0.49"))
										)
									),
									hidden(
									       div(id="colByRank",
											h4("Color Gradient for Nodes by Rank"),
											fluidRow(
												column(4,colourpicker::colourInput("vColT", "Node Color Top Ranked", value="#FF0000")),
												column(4,colourpicker::colourInput("vColB", "Node Color Bottom Ranked", value="#FFFF00"))
											)
									       )
									)
								),
								wellPanel(
									fluidRow(
										column(12,
											h4("Layout"),
											wellPanel(
												selectInput("graphLayout", label="2D Graph Layout",
													choices=c(
														"Nicely"="nicely",
														"Random"="randomly",
														"Circle"="in_circle",
														"Sphere"="on_sphere",
														"Star"="as_star",
														"Tree"="as_tree",
														"Grid"="on_grid",
														"Fruchterman Reingold"="with_fr",
														"Kamada Kawai"="with_kk",
														"DRL"="with_drl",
														"GEM"="with_gem",
														"Graphopt"="with_graphopt",
														"LGL"="with_lgl",
														"MDS"="with_mds"
													), multiple=FALSE, selected="nicely")
											)
										)
									)
								)
							),
							column(6,
								h4("Edge Association Colors"),
								wellPanel(
									fluidRow(
										column(6,colourpicker::colourInput("eColP", "Edge Color +ve", value="#FA8072")),
										column(6,colourpicker::colourInput("eColN", "Edge Color -ve", value="#ADD8E6"))
									),fluidRow(
										column(6,colourpicker::colourInput("ehColP", "Edge Highlight Color +ve", value="#FF0000")),
										column(6,colourpicker::colourInput("ehColN", "Edge Highlight Color -ve", value="#4169E1"))#,
									)
								),
								wellPanel(
									h4("Plot Aesthetics"),
									fluidRow(
										column(6,selectInput("vShape", "Node Shape", choices=c("circle","ellipse","database","box","text"), multiple=FALSE, selected="circle")),
										column(6,colourpicker::colourInput("vBrdCol", "Node Border Color", value="#000000"))
									),fluidRow(
										column(6,numericInput("vSize", "Node Label Size", value="20")),
										column(6,colourpicker::colourInput("vLblCol", "Node Label Color", value="#343434"))
									),fluidRow(
										column(6,numericInput("eWidth", "Edge Width", value="5", width="100%")),
										column(6,colourpicker::colourInput("bgCol", "Background Color", value="#FFFFFF"))
									),fluidRow(
										column(6,numericInput("dDepth", "Depth of Highlighting Nearest Nodes", value="1", width="100%"))
									)
								)
							)
						)
					)),
					tabBox(
						id="net_panel", title="", width=12,
						tabPanel(value="main", title="Main",
							fluidRow(column(12,
								box(title="Download Options", status="success", solidHeader=TRUE, collapsible=TRUE, collapsed=TRUE, width=4,
									wellPanel(
										selectInput("graphExportFormat", NULL, choices=c("edgelist", "pajek", "ncol", "lgl", "graphml", "dimacs", "gml", "dot", "leda"), multiple=FALSE, selected="dot"),
										downloadButton('downloadGraph', label='Export Graph')
									),
									downloadButton('downloadGraphGenes', label='Download Graph Gene List')
								),
								div(id="main_info", class="wrap", actionButton("IGraphVis_info", "Main Graph Info", icon=icon("info-circle")))
							)),fluidRow(
								column(1,
                                                                        checkboxInput(inputId="chkGrpCol", label="Show Modules", value=FALSE)
								),column(2, align="left",
									uiOutput("selFocusMod")
								)
							),fluidRow(
								column(12,
									wellPanel(
										div(id="mainPlot", class="vPlot",
											#visNetworkOutput('IGraphVis', width="100%", height="1000px")
											visNetworkOutput('IGraphVis', width="90%", height="800px")
										)
									)
								)
							)
						),
						tabPanel(value="mods", title="Modules",
							fluidRow(column(12,
								box(title="Download Options", status="success", solidHeader=TRUE, collapsible=TRUE, collapsed=TRUE, width=4,
									wellPanel(
										selectInput("modExportFormat", NULL, choices=c("edgelist", "pajek", "ncol", "lgl", "graphml", "dimacs", "gml", "dot", "leda"), multiple=FALSE, selected="dot", width="50%"),
										downloadButton('downloadMod', label='Export Module')
									),
									downloadButton('downloadModGenes', label='Download Module Gene List')
								),
								div(id="mod_info", class="wrap", actionButton("modVis_info", "Module Info", icon=icon("info-circle")))
							)),fluidRow(
								column(4,
									uiOutput("selMod")
								)
							),fluidRow(
								column(12,
									wellPanel(
										div(id="modPlot", class="vPlot",
											visNetworkOutput('modVis', width="99%", height="800px")
										)
									)
								)
							)
						)
					)))
				),tabPanel(value="radar_chart", title="Module Comparison Chart",
					fluidRow(column(12,
						div(id="radar_info", class="wrap", actionButton("radarChart_info", "Radar Chart Info", icon=icon("info-circle"))),
						hr()
					)),fluidRow(column(12,
						wellPanel(
							div(id="rChart", class="rChart",
								#chartJSRadarOutput('radarChart', width="100%", height="500")
								chartJSRadarOutput('radarChart', width="20%", height="20%")
							)
						)
					))
				),tabPanel(value="enrichment_panel", title="Annotation Enrichment",
					fluidRow(column(12,
                                                tabBox(
                                                        id="ann_tbox", title="", width=12,
                                                        tabPanel(value="enrichment_table", title="Enrichment Table",
                                                                fluidRow(column(12,
                                                                        box(title="Download & Export", status="success", solidHeader=TRUE, collapsible=TRUE, collapsed=TRUE, width=4,
										downloadButton('downloadGO', label='Export Enriched GO Tables')
                                                                        ),
                                                                        div(id="et_info", class="wrap", actionButton("enrichmentDT_info", "Annotation Enrichment Info", icon=icon("info-circle")))
                                                                )),fluidRow(column(4,
                                                                        uiOutput("selModEnrich")
                                                                )),fluidRow(column(12,
                                                                        DT::dataTableOutput('enrichmentDT')
                                                                ))
                                                        ),tabPanel(value="sim_map", title="Module Similarity by GO Heatmap",
                                                                fluidRow(column(12,
									box(title="Download & Export", status="success", solidHeader=TRUE, collapsible=TRUE, collapsed=TRUE, width=4,
										downloadButton('downloadSimMat', label='Download Matrix of GO based Module Similarity')
									),
									div(id="heatmap_info", class="wrap", actionButton("simHeatmap_info", "Similarity Heatmap Info", icon=icon("info-circle")))
								)),fluidRow(column(12,
                                                                        plotOutput('simHeatmap', width="99%", height="800px")
                                                                ))
                                                        )
                                                )
                                        ))
				),tabPanel(value="rank_table", title="Rank Table",
					fluidRow(column(12,
						box(title="Download & Export", status="success", solidHeader=TRUE, collapsible=TRUE, collapsed=TRUE, width=4,
							downloadButton('downloadRankDT', label='Export Ranked Gene Tables')
						),
						div(id="rt_info", class="wrap", actionButton("rankDT_info", "Rank Table Info", icon=icon("info-circle")))
					)),fluidRow(column(12,
						DT::dataTableOutput('rankDT')
					))
				),tabPanel(value="res_mod", title="Response Module Optimization",
                                        fluidRow(column(3,
                                                        uiOutput("selOptMod")
                                                ),column(3,
                                                        uiOutput("selConMod")
                                                ),column(3,
                                                        uiOutput("selIndMod")
                                                )
                                        ),fluidRow(column(3,
                                                        #actionButton("rChartButton", "Plot"),
							shinyBS::bsButton("rChartButton", label="Plot", style="danger", icon=icon("exclamation-circle")),
							hr()
                                                )
					),fluidRow(column(12,
                                                tabBox(
                                                        id="opt_tbox", title="", width=12,
                                                        tabPanel(value="radar_chart_opt", title="Radar Chart",
                                                                fluidRow(column(12,
									box(title="Download & Export", status="success", solidHeader=TRUE, collapsible=TRUE, collapsed=TRUE, width=4,
										downloadButton('downloadOptRankDT', label='Export Ranked Gene Tables (Optimized)')
									),
									div(id="opt_radar_info", class="wrap", actionButton("optRadarChart_info", "Radar Chart Info", icon=icon("info-circle"))),
									hr()
								)),fluidRow(column(12,
                                                                        wellPanel(
                                                                                div(id="rChart", class="rChart",
                                                                                        chartJSRadarOutput('optRadarChart', width="70%", height="70%")
                                                                                )
                                                                        )
                                                                ))
                                                        ),tabPanel(value="tile_plot", title="GO Enrichment Tile Plot",
                                                                fluidRow(column(12,
									box(title="Download & Export", status="success", solidHeader=TRUE, collapsible=TRUE, collapsed=TRUE, width=4,
										downloadButton('downloadTilePlot', label='Export Tileplots')
									),
                                                                        div(id="tile_info", class="wrap", actionButton("tile_plot_panel_info", "Tile Plot Info", icon=icon("info-circle"))),
                                                                        hr(),
                                                                        tabBox(
                                                                               id="tile_plot_panel", title="", width=12,
                                                                               tabPanel(value="BP", title="Biological Process", plotOutput('tilePlotBP', width = "99%", height = "800px")),
                                                                               tabPanel(value="CC", title="Cellular Component", plotOutput('tilePlotCC', width = "99%", height = "800px")),
                                                                               tabPanel(value="MF", title="Molecular Function", plotOutput('tilePlotMF', width = "99%", height = "800px"))
                                                                        )
                                                                ))
                                                        )
                                                )
                                        ))
				)
			)
		#)),fluidRow(column(12,
                #        actionButton("rserve_submit", "START SOCKET DAEMON")
		))
	)
)
