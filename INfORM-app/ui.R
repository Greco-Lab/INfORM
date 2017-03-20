library(shiny)
library(shinyjs)
library(colourpicker)
library(shinyBS)
library(shinydashboard)
library(DT)
library(visNetwork)

jsCode <- 'shinyjs.plot_bg_col = function(col){
        var el = $("." + "vPlot");
        el.css("background-color", col);
}'

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
		extendShinyjs(text = jsCode),
		fluidRow(
			box(
				title="Upload", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, width=12,
				fluidRow(
					column(12,
					    h4("Upload Gene Expression Data")
					)
				),
				fluidRow(
					column(4,
						fileInput(inputId="gx", label="Gene Expression Table")
					),
					column(4,
						fileInput(inputId="dgx", label="Differential Gene Expression Table")
					),
					column(4,
						fileInput(inputId="mat", label="Adjacency Matrix (Pre-Inferred with INfORM) [OPTIONAL]")
					)
				),
				fluidRow(
					column(6,
						#checkboxInput("checkUnconnected", "Remove Unconnected Genes", value=FALSE),
						uiOutput("selOrganism"),
						hr(),
						actionButton("runINfORM", "Run INfORM"),
						bsTooltip("runINfORM", "Run INfORM with provided Gene Expression Table and Differential Expression Table!", placement="right")
					),
					column(6,
						uiOutput("selCores")
					)
				),
				fluidRow(
					column(12,
						a(id = "showAdvancedOptions", "Show/Hide Advanced Options"),
						hidden(div(id="advancedOptions",
						    wellPanel(
							fluidRow(
								column(12,
								    h4("Select Parameters to Run MINET (Mutual Information NETworks)"),
								    uiOutput("selMethod"),
								    uiOutput("selEst"),
								    uiOutput("selDisc")#,
								    #uiOutput("selCores")
								)
							)
						    ),
						    fluidRow(
							column(6,
								h4("Configure Gene Ranking"),
								wellPanel(
									uiOutput("rankAttr") 
								)
							),
							column(6,
								h4("Candidate Selection"),
								wellPanel(
									fluidRow(
										column(3, numericInput("candid_val", "Top Candidate Genes", value="5")),
										column(8, radioButtons("candid_type", NULL, choices=c("Count"="count", "Percent"="precentage"), selected="count", inline=TRUE))
									)
								)
							)
						    ),
						    fluidRow(
							column(12,
							    h4("Sub-Graph Selection by Ranked Genes"),
							    wellPanel(
								fluidRow(
								    column(4, align="center", h4("Small Set"),
									    fluidRow(
										    numericInput("l1_val", NULL, value="5", width="20%")
									    )
								    ),
								    column(4, align="center", h4("Medium Set"),
									    fluidRow(
										    numericInput("l2_val", NULL, value="10", width="20%")
									    )
								    ),
								    column(4, align="center", h4("Large Set"),
									    fluidRow(
										    numericInput("l3_val", NULL, value="20", width="20%")
									    )
								    )
								)
							    )
							)
						    )
						))
					)
				)
			)
		),
		fluidRow(
			tabBox(
				id="display", title="Display Area", width=12,
				tabPanel(value="corMat_display", title="Adjacency Matrix", 
					fluidRow(
						valueBoxOutput('totalGeneBox'),
						valueBoxOutput('connGeneBox'),
						valueBoxOutput('unGeneBox')
					),
					fluidRow(
						column(12,
							fluidRow(
								column(4,
									downloadButton('downloadMat', label='Download Adjacency Matrix'),
									hr()
								)
							),
							fluidRow(
								column(4,
									downloadButton('downloadUnGenes', label='Download Un-Connected Genes List')
								)#,
								#column(4,
								#	downloadButton('downloadMat', label='Download Adjacency Matrix')
								#)
							),
							fluidRow(
								tableOutput('unGeneList')
							)
						)
					)
				),
				tabPanel(value="net_display", title="Network",
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
											column(4,colourpicker::colourInput("vColP", "Node Color +ve", value="salmon")),
											column(4,colourpicker::colourInput("vColN", "Node Color -ve", value="lightblue")),
											column(4,colourpicker::colourInput("vColD", "Node Color Default", value="lightgrey"))
										),
										fluidRow(
											column(4,colourpicker::colourInput("vhColP", "Node Highlight Color +ve", value="red")),
											column(4,colourpicker::colourInput("vhColN", "Node Highlight Color -ve", value="royalblue")),
											column(4,colourpicker::colourInput("vhColD", "Node Highlight Color Default", value="grey"))
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
												column(4,colourpicker::colourInput("vColT", "Node Color Top Ranked", value="red")),
												column(4,colourpicker::colourInput("vColB", "Node Color Bottom Ranked", value="yellow"))
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
										column(6,colourpicker::colourInput("eColP", "Edge Color +ve", value="salmon")),
										column(6,colourpicker::colourInput("eColN", "Edge Color -ve", value="lightblue"))
									),
									fluidRow(
										column(6,colourpicker::colourInput("ehColP", "Edge Highlight Color +ve", value="red")),
										column(6,colourpicker::colourInput("ehColN", "Edge Highlight Color -ve", value="royalblue"))#,
									)
								),
								wellPanel(
									h4("Plot Aesthetics"),
									fluidRow(
										column(6,selectInput("vShape", "Node Shape", choices=c("circle","ellipse","database","box","text"), multiple=FALSE, selected="circle")),
										column(6,colourpicker::colourInput("vBrdCol", "Node Border Color", value="black"))
									),
									fluidRow(
										column(6,numericInput("vSize", "Node Label Size", value="20")),
										column(6,colourpicker::colourInput("vLblCol", "Node Label Color", value="#343434"))
									),
									fluidRow(
										column(6,numericInput("eWidth", "Edge Width", value="5", width="100%")),
										column(6,colourpicker::colourInput("bgCol", "Background Color", value="white"))
									),
									fluidRow(
										column(6,numericInput("dDepth", "Depth of Highlighting Nearest Nodes", value="2", width="100%"))
									)
								)
							)
						)
					)),
					tabBox(
						id="net_panel", title="", width=12,
						tabPanel(value="main", title="Main",
							fluidRow(
								box(title="Download Options", status="success", solidHeader=TRUE, collapsible=TRUE, collapsed=TRUE, width=4,
									wellPanel(
										selectInput("graphExportFormat", NULL, choices=c("edgelist", "pajek", "ncol", "lgl", "graphml", "dimacs", "gml", "dot", "leda"), multiple=FALSE, selected="dot"),
										downloadButton('downloadGraph', label='Export Graph')
									),
									downloadButton('downloadGraphGenes', label='Download Graph Gene List')
								),
								div(id="main_info", class="wrap", actionButton("IGraphVis_info", "Main Graph Info", icon=icon("info-circle")))
							    ),
							    fluidRow(
								    column(12,
									    wellPanel(
										    div(id="mainPlot", class="vPlot",
											    visNetworkOutput('IGraphVis', width = "100%", height = "1000px")
										    )
									    )
								    )
							    )
						    ),
						tabPanel(value="small", title="Small Set",
							fluidRow(
								box(title="Download Options", status="success", solidHeader=TRUE, collapsible=TRUE, collapsed=TRUE, width=4,
									wellPanel(
										selectInput("smallSubNetExportFormat", NULL, choices=c("edgelist", "pajek", "ncol", "lgl", "graphml", "dimacs", "gml", "dot", "leda"), multiple=FALSE, selected="dot", width="50%"),
										downloadButton('downloadSmallSubNet', label='Export Sub-Network')
									),
									downloadButton('downloadSmallSubNetGenes', label='Download Sub-Network Gene List')
								),
								div(id="small_info", class="wrap", actionButton("smallSubNetVis_info", "Small Set Graph Info", icon=icon("info-circle")))
							),
							fluidRow(
								column(12,
									wellPanel(
										div(id="sPlot", class="vPlot",
											visNetworkOutput('smallSubNetVis', width = "99%", height = "800px")
										)
									)
								)
							)
						),
						tabPanel(value="medium", title="Medium Set",
							fluidRow(
								box(title="Download Options", status="success", solidHeader=TRUE, collapsible=TRUE, collapsed=TRUE, width=4,
									wellPanel(
										selectInput("mediumSubNetExportFormat", NULL, choices=c("edgelist", "pajek", "ncol", "lgl", "graphml", "dimacs", "gml", "dot", "leda"), multiple=FALSE, selected="dot", width="50%"),
										downloadButton('downloadMediumSubNet', label='Export Sub-Network')
									),
									downloadButton('downloadMediumSubNetGenes', label='Download Sub-Network Gene List')
								),
								div(id="medium_info", class="wrap", actionButton("mediumSubNetVis_info", "Medium Set Graph Info", icon=icon("info-circle")))
							),
							fluidRow(
								column(12,
									wellPanel(
										div(id="mPlot", class="vPlot",
											visNetworkOutput('mediumSubNetVis', width = "99%", height = "800px")
										)
									)
								)
							)
						),
						tabPanel(value="large", title="Large Set",
							fluidRow(
								box(title="Download Options", status="success", solidHeader=TRUE, collapsible=TRUE, collapsed=TRUE, width=4,
									wellPanel(
										selectInput("largeSubNetExportFormat", NULL, choices=c("edgelist", "pajek", "ncol", "lgl", "graphml", "dimacs", "gml", "dot", "leda"), multiple=FALSE, selected="dot", width="50%"),
										downloadButton('downloadLargeSubNet', label='Export Sub-Network')
									),
									downloadButton('downloadLargeSubNetGenes', label='Download Sub-Network Gene List')
								),
								div(id="large_info", class="wrap", actionButton("largeSubNetVis_info", "Large Set Graph Info", icon=icon("info-circle")))
							),
							fluidRow(
								column(12,
									wellPanel(
										div(id="lPlot", class="vPlot",
											visNetworkOutput('largeSubNetVis', width = "99%", height = "800px")
										)
									)
								)
							)
						)
					)))
				),
				tabPanel(value="enrichment_table", title="Annotation Enrichment", 
					fluidRow(column(12,
						box(title="Download & Export", status="success", solidHeader=TRUE, collapsible=TRUE, collapsed=TRUE, width=4,
							wellPanel(
								downloadButton('downloadGO', label='Download Progressively Enriched GO List')
							)
						),
						div(id="et_info", class="wrap", actionButton("enrichmentDT_info", "Annotation Enrichment Info", icon=icon("info-circle")))
					)),
					fluidRow(column(12,
						DT::dataTableOutput('enrichmentDT')
					))
				),
				tabPanel(value="tile_plot", title="GO Enrichment Tile Plot", 
					fluidRow(column(12,
					div(id="tile_info", class="wrap", actionButton("tile_plot_panel_info", "Tile Plot Info", icon=icon("info-circle"))),
					hr(),
					tabBox(
						id="tile_plot_panel", title="", width=12,
						tabPanel(value="BP", title="Biological Process", plotOutput('tilePlotBP', width = "99%", height = "800px")),
						tabPanel(value="CC", title="Cellular Component", plotOutput('tilePlotCC', width = "99%", height = "800px")),
						tabPanel(value="MF", title="Molecular Function", plotOutput('tilePlotMF', width = "99%", height = "800px"))
					)))
				),
				tabPanel(value="rank_table", title="Rank Table", 
					div(id="rt_info", class="wrap", actionButton("rankDT_info", "Rank Table Info", icon=icon("info-circle"))),
					hr(),
					downloadButton('downloadRankDT', label='Download Ranked Gene Table'),
					DT::dataTableOutput('rankDT')
				)
			)
		)
	)
)

