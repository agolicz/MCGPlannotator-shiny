library(shiny)
library(shinythemes)
library(caret)
library(ggplot2)
library(RColorBrewer)
setwd("/home/ubuntu/mcgplannotator_pipeline")

ui <- fluidPage(theme = shinytheme("spacelab"),
                titlePanel("Rice-MCGPlannotator"),
                sidebarLayout(
                  sidebarPanel(
                    selectInput(inputId = "type",label = "Tissue type (ET):", choices = c("rep","veg")),
                    textAreaInput(inputId = "keypos", label = "Key words or phrases to include (semicolon separated - no spaces surrounding semicolons)", value="pollen;reproduction"),
                    textAreaInput(inputId = "keyneg", label = "Key words or phrases to exclude (semicolon separated) - no spaces surrounding semicolons", value="vegetative"),
                    numericInput(inputId = "alpha", label = "Phenotype (P):", value=0.6),
                    numericInput(inputId = "beta", label = "Homology (H):", value=0.6),
                    numericInput(inputId = "gamma", label = "Community participation (CP):", value=0.4),
                    numericInput(inputId = "delta", label = "Community function (CF):", value=0.3),
                    numericInput(inputId = "epsilon", label = "Sequence diversity (D):", value=0.2),
                    numericInput(inputId = "zeta", label = "Expression value (EV):", value=0.1),
                    numericInput(inputId = "possize", label = "Size of positive training set:", value=100),
                    numericInput(inputId = "negsize", label = "Size of negative training set:", value=250),
                    numericInput(inputId = "negpool", label = "Percetage cut-off between positive and negative set", value=5),
                    selectInput(inputId = "codtype", label = "Coding or non-coding", choices = c("coding", "non-coding")),
                    checkboxGroupInput(inputId = "cfeats", label = "Features to be used in the classifier", choices = c("ET", "P", "H", "CP", "CF", "D", "EV"), selected = c("ET", "P", "H", "CP", "CF", "D", "EV")),
                    actionButton(inputId = "start", label = "Run MCGPlannotator"),
                    downloadButton("downloadData", label = "Download results")
                  ),
                  mainPanel(
                    tabsetPanel(type = "tabs",
                                tabPanel("Analysis", verbatimTextOutput("confusion"), plotOutput(outputId = "rocPlot"), plotOutput(outputId = "scramPlot"), tableOutput("table")),
                                tabPanel("Help",
                                         h5("Expression type (ET, value = 0 or 1) - Is the highest expression recorded in chosen tissue type (rep - reproductive or veg - vegetative)."),
                                         h5("Phenotype category (P, value = 0 or 1) - Is the phenotype category consistent with tissue type."),
                                         h5("Sequence homology (H, value = 0 or 1) - Is the homolog annotated with functions found among key words."),
                                         h5("Community participation (CP, value = 0 or 1) - Is the gene found within a community in co-expression network."),
                                         h5("Community function (CF, value = 0 or 1) - Does the gene belong to a community annotated with functions found among the key words."),
                                         h5("Sequence diversity (D,  value = 0 or 1) - Does the gene display low sequence diversity within the species."),
                                         h5("Expression value (EV, value = log(FPKM)) - The FPKM value for the gene in the tissue with highest expression."),
                                         h5("The weights of contribution of P, H, CP, CF, D and EV to the PI score can be adjusted."),
                                         h5("For more information consult: xxx"),
                                         img(src="Fig1.09122017.png"))
                    )
                  )
                )
)

server <- function(input, output) {
  
  v <- reactiveValues(stats = NULL, roc = NULL, data = NULL, out = NULL, sf = NULL, sid = NULL, scamr = NULL)
  
  observeEvent(input$start, {
    cdate <- format(Sys.time(), "%d%b%Y%H%M%S")
    outdir <- paste("./MCGPlannotator_results", cdate, sep="")
    seid <- paste("MCGPlannotator", cdate, sep="")
    print(seid)
    #prep <- paste("mkdir", outdir,sep=" ")
    #print(prep)
    #system(prep)
    pars <- ""
    pars2 <- ""
    if(input$codtype=="coding") {
      pars <- paste(input$alpha,input$beta,input$gamma,input$delta,input$epsilon,input$zeta,sep=",")
      pars2 <- paste("0","0","0","0","0",sep=",")
    }
    else {
      pars <- paste("0","0","0","0","0","0",sep=",")
      pars2 <- paste(input$alpha,input$gamma,input$delta,input$epsilon,input$zeta,sep=",")
    }
    feats <- paste(input$cfeats,collapse=",")
    print(feats)
    comr <- paste("bash", "runMCGPlannotator.sh", paste('"',input$keypos, '"', sep=""), paste('"', input$keyneg, '"', sep=""), paste('"', input$type, '"', sep=""), paste('"', pars, '"', sep=""), paste('"', pars2, '"', sep=""), paste('"', input$codtype, '"', sep=""), input$possize, input$negsize, input$negpool, paste('"', feats, '"', sep=""), paste('"', outdir, '"', sep=""), sep=" ")
    print(comr)
    system(comr)
    
    ts<-scan(paste(outdir, "error.txt", sep="/"), what="")
    v$sf <- ts
    cat(paste(paste(outdir, "error.txt", sep="/"),"\n",sep=""))
    ts2<-scan(paste(outdir, "error.scram.txt", sep="/"), what="")
    
    if(ts2!="SUCCESS") {
      v$scramr <- NULL
    }
    else {
      load(paste(outdir, "classifer.roc.scram", sep = "/"))
      v$scramr <- pdfc
    }
    
    if(ts!="SUCCESS") {
      v$stats <- "Could not build classifier, please choose different parameters."
      v$data <-data.frame()
      v$roc <- NULL
      v$out <- NULL
      v$sid <- seid
    }
    
    else {
      v$sid <- seid
      v$out <- outdir
      load(paste(outdir, "classifer.stats", sep = "/"))
      v$stats<-sim.res
      load(paste(outdir, "classifer.roc", sep = "/"))
      v$roc <- pdfc
      
      t<-data.frame()
      if(input$codtype=="coding") {
        t<-read.csv(paste(outdir, "pos.cod.table", sep="/"),sep="\t",header=F)
        names(t)<-c("GeneID", "ET", "EV", "P" , "H", "CP", "CF", "D", "PI", "GeneID2", "HighestExpression","MutantLinesCount","MutantIDs","MostCommonPheno","MostCommonPhenoCount","MostCommonPhenoCategory","MostCommonPhenoCategoryCount")
        v$data <- t
      } 
      else {
        t<-read.csv(paste(outdir, "pos.nc.table", sep="/"),sep="\t",header=F)
        names(t)<-c("GeneID", "ET", "EV", "P" , "H", "CP", "CF", "D", "PI", "GeneID2", "HighestExpression","MutantLinesCount","MutantIDs","MostCommonPheno","MostCommonPhenoCount","MostCommonPhenoCategory","MostCommonPhenoCategoryCount")
        v$data <- t
      }
    }
  })
  
  output$confusion <- renderPrint({
    #if (is.null(v$data)) return()
    if(is.null(v$stats))  {
      cat(paste("Press \'Run MCGPlannotator\' to begin...","\n",sep=""))
      cat("Results may take several minutes to appear or update...")
    }
    else{
      cat(paste("Session id: ", v$sid, "\n", "\n", sep=""))
      if(v$sf != "SUCCESS"){
        cat(paste(v$stats,"\n",sep=""))
      }
      else {
        tdf.sn <- subset(v$stats, measure == "Sensitivity")
        tdf.sp <- subset(v$stats, measure == "Specificity")
        tdf.acc <- subset(v$stats, measure == "Accuracy")
        tdf.auc <- subset(v$stats, measure == "AUC")
        tdf.mcc <- subset(v$stats, measure == "MCC")
        cat("Sensitivity: ", mean(tdf.sn$value),"(",sd(tdf.sn$value),")","\n",sep="")
        cat("Specificity: ", mean(tdf.sp$value),"(",sd(tdf.sp$value),")","\n",sep="")
        cat("Accuracy: ", mean(tdf.acc$value),"(",sd(tdf.acc$value),")","\n",sep="")
        cat("AUC: ", mean(tdf.auc$value),"(",sd(tdf.auc$value),")","\n",sep="")
        cat("MCC: ", mean(tdf.mcc$value),"(",sd(tdf.mcc$value),")","\n",sep="")
      }
    }
  })
  
  output$rocPlot <- renderPlot({
    if(is.null(v$roc)) {
      return()
      }
    else {
      ggplot(data= v$roc, aes(x=fpr, y=tpr, color=fold, linetype=fold)) + geom_line(size=1, alpha=0.5)+ ggtitle("ROC curve") +ylab("TPR")+xlab("FPR")+ theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ scale_colour_brewer(palette="Set1")+theme(legend.title=element_blank())
    }
  }, width=400, height=300)
  
  output$scramPlot <- renderPlot({
    if(is.null(v$scramr)) {
      return()
    }
    else {
      ggplot(data= v$scramr, aes(x=fpr, y=tpr, color=fold, linetype=fold)) + geom_line(size=1, alpha=0.5) + ggtitle("ROC curve random control - scrambled labels") +ylab("TPR")+xlab("FPR")+ theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ scale_colour_brewer(palette="Set1")+theme(legend.title=element_blank())
    }
  }, width=400, height=300)
  
  output$table <- renderTable({
    #if (is.null(v2$data)) return()
    if(is.null(v$data))  {
      head(v$data)
    }
    else{
      head(v$data,n=nrow(v$data))
    }
  }, rownames = TRUE)
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("output", "zip", sep=".")
    },
    content = function(fname) {
      cat(paste("Zip check: ",v$out, "\n", sep=""))
      fs <- c()
      if(is.null(v$out)) {
        fs <- c(fs, "download.error.txt")
        write("download.error.txt")
        print(fs)
        zip(zipfile=fname, files=fs)
      }
      #tmpdir <- tempdir()
      #setwd(tempdir())
      else {
        if(input$codtype=="coding") {
          for (i in c(paste(v$out,"bayes.cod.classifier.model",sep="/"), paste(v$out,"PI.scores.cod.top.tsv",sep="/"), paste(v$out,"PI.scores.cod.non.tsv",sep="/"), paste(v$out,"prediction.results.cod.tsv",sep="/"), paste(v$out,"results.readme.txt",sep="/"), paste(v$out,"settings.txt",sep="/"))) {
            fs <- c(fs, i)
            write(i)
          }
        }
        else {
          for (i in c(paste(v$out,"bayes.nc.classifier.model",sep="/"), paste(v$out,"PI.scores.nc.top.tsv",sep="/"), paste(v$out,"PI.scores.nc.non.tsv",sep="/"), paste(v$out,"prediction.results.nc.tsv",sep="/"), paste(v$out,"results.readme.txt",sep="/"), paste(v$out,"settings.txt",sep="/"))) {
            fs <- c(fs, i)
            write(i)
          }
        }
        print(fs)
        zip(zipfile=fname, files=fs)
      }
    },
    contentType = "application/zip"
  )
}

shinyApp(ui = ui, server = server)
