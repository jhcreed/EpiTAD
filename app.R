#######################################################################################
# working directory, options, libraries

setwd("/Volumes/Lab_Gerke/ShinyApps/LD_HIC")

options(stringsAsFactors=FALSE)

library(shiny)
library(shinydashboard)
library(readr)
# library(sqldf)
library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(labeling)
library(ggrepel)
library(shinythemes)
library(httr)
library(ggbio)
library(GenomicRanges)
library(biomaRt) #needed for genes until plotting updated

#######################################################################################
# set up dashboard

ui <- dashboardPage(skin = "black",
                    dashboardHeader(title ="HIC viewer"),
                    dashboardSidebar(
                      sidebarMenu(id="tabs", 
                                  menuItem("SNPs",tabName="tab1"), menuItem("Info",tabName="tab2"))
                    ),
                    dashboardBody(
                      tabItems(
                        tabItem(tabName="tab1",
                                fluidRow(
                                  box(title = "Inputs",
                                      sliderInput("chr","Chromosome",min=1,max=22,value=22),
                                      #numericInput
                                      numericInput("CHRboundary1","Start Position",value=100000),
                                      numericInput("CHRboundary2","End Position",value=200000),
                                      #sliderInput("CHRboundary","Boundaries",min=0,max=247200000,value=c(0,20000)),
                                      textInput("rsids","Specific SNP of Interest?"),
                                      sliderInput("rthresh","R2 Threshold",min=0.1,max=1,value=0.8),
                                      actionButton("go1","Go!")),
                                  tabBox(title="Tables",
                                         tabPanel("HI-C",
                                                  textOutput("HICboundaries"),
                                                  br(),
                                                  tableOutput("HICGenes")),
                                         tabPanel("LD", dataTableOutput("LDresults"))
                                  )
                                ),
                                fluidRow(
                                  column(12,align="center",offset=2,
                                         box(title="Visuals",
                                             sliderInput("HICintensity","HIC intensity",min=1,max=50,value=4),
                                             plotOutput("HICheatmap", height = "600px", width = "750px"),width=8))
                                )
                        ),
                        tabItem(tabName="tab2", textOutput("LDinfo"))
                      )
                    )
)


server <- shinyServer(function(input, output) {
  
  ###do we need both of these datasets and the biomart info we load?
  ####need tads as it gives the TAD boundaries defined by dixon but tad_genes made obsolete by geneannot
  tads <- read.table("./Data/total.combined.domain",sep="")
  ####DELETE BELOW ONCE EVERYTHING IS MOVED:
  tad_genes <- read.table("./Data/tad_genes.txt", sep="\t", header=TRUE)
  ensembl54=useMart("ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl")
  dbsnp = useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
  
  geneannot <- read.table("Data/geneAnnot.txt", sep="\t", header=TRUE)
  ###we should not pre-load this file! rather, query rows on an as-needed basis
  #snpannot <- read.table("Data/snpAnnot.txt", sep="\t", header=TRUE)
  ###here's an example using sqldf
  # mysnp <- "rs765616855"
  # sqltext <- paste0("select*from file where refsnp_id='", mysnp, "'")
  # tempannot <- read.csv.sql("Data/snpAnnot.txt", sep="\t", header=TRUE,
  #                           sql=sqltext)
  ###here's an example using readr which should be faster since uses chunks
  mysnp <- "rs765616855"
  ####callback from DataFrameCallback$new(function(x, pos) x['refsnp_id'==mysnp,]) to DataFrameCallback$new(function(x, pos) x[x$refsnp_id==mysnp,])
  ####made to work in the context of shiny
  # tempannot <- read_delim_chunked(paste0("Data/snpAnnot",input$chr,".txt"), delim="\t",
  #                                 callback=DataFrameCallback$new(function(x, pos) x[x$refsnp_id==input$rsids,]),
  #                                 progress=FALSE, chunk_size=50000)
  ### I don't think we need these functions. Let's just perform the subset inline with 
  ### something like snp[loc1<=x & x<=loc2, 2,3]
  getBoundaries <- function(x, data) {
    tmp <- data %>%
      filter(V2 <= x, x <= V3)
    return(tmp[,2:3])
  }
  
  getGenes <- function(x, data) {
    tmp <- data %>%
      filter(V3 <= x, x <= V4)
    return(tmp$V5)
  }
  
  ####seeing if read_delim_chunked works better than R data frames then will return to the dropbox issue
    # observeEvent(input$go1,{
  #   if (input$chr==1){
  #     response<-GET("https://www.dropbox.com/s/4z1cga3mmijlzg5/chr1.RData?dl=1")
  #     writeBin(response$content,"chr_dat.RData")
  #     load("chr_dat.RData")
  #   }
  #   if (input$chr==2){
  #     response<-GET("https://www.dropbox.com/s/2p59kghyjds3j8d/chr2.RData?dl=1")
  #     writeBin(response$content,"chr_dat.RData")
  #     load("chr_dat.RData")
  #   }
  #   if (input$chr==3){
  #     response<-GET("https://www.dropbox.com/s/j8vgztz6ronqmy0/chr3.RData?dl=1")
  #     writeBin(response$content,"chr_dat.RData")
  #     load("chr_dat.RData")
  #   }
  #   if (input$chr==4){
  #     response<-GET("https://www.dropbox.com/s/bgwuutleldol6x2/chr4.RData?dl=1")
  #     writeBin(response$content,"chr_dat.RData")
  #     load("chr_dat.RData")
  #   }
  #   if (input$chr==5){
  #     response<-GET("https://www.dropbox.com/s/dix4998uk1yy6m6/chr5.RData?dl=1")
  #     writeBin(response$content,"chr_dat.RData")
  #     load("chr_dat.RData")
  #   }
  #   if (input$chr==6){
  #     response<-GET("https://www.dropbox.com/s/lmw98zkoolypkt6/chr6.RData?dl=1")
  #     writeBin(response$content,"chr_dat.RData")
  #     load("chr_dat.RData")
  #   }
  #   if (input$chr==7){
  #     response<-GET("https://www.dropbox.com/s/n6wsj3gdq05azl0/chr7.RData?dl=1")
  #     writeBin(response$content,"chr_dat.RData")
  #     load("chr_dat.RData")
  #   }
  #   if (input$chr==8){
  #     response<-GET("https://www.dropbox.com/s/8ryjzgo4ife2byj/chr8.RData?dl=1")
  #     writeBin(response$content,"chr_dat.RData")
  #     load("chr_dat.RData")
  #   }
  #   if (input$chr==9){
  #     response<-GET("https://www.dropbox.com/s/w8o3s4spwsdqcw9/chr9.RData?dl=1")
  #     writeBin(response$content,"chr_dat.RData")
  #     load("chr_dat.RData")
  #   }
  #   if (input$chr==10){
  #     response<-GET("https://www.dropbox.com/s/6o434hn5sc7kibb/chr10.RData?dl=1")
  #     writeBin(response$content,"chr_dat.RData")
  #     load("chr_dat.RData")
  #   }
  # })
  
  data.hic.genes<-eventReactive(input$go1,{
    region<-paste0(input$chr,":",input$CHRboundary1,":",input$CHRboundary2)
    results<-getBM(attributes = c("hgnc_symbol","start_position","end_position"),
                   filters=c("chromosomal_region"), values=region,mart = ensembl54)
    results<-results[results$hgnc_symbol!="",]
    for(i in 1:nrow(results)){results$y[i]<-i}
    return(results)
  })
  
  snp.position<-eventReactive(input$go1,{
    pos1<-getBM(attributes = c("refsnp_id","chr_name","chrom_start"),
                filters=c("snp_filter"), values=input$rsids,mart = dbsnp)
    return(pos1)
  })
  
  ld.data2<-eventReactive(input$go1,{
    ####QUICK UPLOAD
    #load(paste0("./Data/ld_chr",input$chr,".RData"))
    # data<-ld.data
    # data<-ld.data[ld.data$V4==input$rsids | 
    #                 (ld.data$V1>=input$CHRboundary[1] & ld.data$V1<=input$CHRboundary[2] & 
    #                    ld.data$V2>=input$CHRboundary[1] & ld.data$V2>=input$CHRboundary[1]),]
    #data<-data[V7>=input$rthresh,]
    ldFun<-function(x,pos) x[x$V4==input$rsids | 
                               (x$V1>=input$CHRboundary1 & x$V1<=input$CHRboundary2 & 
                                  x$V2>=input$CHRboundary1 & x$V2>=input$CHRboundary2),]
    ld.matrix<-read_delim_chunked("https://www.dropbox.com/s/2pxywhotssditcb/ld_chr10_small.txt?dl=1", delim="\t",
                                  callback=DataFrameCallback$new(ldFun),
                                  progress=FALSE, chunk_size=50000)
    return(ld.matrix)
  })
  
  # #using plyr instead?
  # ld.data3<-eventReactive(input$go1,{
  #   #QUICK UPLOAD
  #   #load(paste0("./Data/ldvect_chr",input$chr,".RData"))
  #   # data<-ld
  #   # data<-ld[ld$V1>=input$CHRboundary[1] & ld$V1<=input$CHRboundary[2] ,]
  #   ldFun2<-function(x,pos) x[x$V1>=input$CHRboundary[1] & x$V1<=input$CHRboundary[2] ,]
  #   ld.vect<-read_delim_chunked(paste0("./Data2/ldvect_chr",input$chr,".txt"), delim="\t",
  #                                 callback=DataFrameCallback$new(ldFun2),
  #                                 progress=FALSE, chunk_size=50000)
  #   return(ld.vect)
  # })
   
  hic.data<-eventReactive(input$go1,{
    ####QUICK UPLOAD
    #load(paste0("./Data/chr",input$chr,"_hic.RData"))
    # data_hic<-as.data.frame(raw_hic)
    # data_hic<-raw_hic[(raw_hic$loc1>=input$CHRboundary[1] & raw_hic$loc1<=input$CHRboundary[2] &
    #                      raw_hic$loc2>=input$CHRboundary[1] & raw_hic$loc2<=input$CHRboundary[2]), ]
    hicFun<-function(x,pos) x[(x$loc1>=input$CHRboundary1 & x$loc1<=input$CHRboundary2 &
                              x$loc2>=input$CHRboundary1 & x$loc2<=input$CHRboundary2), ]
    hic.matrix<-read_delim_chunked("https://www.dropbox.com/s/vlov9d3pb3i1lhv/chr10_hic.txt?dl=1", delim="\t",
                                callback=DataFrameCallback$new(hicFun),
                                progress=FALSE, chunk_size=50000)
    return(hic.matrix)
  })
  
  ctcf.data<-eventReactive(input$go1,{
    load("./Data/CFCTsites.RData")
    dat<-ctcf_loc[ctcf_loc$chr==input$chr & 
                    (ctcf_loc$start>=input$CHRboundary1 & ctcf_loc$start<=input$CHRboundary2 &
                       ctcf_loc$end>=input$CHRboundary1 & ctcf_loc$end<=input$CHRboundary2),]
    return(dat)
  })
  
  output$LDresults<-renderDataTable({
    ld.data3<-ld.data2()
    ld.data3<-ld.data3[ld.data3$V7>=input$rthresh,]
    colnames(ld.data3)<-c("Position 1","Position 2","RSid 1", "RSid 2","r^2")
    return(ld.data3)
  })
  
  # output$LDheatmap<-renderPlot({
  #   ld.heat<-ld.data2()
  #   ld.heat<-ld.heat[order(ld.heat$V1,ld.heat$V2),]
  #   ggplot(ld.heat,aes(V1,V2))+geom_tile(aes(fill=V7),colour="white")+scale_fill_gradient(low="white",high="red")
  # })
  
  output$HICboundaries<-renderText({
    snp<-ld.data2()
    snp<-snp[snp$V4==input$rsids,]
    snp<-snp[1,1]
    tad_data<-tads
    tad_data<-tad_data[tad_data$V1==paste0("chr",input$chr),]
    #getBoundaries(snp,data=tad_data)
    tad_b<-getBoundaries(snp,data=tad_data)
    if (is.na(tad_b[1,1])){
      return(paste0("SNP not located within a known TAD"))
    }else{
      return(paste0("In TAD: boundary from ",tad_b[1,1]," to ",tad_b[1,2]))
    }
  })
  
  output$HICGenes<-renderTable({
    snp<-ld.data2()
    snp<-snp[snp$V4==input$rsids,]
    snp<-snp[1,1]
    tad_data<-tads
    tad_data<-tad_data[tad_data$V1==paste0("chr",input$chr),]
    boundaries<-getBoundaries(snp,data=tad_data)
    location<-paste0(input$chr,":",boundaries[1,1],":",boundaries[1,2])
    genes<-tad_genes
    genes<-genes[genes$region %in% location,]
    return(genes[,1:3])
    
  })
  
  output$HICheatmap<-renderPlot({
    hic<-hic.data()
    #need to call plyr specifically within shiny
    hic2<-plyr::ddply(hic, .(loc1), summarize, contact_value=sum(raw))
    genes<-data.hic.genes()
    ctcf<-ctcf.data()
    #tesing out the quickness of using plyr instead of loading ld vectors 
    #ld<-ld.data3()
    ld<-ld.data2()
    ld<-plyr::ddply(ld, .(V1), summarize, contact_value=sum(V7))
    snp.pos<-snp.position()
    hic$raw_max<-ifelse(hic$raw>=input$HICintensity,input$HICintensity,hic$raw)
    chr_view<-Ideogram(genome="hg19",subchr = paste0("chr",input$chr))# + xlim(GRanges(paste0("chr",input$chr), IRanges(input$CHRboundary[1],input$CHRboundary[2])))
    
    heat_view<-ggplot(hic, aes(loc1, loc2)) + geom_tile(aes(fill = raw_max),colour="white") +
      scale_fill_gradient(low = "white",high = "purple")+ theme(legend.position = "none") +labs(y="loci") +
      geom_vline(xintercept = snp.pos[1,3])
    
    gene_view<-ggplot(genes)+geom_segment(data=genes,aes(x=start_position,y=1,xend=end_position,yend=1,color=hgnc_symbol,size=3,alpha=0.1))+
      ylim(0,2)+geom_label_repel(data=genes, mapping=aes(x=start_position, y=1, label=hgnc_symbol,color=hgnc_symbol), size=3)+
      theme(legend.position = "none",axis.text.y=element_blank(),axis.title.y=element_blank())+
      geom_vline(xintercept = snp.pos[1,3])
    
    if(nrow(ctcf)>=1){
      ctcf_view<-ggplot(ctcf)+geom_segment(data=ctcf,aes(x=start,y=1,xend=end,yend=1),size=3)+ylim(0.5,1.5)+ 
        theme(legend.position = "none",axis.text.y=element_blank(),axis.title.y=element_blank())+
        geom_vline(xintercept = snp.pos[1,3])
    }else{ctcf_view<-ggplot(ctcf)+ geom_blank()}
    
    #editing ld for plyr version
    ld_view<-ggplot(ld)+geom_smooth(data=ld,aes(x=V1,y=contact_value),span=0.2)+
      geom_vline(xintercept = snp.pos[1,3])
    
    hc_view<-ggplot(hic2)+geom_smooth(data=hic2,aes(x=loc1,y=contact_value),span=0.2,color="purple")+
      geom_vline(xintercept = snp.pos[1,3])
    
    tks1<-tracks(CHR=chr_view, HIC=heat_view, Genes=gene_view, CTCF=ctcf_view,LD=ld_view,HIC2=hc_view,heights = c(2,3,1,1,2,2))
    
    return(tks1)
  })
  
})

# Run the application 
shinyApp(ui = ui, server = server)