
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#App created for visualisation of SC Live-cell imaging data. It takes as output the files from
#stitching script

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("/tungstenfs/scratch/ggiorget/zhan/2019/github/Imaging_tania_project/functions_v0.R")
library(shiny)
library(ggplot2)

##global options
timeColumn = 1 # corresponds to column with time position information
tactual = 0.5 #minutes
sigma_threshold = 4 #how many standard deviation away from the mean of the derivative of distance with respect to time to exclude because of wrong stitching
dir = "/tungstenfs/scratch/ggiorget/_LIVECELL/Analysis_Data/PGK_G8_A11_B3_F6_SingleCells/results_analysis_corrected_traj" 
dir_uncorrected = "/tungstenfs/scratch/ggiorget/_LIVECELL/Analysis_Data/PGK_G8_A11_B3_F6_SingleCells/results_analysis_uncorrected_traj/" 
#end

files = list.files(paste0(dir,"/csv/"),pattern = "_reconstructed_tracks.csv")
files_uncorrected = list.files(paste0(dir_uncorrected,"/csv/"),pattern = "_reconstructed_tracks.csv")


##filtering, deleting too sharping increase and decrease  (MUST BE MOFIDIED IF NUMBER OF CHANNELS CHANGES)
#d = 2 column matrix (1 column denominator, second column numerator of derivation), threshold on sd to filter
#output the time point you have to filter
derivative_filter = function(d,threshold){
  derivative = diff(d[,2])/diff(d[,1])
  plot(derivative)
  abline(h=c((mean(derivative))+threshold*sd(derivative),(mean(derivative))-threshold*sd(derivative)))
  return(d[(which(abs(derivative)>(mean(derivative))+threshold*sd(derivative))+1),1])
}



# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Visualisation of Live Cell Imaging data"),
  helpText(p("Analysis pipeline"),
           HTML("<ul><li>Images are cropped and trajectories are constructed using TrackMate</li>"),
           HTML("<li>Trajectories are stitched together using costum script</li>"),
           HTML("<li>Trajectories are stitched together using costum script</li>"),
           HTML("<li>Filtering on pair-wise distance is based on speed of movement</li>"),
           HTML("<li>Pairs: upper=ccf; diagonal = acf; low=scatter</li><ul>")),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      selectInput("d_filter", "Filter (movement)",  choices=c("ON","OFF")),
      selectInput("type", "Which Plot",  choices=c("Pair-wise_dist","ECDF_all_data","Gyration_radius","Auto_cross_pairs")),
      selectInput("cell", "Which cell:",  choices=files),
      selectInput("cell_u", "Which cell (Uncorrected):",  choices=files_uncorrected),
      selectInput("dimension","2 or 3D", choices = c("3D","2D")),
      selectInput("ymax", "Max value in y axis:",  choices=c(1.5,(1:10)/10,1,2,2.5,3,4,5)),
      downloadButton("downloadData", "Reconstructed trajectory"),
      downloadButton("downloadPlot", "Plot"),
      width = 15
    ),
    
    # Show a plot of the generated distribution
    mainPanel( fluidRow(
      column(2, align="right",
             plotOutput(outputId = "plotgraph1", width  = "500px",height = "400px"),
             plotOutput(outputId = "plotgraph2", width  = "500px",height = "400px")
      ))
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  
  pt1 <- reactive({
    if (!(input$type=="Pair-wise_dist")) return(NULL)
    if(input$dimension == "3D"){ posColumns = 2:4}else{posColumns = 2:3 } #columns correpsonding to the x,y,z position
    
    #read data
    reconstructed_tracks = read.csv(paste0(dir,"/csv/",input$cell))
    channels = names(table(reconstructed_tracks$channel))
    
    ####### separate tracks by channels  (MUST BE MOFIDIED IF NUMBER OF CHANNELS CHANGES)
    channel1 = reconstructed_tracks[reconstructed_tracks$channel==channels[1],c(timeColumn,posColumns)]
    channel2 = reconstructed_tracks[reconstructed_tracks$channel==channels[2],c(timeColumn,posColumns)]
    channel3 = reconstructed_tracks[reconstructed_tracks$channel==channels[3],c(timeColumn,posColumns)]
    
    ##distances between channels  (MUST BE MOFIDIED IF NUMBER OF CHANNELS CHANGES)
    dist_c12 = dist_channels(channel1,channel2)
    dist_c23 = dist_channels(channel2,channel3)
    dist_c13 = dist_channels(channel1,channel3)
    dist_c12$t = dist_c12$t*tactual
    dist_c23$t = dist_c23$t*tactual
    dist_c13$t = dist_c13$t*tactual
    
    
    #filter
    if(input$d_filter=="ON"){
      filter = NULL
      filter = c(filter,derivative_filter(dist_c12,sigma_threshold))
      filter = c(filter,derivative_filter(dist_c23,sigma_threshold))
      filter = c(filter,derivative_filter(dist_c13,sigma_threshold))
      filter = unique(filter)
      
      dist_c12 = dist_c12[!(dist_c12$t %in% filter), ]
      dist_c23 = dist_c23[!(dist_c23$t %in% filter), ]
      dist_c13 = dist_c13[!(dist_c13$t %in% filter), ]
    }
    plot(dist_c12,type="l",ylim=c(0,as.numeric(input$ymax)), ylab="Distance in µm", xlab="t in min", main=paste0(input$dimension, " Pair-wise_distances filter ",input$d_filter), col=3)
    lines(dist_c23,type="l",col=4)
    lines(dist_c13,type="l",col=2)
    legend(x="topright",legend=c("Chic vs Tsix","Tsix vs Linx","Chic vs Linx"),col=c(3,4,2),lty=c(1,1))
    
  })
  
  pt2 <- reactive({
    if (!(input$type=="ECDF_all_data")) return(NULL)
    
    if(input$d_filter=="ON"){
      files_distance12 = list.files(paste0(dir,"/csv/"),pattern = "12distances_filtered.csv",full.names = T)
      files_distance23 = list.files(paste0(dir,"/csv/"),pattern = "23distances_filtered.csv",full.names = T)
      files_distance13 = list.files(paste0(dir,"/csv/"),pattern = "13distances_filtered.csv",full.names = T)
    }else{
      files_distance12 = list.files(paste0(dir,"/csv/"),pattern = "12distances.csv",full.names = T)
      files_distance23 = list.files(paste0(dir,"/csv/"),pattern = "23distances.csv",full.names = T)
      files_distance13 = list.files(paste0(dir,"/csv/"),pattern = "13distances.csv",full.names = T)
    }
    
    
    data12 = read.csv(files_distance12[1])
    for(i in 2:length(files_distance12)){
      a = read.csv(files_distance12[i])
      data12 = rbind(data12,a)
    }
    
    data23 = read.csv(files_distance23[1])
    for(i in 2:length(files_distance23)){
      a = read.csv(files_distance23[i])
      data23 = rbind(data23,a)
    }
    data13 = read.csv(files_distance13[1])
    for(i in 2:length(files_distance13)){
      a = read.csv(files_distance13[i])
      data13 = rbind(data13,a)
    }
    
    
    data12$type = "Chic vs Tsix"
    data23$type = "Tsix vs Linx"
    data13$type = "Chic vs Linx"
    
    data = rbind(data12,data23,data13)
    
    ggplot(data,aes(dist,colour = type))+
      stat_ecdf()+
      xlim(0,as.numeric(input$ymax)) +
      labs("Only in 3D by now")
  })
  
  
  
  pt3 <- reactive({
    if (!(input$type=="Gyration_radius")) return(NULL)
    
    if(input$dimension == "3D"){ posColumns = 2:4}else{posColumns = 2:3 } #columns correpsonding to the x,y,z position
    
    reconstructed_tracks = read.csv(paste0(dir,"/csv/",input$cell))
    channels = names(table(reconstructed_tracks$channel))
    
    ####### separate tracks by channels  (MUST BE MOFIDIED IF NUMBER OF CHANNELS CHANGES)
    channel1 = reconstructed_tracks[reconstructed_tracks$channel==channels[1],c(timeColumn,posColumns)]
    channel2 = reconstructed_tracks[reconstructed_tracks$channel==channels[2],c(timeColumn,posColumns)]
    channel3 = reconstructed_tracks[reconstructed_tracks$channel==channels[3],c(timeColumn,posColumns)]
    
    if(input$d_filter=="ON"){
      ##distances between channels  (MUST BE MOFIDIED IF NUMBER OF CHANNELS CHANGES)
      dist_c12 = dist_channels(channel1,channel2)
      dist_c23 = dist_channels(channel2,channel3)
      dist_c13 = dist_channels(channel1,channel3)
      filter = NULL
      filter = c(filter,derivative_filter(dist_c12,sigma_threshold))
      filter = c(filter,derivative_filter(dist_c23,sigma_threshold))
      filter = c(filter,derivative_filter(dist_c13,sigma_threshold))
      filter = unique(filter)
      reconstructed_tracks[!(reconstructed_tracks$t %in% filter),]
      channel1 = reconstructed_tracks[reconstructed_tracks$channel==channels[1],c(timeColumn,posColumns)]
      channel2 = reconstructed_tracks[reconstructed_tracks$channel==channels[2],c(timeColumn,posColumns)]
      channel3 = reconstructed_tracks[reconstructed_tracks$channel==channels[3],c(timeColumn,posColumns)]
    }
    
    #calculate Gyration radius (mean((rk-rcm)^2))
    alldata = merge(merge(channel1,channel2,by=1),channel3,by=1)
    t = alldata[,1] * tactual
    alldata = alldata[,-1]
    
    if(input$dimension == "3D"){
      allx = alldata[,c(1,4,7)]
      ally = alldata[,c(1,4,7)+1]
      allz = alldata[,c(1,4,7)+2]
      xmean = rowMeans(allx)
      ymean = rowMeans(ally)
      zmean = rowMeans(allz)
      
      rg2 = rowMeans((allx-xmean)**2 + (ally-ymean)**2 + (allz-zmean)**2)
    }else{
      allx = alldata[,c(1,3,5)]
      ally = alldata[,c(1,3,5)+1]
      xmean = rowMeans(allx)
      ymean = rowMeans(ally)
      rg2 = rowMeans((allx-xmean)**2 + (ally-ymean)**2)
    }
    rg = sqrt(rg2)
    
    plot(t,rg2,xlab="time (minutes)",ylab = "Gyration radius (um)", type = "l", main = paste0("filter ",input$d_filter))
    
  })
  
  
  pt4 <- reactive({
    if (!(input$type=="Auto_cross_pairs")) return(NULL)
    if(input$dimension == "3D"){ posColumns = 2:4}else{posColumns = 2:3 } #columns correpsonding to the x,y,z position
    
    #read data
    reconstructed_tracks = read.csv(paste0(dir,"/csv/",input$cell))
    channels = names(table(reconstructed_tracks$channel))
    
    ####### separate tracks by channels  (MUST BE MOFIDIED IF NUMBER OF CHANNELS CHANGES)
    channel1 = reconstructed_tracks[reconstructed_tracks$channel==channels[1],c(timeColumn,posColumns)]
    channel2 = reconstructed_tracks[reconstructed_tracks$channel==channels[2],c(timeColumn,posColumns)]
    channel3 = reconstructed_tracks[reconstructed_tracks$channel==channels[3],c(timeColumn,posColumns)]
    
    ##distances between channels  (MUST BE MOFIDIED IF NUMBER OF CHANNELS CHANGES)
    dist_c12 = dist_channels(channel1,channel2)
    dist_c23 = dist_channels(channel2,channel3)
    dist_c13 = dist_channels(channel1,channel3)
    dist_c12$t = dist_c12$t*tactual
    dist_c23$t = dist_c23$t*tactual
    dist_c13$t = dist_c13$t*tactual
    
    
    #filter
    if(input$d_filter=="ON"){
      filter = NULL
      filter = c(filter,derivative_filter(dist_c12,sigma_threshold))
      filter = c(filter,derivative_filter(dist_c23,sigma_threshold))
      filter = c(filter,derivative_filter(dist_c13,sigma_threshold))
      filter = unique(filter)
      
      dist_c12 = dist_c12[!(dist_c12$t %in% filter), ]
      dist_c23 = dist_c23[!(dist_c23$t %in% filter), ]
      dist_c13 = dist_c13[!(dist_c13$t %in% filter), ]
    }
    
    a = merge(merge(dist_c12,dist_c13,by=1),dist_c23,by=1)
    time = a[,1]*tactual
    a = a[,-1]
    colnames(a)= c("Chic vs Tsix","Chic vs Linx","Tsix vs Linx")
    
    #auto correlation
    diagonal.acf = function(x){
      par(new = TRUE)
      auto = acf(x,plot=FALSE)
      plot(time[1:length(auto[[1]])]-min(time),auto[[1]],xlab="lag (Minutes)",ylim=c(0,1),ylab="Auto correlation",type="l")
    }
    
    #cross correlation
    upper.ccf = function(x,y){
      par(new = TRUE)
      auto = ccf(x,y,plot=FALSE)
      plot(c(sort(-time[1:(as.integer(length(auto[[1]])/2)+1)]),
             time[(1:as.integer(length(auto[[1]])/2))+1]-min(time)),
           auto[[1]],xlab="lag (Minutes)",ylab="Cross correlation",type="l",ylim=c(0,1))
    }
    
    panel.smooth <- function(x,y) {par(new=TRUE);
      smoothScatter(x,y, nrpoints=0);
      abline(a=0,b=1,col=2)}
    
    pairs(a,diag.panel = diagonal.acf,
          upper.panel = upper.ccf,
          lower.panel = panel.smooth)  
    
    
  })
  
  
  
  
  
  
  
  
  
  
  
  
  pt1_u <- reactive({
    if (!(input$type=="Pair-wise_dist")) return(NULL)
    if(input$dimension == "3D"){ posColumns = 2:4}else{posColumns = 2:3 } #columns correpsonding to the x,y,z position
    
    #read data
    reconstructed_tracks = read.csv(paste0(dir_uncorrected,"/csv/",input$cell_u))
    channels = names(table(reconstructed_tracks$channel))
    
    ####### separate tracks by channels  (MUST BE MOFIDIED IF NUMBER OF CHANNELS CHANGES)
    channel1 = reconstructed_tracks[reconstructed_tracks$channel==channels[1],c(timeColumn,posColumns)]
    channel2 = reconstructed_tracks[reconstructed_tracks$channel==channels[2],c(timeColumn,posColumns)]
    channel3 = reconstructed_tracks[reconstructed_tracks$channel==channels[3],c(timeColumn,posColumns)]
    
    ##distances between channels  (MUST BE MOFIDIED IF NUMBER OF CHANNELS CHANGES)
    dist_c12 = dist_channels(channel1,channel2)
    dist_c23 = dist_channels(channel2,channel3)
    dist_c13 = dist_channels(channel1,channel3)
    dist_c12$t = dist_c12$t*tactual
    dist_c23$t = dist_c23$t*tactual
    dist_c13$t = dist_c13$t*tactual
    
    
    #filter
    if(input$d_filter=="ON"){
      filter = NULL
      filter = c(filter,derivative_filter(dist_c12,sigma_threshold))
      filter = c(filter,derivative_filter(dist_c23,sigma_threshold))
      filter = c(filter,derivative_filter(dist_c13,sigma_threshold))
      filter = unique(filter)
      
      dist_c12 = dist_c12[!(dist_c12$t %in% filter), ]
      dist_c23 = dist_c23[!(dist_c23$t %in% filter), ]
      dist_c13 = dist_c13[!(dist_c13$t %in% filter), ]
    }
    plot(dist_c12,type="l",ylim=c(0,as.numeric(input$ymax)), ylab="Distance in µm", xlab="t in min", main=paste0(input$dimension, " Pair-wise_distances filter ",input$d_filter," uncorrected"), col=3)
    lines(dist_c23,type="l",col=4)
    lines(dist_c13,type="l",col=2)
    legend(x="topright",legend=c("Chic vs Tsix","Tsix vs Linx","Chic vs Linx"),col=c(3,4,2),lty=c(1,1))
    
  })
  
  pt2_u <- reactive({
    if (!(input$type=="ECDF_all_data")) return(NULL)
    
    if(input$d_filter=="ON"){
      files_distance12 = list.files(paste0(dir_uncorrected,"/csv/"),pattern = "12distances_filtered.csv",full.names = T)
      files_distance23 = list.files(paste0(dir_uncorrected,"/csv/"),pattern = "23distances_filtered.csv",full.names = T)
      files_distance13 = list.files(paste0(dir_uncorrected,"/csv/"),pattern = "13distances_filtered.csv",full.names = T)
    }else{
      files_distance12 = list.files(paste0(dir_uncorrected,"/csv/"),pattern = "12distances.csv",full.names = T)
      files_distance23 = list.files(paste0(dir_uncorrected,"/csv/"),pattern = "23distances.csv",full.names = T)
      files_distance13 = list.files(paste0(dir_uncorrected,"/csv/"),pattern = "13distances.csv",full.names = T)
    }
    
    
    data12 = read.csv(files_distance12[1])
    for(i in 2:length(files_distance12)){
      a = read.csv(files_distance12[i])
      data12 = rbind(data12,a)
    }
    
    data23 = read.csv(files_distance23[1])
    for(i in 2:length(files_distance23)){
      a = read.csv(files_distance23[i])
      data23 = rbind(data23,a)
    }
    data13 = read.csv(files_distance13[1])
    for(i in 2:length(files_distance13)){
      a = read.csv(files_distance13[i])
      data13 = rbind(data13,a)
    }
    
    
    data12$type = "Chic vs Tsix"
    data23$type = "Tsix vs Linx"
    data13$type = "Chic vs Linx"
    
    data = rbind(data12,data23,data13)
    
    ggplot(data,aes(dist,colour = type))+
      stat_ecdf()+
      xlim(0,as.numeric(input$ymax)) +
      labs("Only in 3D by now")
  })
  
  
  
  pt3_u <- reactive({
    if (!(input$type=="Gyration_radius")) return(NULL)
    
    if(input$dimension == "3D"){ posColumns = 2:4}else{posColumns = 2:3 } #columns correpsonding to the x,y,z position
    
    reconstructed_tracks = read.csv(paste0(dir_uncorrected,"/csv/",input$cell_u))
    channels = names(table(reconstructed_tracks$channel))
    
    ####### separate tracks by channels  (MUST BE MOFIDIED IF NUMBER OF CHANNELS CHANGES)
    channel1 = reconstructed_tracks[reconstructed_tracks$channel==channels[1],c(timeColumn,posColumns)]
    channel2 = reconstructed_tracks[reconstructed_tracks$channel==channels[2],c(timeColumn,posColumns)]
    channel3 = reconstructed_tracks[reconstructed_tracks$channel==channels[3],c(timeColumn,posColumns)]
    
    if(input$d_filter=="ON"){
      ##distances between channels  (MUST BE MOFIDIED IF NUMBER OF CHANNELS CHANGES)
      dist_c12 = dist_channels(channel1,channel2)
      dist_c23 = dist_channels(channel2,channel3)
      dist_c13 = dist_channels(channel1,channel3)
      filter = NULL
      filter = c(filter,derivative_filter(dist_c12,sigma_threshold))
      filter = c(filter,derivative_filter(dist_c23,sigma_threshold))
      filter = c(filter,derivative_filter(dist_c13,sigma_threshold))
      filter = unique(filter)
      reconstructed_tracks[!(reconstructed_tracks$t %in% filter),]
      channel1 = reconstructed_tracks[reconstructed_tracks$channel==channels[1],c(timeColumn,posColumns)]
      channel2 = reconstructed_tracks[reconstructed_tracks$channel==channels[2],c(timeColumn,posColumns)]
      channel3 = reconstructed_tracks[reconstructed_tracks$channel==channels[3],c(timeColumn,posColumns)]
    }
    
    #calculate Gyration radius (mean((rk-rcm)^2))
    alldata = merge(merge(channel1,channel2,by=1),channel3,by=1)
    t = alldata[,1] * tactual
    alldata = alldata[,-1]
    
    if(input$dimension == "3D"){
      allx = alldata[,c(1,4,7)]
      ally = alldata[,c(1,4,7)+1]
      allz = alldata[,c(1,4,7)+2]
      xmean = rowMeans(allx)
      ymean = rowMeans(ally)
      zmean = rowMeans(allz)
      
      rg2 = rowMeans((allx-xmean)**2 + (ally-ymean)**2 + (allz-zmean)**2)
    }else{
      allx = alldata[,c(1,3,5)]
      ally = alldata[,c(1,3,5)+1]
      xmean = rowMeans(allx)
      ymean = rowMeans(ally)
      rg2 = rowMeans((allx-xmean)**2 + (ally-ymean)**2)
    }
    rg = sqrt(rg2)
    
    plot(t,rg2,xlab="time (minutes)",ylab = "Gyration radius (um)", type = "l", main = paste0("filter ",input$d_filter))
    
  })
  
  
  pt4_u <- reactive({
    if (!(input$type=="Auto_cross_pairs")) return(NULL)
    if(input$dimension == "3D"){ posColumns = 2:4}else{posColumns = 2:3 } #columns correpsonding to the x,y,z position
    
    #read data
    reconstructed_tracks = read.csv(paste0(dir_uncorrected,"/csv/",input$cell_u))
    channels = names(table(reconstructed_tracks$channel))
    
    ####### separate tracks by channels  (MUST BE MOFIDIED IF NUMBER OF CHANNELS CHANGES)
    channel1 = reconstructed_tracks[reconstructed_tracks$channel==channels[1],c(timeColumn,posColumns)]
    channel2 = reconstructed_tracks[reconstructed_tracks$channel==channels[2],c(timeColumn,posColumns)]
    channel3 = reconstructed_tracks[reconstructed_tracks$channel==channels[3],c(timeColumn,posColumns)]
    
    ##distances between channels  (MUST BE MOFIDIED IF NUMBER OF CHANNELS CHANGES)
    dist_c12 = dist_channels(channel1,channel2)
    dist_c23 = dist_channels(channel2,channel3)
    dist_c13 = dist_channels(channel1,channel3)
    dist_c12$t = dist_c12$t*tactual
    dist_c23$t = dist_c23$t*tactual
    dist_c13$t = dist_c13$t*tactual
    
    
    #filter
    if(input$d_filter=="ON"){
      filter = NULL
      filter = c(filter,derivative_filter(dist_c12,sigma_threshold))
      filter = c(filter,derivative_filter(dist_c23,sigma_threshold))
      filter = c(filter,derivative_filter(dist_c13,sigma_threshold))
      filter = unique(filter)
      
      dist_c12 = dist_c12[!(dist_c12$t %in% filter), ]
      dist_c23 = dist_c23[!(dist_c23$t %in% filter), ]
      dist_c13 = dist_c13[!(dist_c13$t %in% filter), ]
    }
    
    a = merge(merge(dist_c12,dist_c13,by=1),dist_c23,by=1)
    time = a[,1]*tactual
    a = a[,-1]
    colnames(a)= c("Chic vs Tsix","Chic vs Linx","Tsix vs Linx")
    
    #auto correlation
    diagonal.acf = function(x){
      par(new = TRUE)
      auto = acf(x,plot=FALSE)
      plot(time[1:length(auto[[1]])]-min(time),auto[[1]],xlab="lag (Minutes)",ylim=c(0,1),ylab="Auto correlation",type="l")
    }
    
    #cross correlation
    upper.ccf = function(x,y){
      par(new = TRUE)
      auto = ccf(x,y,plot=FALSE)
      plot(c(sort(-time[1:(as.integer(length(auto[[1]])/2)+1)]),
             time[(1:as.integer(length(auto[[1]])/2))+1]-min(time)),
           auto[[1]],xlab="lag (Minutes)",ylab="Cross correlation",type="l",ylim=c(0,1))
    }
    
    panel.smooth <- function(x,y) {par(new=TRUE);
      smoothScatter(x,y, nrpoints=0);
      abline(a=0,b=1,col=2)}
    
    pairs(a,diag.panel = diagonal.acf,
          upper.panel = upper.ccf,
          lower.panel = panel.smooth)  
    
    
  })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Return the requested graph
  graphInput <- reactive({
    switch(input$type,
           "Pair-wise_dist" = pt1(),
           "ECDF_all_data" = pt2(),
           "Gyration_radius" = pt3(),
           "Auto_cross_pairs" = pt4()
    )
  })
  
  graphInput1 <- reactive({
    switch(input$type,
           "Pair-wise_dist" = pt1_u(),
           "ECDF_all_data" = pt2_u(),
           "Gyration_radius" = pt3_u(),
           "Auto_cross_pairs" = pt4_u()
    )
  })
  
  output$plotgraph1 <- renderPlot({ 
    graphInput()
  })
  
  output$plotgraph2 <- renderPlot({ 
    graphInput1()
  })
  
  output$downloadData <- downloadHandler(
    filename = function(){input$cell},
    content = function(file) {
      reconstructed_tracks = read.csv(paste0(dir,"/csv/",input$cell))
      write.csv(reconstructed_tracks, file, row.names = F)
    }
  )
  
  
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste(input$cell, '.png', sep='') },
    content = function(file) {
      ggsave(file,graphInput(),device = "png")
    }
  )
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

