
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
library(plotly)
library(reshape2)
##global options
timeColumn = 1 # corresponds to column with time position information
tactual = 1 #minutes
melements = 10 #minimum timepoints for msd calculation 
mdelay = 150 # maximum delay
sigma_threshold = 4 #how many standard deviation away from the mean of the derivative of distance with respect to time to exclude because of wrong stitching
dir = "/tungstenfs/scratch/ggiorget/_LIVECELL/Analysis_Data/F6_FullSeq_SingleCells//results_analysis_corrected_traj/2h30s_movies/" 
#dir_uncorrected = "/tungstenfs/scratch/ggiorget/_LIVECELL/Analysis_Data/PGK_G8_A11_B3_F6_SingleCells/results_analysis_corrected_traj" 
#end

files = list.files(paste0(dir,"/csv/"),pattern = "_reconstructed_tracks.csv")
#files_uncorrected = list.files(paste0(dir_uncorrected,"/csv/"),pattern = "_reconstructed_tracks.csv")


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
           HTML("<li>Filtering on pair-wise distance is based on speed of movement</li>"),
           HTML("<li>MSD can be done using single channel (position) or using distance between channels (distance)</li>"),
           HTML("<li>MSD calculation can be done using lower in silico time resolution</li>"),
           HTML("<li>MSD calculation for all the cells is done using the actual time resolution of imaging experiment</li>"),
           HTML("<li>Velocity autocorrelation is based on radial velocity (derivative of distance between two channels)</li>"),
           HTML("<li>First passage time, duration (non) contact are based on all the cells </li>"),
           HTML("<li>Pairs: upper=ccf; diagonal = acf; low=scatter</li><ul>")),
          
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        selectInput("d_filter", "Filter (movement)",  choices=c("ON","OFF")),
        selectInput("dimension","2 or 3D", choices = c("3D","2D")),
        selectInput("cell", "Which cell:",  choices=files),
        selectInput("type", "Which Plot",  choices=c("Pair-wise_dist","ECDF_all_data","Gyration_radius","Auto_cross_pairs","MSD","MSD_allCells","Autocorr_velocity_all_cells","first_passage_time_distribution","duration_contact")),
        selectInput("msd", "Type of MSD",  choices=c("distance","position")),
        selectInput("tres", "time resolution for MSD:",  choices=c((1:50)*tactual)),
        selectInput("dist_thresh", "Threshold on distance for contact (um):",  choices=c((10:50)/100)),
        selectInput("autocorrelation", "Type of autocorrelation:",  choices=c("distance","velocity")),
        selectInput("ymax", "Max value in y axis:",  choices=c(1.5,(1:10)/10,1,2,2.5,3,4,5)),
        downloadButton("downloadData", "Reconstructed trajectory"),
        downloadButton("downloadPlot", "Plot"),
        width = 15
      ),
      
      # Show a plot of the generated distribution
      mainPanel(plotOutput("plotgraph"),plotlyOutput("plot"))
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
    
    dist_c12$type = "Chic vs Xite"
    dist_c23$type = "Xite vs Linx"
    dist_c13$type = "Chic vs Linx"
    
    data = rbind(dist_c12,dist_c23,dist_c13)
    ggplot(data,aes(x = t,y=dist,colour=type)) +
      geom_line(size=0.5)+
      geom_point(size=0.2) +
      ylim(0,as.numeric(input$ymax)) +
      labs(paste0(input$dimension, " Pair-wise_distances filter ",input$d_filter)) +
      xlab("t in min")+
      ylab("Distance in µm")
    
    
    # plot(dist_c12[,c(1,2)],type="l",ylim=c(0,2), ylab="Distance in µm", xlab="t in min", main=paste0(input$dimension, " Pair-wise_distances filter ",input$d_filter), col=3)
    # lines(dist_c23[,c(1,2)],type="l",col=4)
    # lines(dist_c13[,c(1,2)],type="l",col=2)
    # legend(x="topright",legend=c("Chic vs Xite","Xite vs Linx","Chic vs Linx"),col=c(3,4,2),lty=c(1,1))
    # 
  })
  
  
  ##############################################################################################################################################################################################################################################################################################
  
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
    
    
    data12$type = "Chic vs Xite"
    data23$type = "Xite vs Linx"
    data13$type = "Chic vs Linx"
    
    data = rbind(data12,data23,data13)
    
    ggplot(data,aes(dist,colour = type))+
      stat_ecdf()+
      xlim(0,as.numeric(input$ymax)) +
      labs("Only in 3D by now")
  })
  
  ##############################################################################################################################################################################################################################################################################################
  
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
    t=alldata[,1]
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
    
    df = data.frame(t = t, rg = rg)
    ggplot(df,aes(x=t,y=rg)) + 
      geom_line() +
      xlab("time (minutes)") +
      ylab("Gyration radius (um)") +
      labs(paste0("filter ",input$d_filter))
    #plot(t,rg2,xlab="time (minutes)",ylab = "Gyration radius (um)", type = "l", main = paste0("filter ",input$d_filter))
    
    })

  ##############################################################################################################################################################################################################################################################################################
  
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
    time = a[,1]
    a = a[,-1]
    colnames(a)= c("Chic vs Xite","Chic vs Linx","Xite vs Linx")
    
    #auto correlation
    diagonal.acf = function(x){
      par(new = TRUE)
      auto = acf(x,plot=FALSE)
      plot(time[1:length(auto[[1]])]-min(time),auto[[1]],xlab="lag (Minutes)",ylim=c(0,1),ylab="Auto correlation",type="l")
    }
    
    #velocity autocorrelation
    diagonal.acf_velocity = function(x){
      par(new = TRUE)
      y = diff(x)/diff(time)
      auto = acf(y,plot=FALSE)
      t = time[-length(time)]
      plot(t[1:length(auto[[1]])]-min(t),auto[[1]],xlab="lag (Minutes)",ylim=c(-0.5,1),ylab="Auto correlation Velocity",type="l")
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
    
    if(input$autocorrelation == "distance"){
      pairs(a,diag.panel = diagonal.acf,
          upper.panel = upper.ccf,
          lower.panel = panel.smooth)  
    }
    if(input$autocorrelation == "velocity"){
      pairs(a,diag.panel = diagonal.acf_velocity,
            upper.panel = upper.ccf,
            lower.panel = panel.smooth)  
    }
    
      
  })
  
  ##############################################################################################################################################################################################################################################################################################
  
  ##MSD on distances btw loci
  pt5 <- reactive({
    if (!(input$type=="MSD")) return(NULL)
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
    
    dist_c12$type = "Chic vs Xite"
    dist_c23$type = "Xite vs Linx"
    dist_c13$type = "Chic vs Linx"
    
    if(input$msd=="position"){
      poscolumns=c(2:4)
      df1 = calculate_msd_trajectory(channel1,poscolumns,input$tres,mdelay,melements)
      df1$type = "Chic1"
      df2 = calculate_msd_trajectory(channel2,poscolumns,input$tres,mdelay,melements) 
      df2$type = "Xite"
      df3 = calculate_msd_trajectory(channel3,poscolumns,input$tres,mdelay,melements) 
      df3$type = "Linx"
    }
    if(input$msd=="distance"){
      poscolumns=2
      df1 = calculate_msd_trajectory(dist_c12,poscolumns,input$tres,mdelay,melements) 
      df1$type = "Chic vs Xite"
      df2 = calculate_msd_trajectory(dist_c23,poscolumns,input$tres,mdelay,melements) 
      df2$type = "Xite vs Linx"
      df3 = calculate_msd_trajectory(dist_c13,poscolumns,input$tres,mdelay,melements) 
      df3$type = "Chic vs Linx"
    }
    
    df = rbind(df1,df2,df3)
    df$time = log10(df$time)
    df$disp = log10(df$disp)
    
    ggplot(df,aes(x=time,y=disp,col=type)) + 
      geom_line() + geom_point()+
      xlab("time (minutes, log10)") +
      ylab(paste0("MSD on ",input$msd, " um (log10)")) +
      labs(paste0("filter ",input$d_filter))
    
    
    })
  
  ##############################################################################################################################################################################################################################################################################################
  
  pt6 <- reactive({
    if (!(input$type=="MSD_allCells")) return(NULL)
    
    files_distance12 = list.files(paste0(dir,"/csv/"),pattern = "distance_MSD_c12.csv",full.names = T)
    files_distance23 = list.files(paste0(dir,"/csv/"),pattern = "distance_MSD_c23.csv",full.names = T)
    files_distance13 = list.files(paste0(dir,"/csv/"),pattern = "distance_MSD_c13.csv",full.names = T)
    
    
    
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
    
    data12 = aggregate(data12$disp,list(data12$time),median)
    data23 = aggregate(data23$disp,list(data23$time),median)
    data13 = aggregate(data13$disp,list(data13$time),median)
    
    data12$type = "Chic vs Xite"
    data23$type = "Xite vs Linx"
    data13$type = "Chic vs Linx"
    
    data = rbind(data12,data23,data13)
    
    data$Group.1 = log10(data$Group.1)
    data$x = log10(data$x)
    
    ggplot(data,aes(x=Group.1, y=x,colour = type))+
      geom_line() + labs("MSD all cells on distance filter ON") + geom_point()+
      xlab("time (minutes, log10) ") + ylab("MSD (log10, um)")  + labs("MSD calculated using imaging resolution time")
    
    #ggplot(data,aes(x=time, y=disp,colour = type))+
    #  stat_smooth(method="loess", span=0.1, se=TRUE, aes(fill=type), alpha=0.3) + labs("MSD all cells on distance filter ON") + 
    #  xlab("time (minutes)") + ylab("MSD")
    
  })
  
  ##############################################################################################################################################################################################################################################################################################
  
  pt7 <- reactive({
    if (!(input$type=="Autocorr_velocity_all_cells")) return(NULL)
    
    files_distance12 = list.files(paste0(dir,"/csv/"),pattern = "12velocity_autocorrelation.csv",full.names = T)
    files_distance13 = list.files(paste0(dir,"/csv/"),pattern = "13velocity_autocorrelation.csv",full.names = T)
    files_distance23 = list.files(paste0(dir,"/csv/"),pattern = "23velocity_autocorrelation.csv",full.names = T)
    
    
    
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
    
    #data12 = aggregate(data12$autocorr,list(data12$t),mean)
    #data23 = aggregate(data23$autocorr,list(data23$t),mean)
    #data13 = aggregate(data13$autocorr,list(data13$t),mean)
    
    data12$type = "Chic vs Xite"
    data23$type = "Xite vs Linx"
    data13$type = "Chic vs Linx"
    
    data = rbind(data12,data23,data13)
    
    #ggplot(data,aes(x=Group.1, y=x,colour = type))+
    #  geom_line() + labs("All cells on distance filter ON") + 
    #  xlab("time delay (minutes)") + ylab("Velocity Autocorrelation")
    ggplot(data,aes(x=t,y=autocorr,colour=type)) + 
      stat_smooth(method="loess", span=0.1, se=TRUE, aes(fill=type), alpha=0.3) +
      theme_bw() + xlab("time delay (minutes)") + ylab("Velocity Autocorrelation") + labs("velocity calculated using imaging resolution time")
    })
  
  ##############################################################################################################################################################################################################################################################################################
  
  pt8 <- reactive({
    if (!(input$type=="first_passage_time_distribution")) return(NULL)
    
    if(input$d_filter=="ON"){
      files_distance12 = list.files(paste0(dir,"/csv/"),pattern = "12distances_filtered.csv",full.names = T)
      files_distance23 = list.files(paste0(dir,"/csv/"),pattern = "23distances_filtered.csv",full.names = T)
      files_distance13 = list.files(paste0(dir,"/csv/"),pattern = "13distances_filtered.csv",full.names = T)
    }else{
      files_distance12 = list.files(paste0(dir,"/csv/"),pattern = "12distances.csv",full.names = T)
      files_distance23 = list.files(paste0(dir,"/csv/"),pattern = "23distances.csv",full.names = T)
      files_distance13 = list.files(paste0(dir,"/csv/"),pattern = "13distances.csv",full.names = T)
    }
    
    data12=NULL
    for(i in 1:length(files_distance12)){
      a = read.csv(files_distance12[i])
      #find timepoints below threshold
      sel = which(a[,2]<input$dist_thresh)
      if(length(sel)>0){
        #select the first timepoint and calculate the difference with the time 0
        data12 = c(data12,(a[sel[1],1]-a[1,1]))
      }
    }
    
    data23=NULL
    for(i in 1:length(files_distance23)){
      a = read.csv(files_distance23[i])
      #find timepoints below threshold
      sel = which(a[,2]<input$dist_thresh)
      if(length(sel)>0){
        #select the first timepoint and calculate the difference with the time 0
        data23 = c(data23,(a[sel[1],1]-a[1,1]))
      }
    }
    
    data13=NULL
    for(i in 1:length(files_distance13)){
      a = read.csv(files_distance13[i])
      sel = which(a[,2]<input$dist_thresh)
      if(length(sel)>0){
        data13 = c(data13,(a[sel[1],1]-a[1,1]))
      }
    }
    
    #create dataframe
    data12 = data.frame(fpt = data12)
    data23 = data.frame(fpt = data23)
    data13 = data.frame(fpt = data13)
    
    data12$type = "Chic vs Xite"
    data23$type = "Xite vs Linx"
    data13$type = "Chic vs Linx"
    
    data = rbind(data12,data23,data13)
    
    ggplot(data,aes(fpt,colour = type,fill=type,alpha=0.3))+
      geom_histogram(alpha=0.3)+
      labs("First passage time") + xlab("time (minutes)")
  })
  
  ##############################################################################################################################################################################################################################################################################################
  
  pt9 <- reactive({
    if (!(input$type=="duration_contact")) return(NULL)
    
    if(input$d_filter=="ON"){
      files_distance12 = list.files(paste0(dir,"/csv/"),pattern = "12distances_filtered.csv",full.names = T)
      files_distance23 = list.files(paste0(dir,"/csv/"),pattern = "23distances_filtered.csv",full.names = T)
      files_distance13 = list.files(paste0(dir,"/csv/"),pattern = "13distances_filtered.csv",full.names = T)
    }else{
      files_distance12 = list.files(paste0(dir,"/csv/"),pattern = "12distances.csv",full.names = T)
      files_distance23 = list.files(paste0(dir,"/csv/"),pattern = "23distances.csv",full.names = T)
      files_distance13 = list.files(paste0(dir,"/csv/"),pattern = "13distances.csv",full.names = T)
    }
    
    data12=NULL
    for(i in 1:length(files_distance12)){
      a = read.csv(files_distance12[i])
      #find timepoints below threshold
      sel = which(a[,2]<input$dist_thresh)
      
      if(length(sel)>1){
        #extract timepoints below threhsolds
        a = a[sel,]
        #count number of consecutive timepoints below threhsolds (a true correspond to a pair of consecutive timepoints where distance is below threshold. a false correspond to a single timepoint where it is close)
        b = rle(diff(a[,1])==tactual)
        #calculating for each strech of close, the total time being close
        time_close = (b$lengths[b$value==TRUE]+1)
        #number of time they are close for below tactual time (# of falses in b)
        n = length(sel)-sum(b$lengths[b$value==TRUE])
        data12 = c(data12,c(time_close,rep(tactual,n)))
      }
    }
    
    
    data23=NULL
    for(i in 1:length(files_distance23)){
      a = read.csv(files_distance23[i])
      #find the first timepoint below threshold
      sel = which(a[,2]<input$dist_thresh)
      
      if(length(sel)>1){
        #extract timepoints below threhsolds
        a = a[sel,]
        #count number of consecutive timepoints below threhsolds (a true correspond to a pair of consecutive timepoints where distance is below threshold. a false correspond to a single timepoint where it is close)
        b = rle(diff(a[,1])==tactual)
        #calculating for each strech of close, the total time being close
        time_close = (b$lengths[b$value==TRUE]+1)
        #number of time they are close for below tactual time (# of falses in b)
        n = length(sel)-sum(b$lengths[b$value==TRUE])
        data23 = c(data23,c(time_close,rep(tactual,n)))
      }
    }
    
    data13=NULL
    for(i in 1:length(files_distance13)){
      a = read.csv(files_distance13[i])
      #find the first timepoint below threshold
      sel = which(a[,2]<input$dist_thresh)
      
      if(length(sel)>1){
        #extract timepoints below threhsolds
        a = a[sel,]
        #count number of consecutive timepoints below threhsolds (a true correspond to a pair of consecutive timepoints where distance is below threshold. a false correspond to a single timepoint where it is close)
        b = rle(diff(a[,1])==tactual)
        #calculating for each strech of close, the total time being close
        time_close = (b$lengths[b$value==TRUE]+1)
        #number of time they are close for below tactual time (# of falses in b)
        n = length(sel)-sum(b$lengths[b$value==TRUE])
        data13 = c(data13,c(time_close,rep(tactual,n)))
      }
    }
    
    data12 = data.frame(duration_contact = data12)
    data23 = data.frame(duration_contact = data23)
    data13 = data.frame(duration_contact = data13)
    
    data12$type = "Chic vs Xite"
    data23$type = "Xite vs Linx"
    data13$type = "Chic vs Linx"
    
    data = rbind(data12,data23,data13)
    
    ggplot(data,aes(duration_contact,colour = type,fill=type,alpha=0.3))+
      geom_histogram(alpha=0.3)+
      labs("Duration_contact") +xlab("duration (in minutes)")
  })
  
  ##############################################################################################################################################################################################################################################################################################
  
  # Return the requested graph
  graphInput <- reactive({
    switch(input$type,
           "Pair-wise_dist" = pt1(),
           "ECDF_all_data" = pt2(),
           "Gyration_radius" = pt3(),
           "Auto_cross_pairs" = pt4(),
           "MSD" = pt5(),
           "MSD_allCells" = pt6(),
           "Autocorr_velocity_all_cells" = pt7(),
           "first_passage_time_distribution" = pt8(),
           "duration_contact" = pt9()
    )
  })
  
  output$plotgraph <- renderPlot({ 
    graphInput() 
  })
  
  
  output$plot <- renderPlotly({
    if(input$dimension == "3D"){ posColumns = 2:4}else{posColumns = 2:3 }
    reconstructed_tracks = read.csv(paste0(dir,"/csv/",input$cell))
    ##plot 3D trajectories
    
    channels = names(table(reconstructed_tracks$channel))
    
    ####### separate tracks by channels  (MUST BE MOFIDIED IF NUMBER OF CHANNELS CHANGES)
    channel1 = reconstructed_tracks[reconstructed_tracks$channel==channels[1],c(timeColumn,posColumns)]
    channel2 = reconstructed_tracks[reconstructed_tracks$channel==channels[2],c(timeColumn,posColumns)]
    channel3 = reconstructed_tracks[reconstructed_tracks$channel==channels[3],c(timeColumn,posColumns)]
    
    ##distances between channels  (MUST BE MOFIDIED IF NUMBER OF CHANNELS CHANGES)
    dist_c12 = dist_channels(channel1,channel2)
    dist_c23 = dist_channels(channel2,channel3)
    dist_c13 = dist_channels(channel1,channel3)

    
    
    #filter
    if(input$d_filter=="ON"){
      filter = NULL
      filter = c(filter,derivative_filter(dist_c12,sigma_threshold))
      filter = c(filter,derivative_filter(dist_c23,sigma_threshold))
      filter = c(filter,derivative_filter(dist_c13,sigma_threshold))
      filter = unique(filter)
      
      reconstructed_tracks = reconstructed_tracks[!(reconstructed_tracks[,1] %in% filter),]
    }
    
    
    plot_ly(reconstructed_tracks, x = ~x, y = ~y, z = ~z, split = ~color ,type = 'scatter3d', mode = 'lines',
            line = list(width = 2))%>%
      layout(
        title = channels,
        scene = list(
          xaxis = list(title = "X position"),
          yaxis = list(title = "Y position"),
         zaxis = list(title = "Z position")
        ))
  })
  
  output$downloadData <- downloadHandler(
    filename = function(){input$cell},
    content = function(file) {
      reconstructed_tracks = read.csv(paste0(dir,"/csv/",input$cell))
      write.csv(reconstructed_tracks, file, row.names = F)
    }
  )
  

  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste(input$cell, '.pdf', sep="") },
    content = function(file) {
      ggsave(filename = file,plot = graphInput(),device = "pdf")
    },
    contentType	= "image/pdf"
  )
  
    
}

# Run the application 
shinyApp(ui = ui, server = server)

