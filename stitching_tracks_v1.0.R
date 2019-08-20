####To Do List for this script:
#better define input data based on filenames and clean up the data (renaming, wait for Jan)
#figure out a smarter way to "weigh" the scores for time intervals that are longer (done)
#write scripts for plots that we mentioned (correlation of pairwise 3D distances, etc.)
#make Shiny app (done)

###Set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("functions_v0.R")
library(plot3D)
library(reshape2)
library(plotly)

###Global options (you can touch)
####CHOOSE DATASET HERE BY CHOOSING FILEPATH:
dirpath = "/tungstenfs/scratch/ggiorget/_LIVECELL/Analysis_Data/PGK_G8_A11_B3_F6_SingleCells/Corrected_traj/2h_movies/" #set the filepath to files you are intested in
##directory containing analysis result, must contain csv directory and pdf directory
outputdir = "/tungstenfs/scratch/ggiorget/_LIVECELL/Analysis_Data/PGK_G8_A11_B3_F6_SingleCells/results_analysis_corrected_traj/" 


trackColumn = 7 #corresponds to column with TRACKID information 
timeColumn = 1 # corresponds to column with time position information
posColumns = 2:4 #columns correpsonding to the x,y,z position
channelColumn = 6 #corresponds to the column assigning the channel
#### SET FOR EVERY NEW DATASET!!! ####
tactual = 0.5 #what is the time interval in min SET FOR EVERY NEW DATASET!!!!
overlap = 0.8 #how much of a new track must be unique (not overlapping with any other track) as fraction
#distances in microns!!!
distthreshold = 4 #max distance between 2 tracks
#Somehow the assignments here are wrong, figure this out later
t_dist_chan = 1 #threshold to penalize distance between different channels at certain timepoints
t_dist_inchan = 1 # threshold to penalize distance within the same channel between different tracks
t_dist_time = 1 # threshold to penalize distance in time between tracks in the same channel
numtrack = 0 # threshold to decide on the minimum length of track that should be included


sigma_threshold = 4 #how many standard deviation away from the mean of the derivative of distance with respect to time to exclude because of wrong stitching
##end global options


##do not touch below if you are not aware!

###Import data file and make it useful for working in R (this could be prettier)
filelist = list.files(path=dirpath, pattern="\\.csv$",recursive = TRUE)
for(fileno in 1:length(filelist)){
  filename = filelist[fileno]
  traj <- read.csv(file = paste(dirpath,filename, sep=""), sep = ",", header = TRUE)
  #clean up some of the dataset columns
  traj$RowID <- gsub("_w[0-9]_[0-9]+_[0-9]+", "", x=traj$RowID)
  traj$Original.RowID  = gsub("cell[0-9]_w", "C", x=as.character(traj$Original.RowID))
  names(traj)[names(traj) == "Original.RowID"] <- "channel"
  names(traj)[names(traj) == "particle"] <- "TRACK_ID"
  names(traj)[names(traj) == "RowID"] <- "Cell_no"
  
  
  
  
  ###Analysis
  time = names(table(traj[,timeColumn]))
  channels = names(table(traj[,channelColumn]))
  
  tracks = NULL
  ##find the closest pair at given time
  for (i in 1:length(time)) {
    #extract data for the timepoint corresponding to i
    data_c1 = traj[traj$t==time[i] & traj$channel==channels[1],]
    data_c2 = traj[traj$t==time[i] & traj$channel==channels[2],]
    data_c3 = traj[traj$t==time[i] & traj$channel==channels[3],]
    
    if(sum(nrow(data_c1), nrow(data_c2), nrow(data_c3))<=1){
      track_3d <- as.vector(rep("NA",6))
    }else{
      #Calculate pair-wise distances and return the track number corresponding to the minimum distance
      if(nrow(data_c1)==0 || nrow(data_c2)==0){x1 <- as.vector(rep("NA", 2))}else{x1 <- as.vector(pair_dis(data_c1, data_c2))}
      if(nrow(data_c1)==0 || nrow(data_c3)==0){x2 <- as.vector(rep("NA", 2))}else{x2 <- as.vector(pair_dis(data_c1, data_c3))}
      if(nrow(data_c2)==0 || nrow(data_c3)==0){x3 <- as.vector(rep("NA", 2))}else{x3 <- as.vector(pair_dis(data_c2, data_c3))}
      #combine to trackset for all three channels and make array out of 
      track_3d <- c(x1, x2, x3) 
    }
    tracks <- c(tracks,track_3d)
  }
  
  if(!file.exists(paste0(outputdir,"/csv/",gsub("/","_",filename),"_reconstructed_tracks.csv")) | !file.exists(paste0(outputdir,"/csv/",gsub("/","_",filename),"_tracks_used.csv"))){
    reconstructed_tracks = NULL
    tracks_used = NULL
    #how many times a track is found to be the closest to another channel
    count_closest = table(c(tracks))
    
    
    for(ch in 1:length(channels)){ #run on channel
      timeappo = time
      selected_channel=channels[ch] #choose channel
      ##sort based on closeness with other channels
      selected_channel_count_closest = sort(count_closest[grep(selected_channel,names(count_closest))],decreasing=T)
      ##choose the best track to start with (the one that is most often clostest to another track in other channels), then use the sorting based on both differences between cms and the sorting based on distances in channels (this to avoid jumps when reconstructing)
      #sorting also based on "distance" in time
      sorted_channels_and_distanceWithinChannel = selected_channel_count_closest 
      k=1
      while(k<=length(sorted_channels_and_distanceWithinChannel)){ #run on tracks
        closest_id = strsplit(names(sorted_channels_and_distanceWithinChannel[k]),"_")[[1]][2] #current longest track
        df = traj[traj$channel==selected_channel & traj$TRACK_ID==closest_id,] #extracting current longest track
        if(sum(df$t %in% timeappo)/nrow(df)>=overlap & nrow(df)>=numtrack) { #if current longest track is not overlapping with previous longer tracks then save it, also select for tracks that are longer than three timepoints
          df = df[df$t %in% timeappo,]
          reconstructed_tracks = rbind(reconstructed_tracks,df) #save current longest track
          timeappo = timeappo[!(timeappo %in% df$t)] #eliminate time already taken by previous tracks
          selected_channel_count_closest = selected_channel_count_closest[-k] #eliminate track
          if(length(selected_channel_count_closest)==0){break}
          #resort order of tracks based on closeness with previous track and the closeness with other channels: this is to avoid jumps when reconstructing tracks
          #order_based_on_channel_distance = 1:length(selected_channel_count_closest) #rank based on distance with other channel
          
          shifted = selected_channel_count_closest - min(selected_channel_count_closest)
          order_based_on_channel_distance = as.vector(1 - (shifted)/max(shifted))
          
          names(order_based_on_channel_distance) = names(selected_channel_count_closest) #rank based on distance with other channel
          
          remaining_tracks = traj[traj$channel==selected_channel & traj$TRACK_ID!=closest_id,] #extract remaining tracks to find rank based of distance between centre of mass (top rank = closest with the last track): this is to avoid jumps when reconstructing tracks
          distance_cms = NULL
          distance_time = NULL
          for(tids in unique(remaining_tracks$TRACK_ID)){
            d = NULL
            dt = NULL
            for(tidsrecorder in unique(reconstructed_tracks$TRACK_ID)){
              dt = abs(mean(reconstructed_tracks$t[reconstructed_tracks$TRACK_ID == tidsrecorder])-mean((remaining_tracks[remaining_tracks$TRACK_ID==tids,"t"]))) #distance in time
              d = c(d,dist(rbind(cm(reconstructed_tracks[reconstructed_tracks$TRACK_ID == tidsrecorder,posColumns]),cm(remaining_tracks[remaining_tracks$TRACK_ID==tids,posColumns]))))#calculate distance btw centre of mass between last saved track and all the remaining
            }
            d = min(d)
            names(d)= paste(selected_channel,tids,sep="_")
            distance_cms = c(distance_cms,d) #save distance btw centre of mass
            
            dt = min(dt)
            names(dt) =  paste(selected_channel,tids,sep="_")
            distance_time = c(distance_time,dt) 
          }
          distance_cms = sort(distance_cms) #sort from lowest to highest 
          distance_cms = distance_cms[distance_cms<distthreshold] #filter based on treshold
          distance_time = sort(distance_time) #sort from lowest to highest 
          
          #order_based_on_dist_cm_same_channel = 1:length(distance_cms) #order based on distances btw centre of mass
          shifted = distance_cms - min(distance_cms)
          order_based_on_dist_cm_same_channel = shifted/max(shifted)
          
          names(order_based_on_dist_cm_same_channel) = names(distance_cms)
          
          #order_based_on_dist_t_same_channel = 1:length(distance_time) #order based on distances btw time
          shifted = distance_time - min(distance_time)
          order_based_on_dist_t_same_channel = shifted/max(shifted)
          
          names(order_based_on_dist_t_same_channel) = names(distance_time)
         
          combine_order = merge((merge(order_based_on_dist_cm_same_channel,order_based_on_channel_distance,by=0,all=TRUE)),order_based_on_dist_t_same_channel,by.x=1,by.y=0) #combine the two ordering
          na_replacement = max(combine_order[,c(2,3)],na.rm = T)
          combine_order[is.na(combine_order)]=na_replacement 
          #x = distance channels; y.x = distance within channel; y.y = distance in time within channel, can be used to imply penalties on specific measurements, I think the assignments here are wrong
          final_order = t_dist_chan*combine_order$x+ t_dist_inchan*combine_order$y.x + t_dist_time*combine_order$y.y #order that consider both distance between cms and distance between channels
          names(final_order)=combine_order$Row.names
          sorted_channels_and_distanceWithinChannel = sort(final_order)
          ##end resort
          
          tracks_used=rbind(tracks_used,c(selected_channel,closest_id))
          k=0 #restart
        }
        k=k+1
      }
    }
    
    
    reconstructed_tracks = reconstructed_tracks[order( reconstructed_tracks[,5], reconstructed_tracks[,1] ),]
    reconstructed_tracks$color = as.numeric(gsub("C","",reconstructed_tracks$channel))
    write.csv(file = paste0(outputdir,"/csv/",gsub("/","_",filename),"_reconstructed_tracks.csv"),x=reconstructed_tracks,row.names = F)
    write.csv(file = paste0(outputdir,"/csv/",gsub("/","_",filename),"_tracks_used.csv"),x=tracks_used,row.names = F)
  }else{
    reconstructed_tracks = read.csv(paste0(outputdir,"/csv/",gsub("/","_",filename),"_reconstructed_tracks.csv"))
    tracks_used = read.csv(paste0(outputdir,"/csv/",gsub("/","_",filename),"_tracks_used.csv"))
  }
  
  

  
  # ##plot 3D trajectories
  # plot_ly(reconstructed_tracks, x = ~x, y = ~y, z = ~z, split = ~color ,type = 'scatter3d', mode = 'lines',
  #         line = list(width = 2))%>%
  #   layout(
  #     title = channels[ch],
  #     scene = list(
  #       xaxis = list(title = "X position"),
  #       yaxis = list(title = "Y position"),
  #      zaxis = list(title = "Z position")
  #     ))
  
  
  
  pdf(paste0(outputdir,"/pdf/",gsub("/","_",filename),".pdf"))
  
  ######quality control
  for(ch in channels){
    i=1
    traj_check = traj[traj$channel==ch,]
    for(id in unique(traj_check$TRACK_ID)){
      subset = traj_check[traj_check$TRACK_ID==id,]
      if(id %in% reconstructed_tracks$TRACK_ID[reconstructed_tracks$channel==ch]){col=2}else{col=1}
      if(i==1){
        plot(subset$t,rep(1+i/80,nrow(subset)),xlim=c(min(as.numeric((time))),max(as.numeric(time))),col=col,ylim=c(0.8,2),type="l",xlab="time",ylab="",main=paste0("channel = ",ch))
      }else{
        lines(subset$t,rep(1+i/80,nrow(subset)),col=col)
      }
      i=i+1
    }
    legend(x="topleft",legend=c("discarded","kept"),col=1:2,lty=c(1,1))
  }
  
  ####### separate tracks by channels  (MUST BE MOFIDIED IF NUMBER OF CHANNELS CHANGES)
  channel1 = reconstructed_tracks[reconstructed_tracks$channel==channels[1],c(timeColumn,posColumns)]
  channel2 = reconstructed_tracks[reconstructed_tracks$channel==channels[2],c(timeColumn,posColumns)]
  channel3 = reconstructed_tracks[reconstructed_tracks$channel==channels[3],c(timeColumn,posColumns)]
  
  
  #Plot the tracks
  scatter3D(channel1$x, channel1$y, channel1$z, type="l", col=3, ticktype="detailed", xlab="Position in X in µm", ylab="Position in Y in µm", zlab="Position in Z in µm", main="Trajectories of three loci over time (complete data set)")
  scatter3D(channel2$x, channel2$y, channel2$z, type="l", col=4, add = TRUE)
  scatter3D(channel3$x, channel3$y, channel3$z, type="l", col=2, add = TRUE)
  legend(x="topright",legend=c("Chic","Tsix","Linx"),col=c(3,4,2),lty=c(1,1))
  
  ##distances between channels  (MUST BE MOFIDIED IF NUMBER OF CHANNELS CHANGES)
  dist_c12 = dist_channels(channel1,channel2)
  dist_c23 = dist_channels(channel2,channel3)
  dist_c13 = dist_channels(channel1,channel3)
  dist_c12$t = dist_c12$t*tactual
  dist_c23$t = dist_c23$t*tactual
  dist_c13$t = dist_c13$t*tactual
  
  #saving distances
  write.csv(file = paste0(outputdir,"/csv/",gsub("/","_",filename),"_12distances.csv"),x=dist_c12,row.names = F)
  write.csv(file = paste0(outputdir,"/csv/",gsub("/","_",filename),"_23distances.csv"),x=dist_c23,row.names = F)
  write.csv(file = paste0(outputdir,"/csv/",gsub("/","_",filename),"_13distances.csv"),x=dist_c13,row.names = F)
  
  
  
  ##filtering, deleting too sharping increase and decrease  (MUST BE MOFIDIED IF NUMBER OF CHANNELS CHANGES)
  #d = 2 column matrix (1 column denominator, second column numerator of derivation), threshold on sd to filter
  #output the time point you have to filter
  derivative_filter = function(d,threshold){
    derivative = diff(d[,2])/diff(d[,1])
    plot(derivative)
    abline(h=c((mean(derivative))+threshold*sd(derivative),(mean(derivative))-threshold*sd(derivative)))
    return(d[(which(abs(derivative)>(mean(derivative))+threshold*sd(derivative))+1),1])
  }
  
  filter = NULL
  filter = c(filter,derivative_filter(dist_c12,sigma_threshold))
  filter = c(filter,derivative_filter(dist_c23,sigma_threshold))
  filter = c(filter,derivative_filter(dist_c13,sigma_threshold))
  filter = unique(filter)
  
  dist_c12_filtered = dist_c12[!(dist_c12$t %in% filter), ]
  dist_c23_filtered = dist_c23[!(dist_c23$t %in% filter), ]
  dist_c13_filtered = dist_c13[!(dist_c13$t %in% filter), ]
  
  write.csv(file = paste0(outputdir,"/csv/",gsub("/","_",filename),"_12distances_filtered.csv"),x=dist_c12_filtered,row.names = F)
  write.csv(file = paste0(outputdir,"/csv/",gsub("/","_",filename),"_23distances_filtered.csv"),x=dist_c23_filtered,row.names = F)
  write.csv(file = paste0(outputdir,"/csv/",gsub("/","_",filename),"_13distances_filtered.csv"),x=dist_c13_filtered,row.names = F)
  
  
  
  
  ##plotting 3D pairwise distances distances
  plot(dist_c12,type="l",ylim=c(0,5), ylab="Distance in µm", xlab="t in min", main="2D pair-wise distances (no filter)", col=3)
  lines(dist_c23,type="l",col=4)
  lines(dist_c13,type="l",col=2)
  legend(x="topright",legend=c("Chic vs Tsix","Tsix vs Linx","Chic vs Linx"),col=c(3,4,2),lty=c(1,1))
  
  plot(dist_c12_filtered,type="l",ylim=c(0,2), ylab="Distance in µm", xlab="t in min", main="2D pair-wise distances (filtering)", col=3)
  lines(dist_c23_filtered,type="l",col=4)
  lines(dist_c13_filtered,type="l",col=2)
  legend(x="topright",legend=c("Chic vs Tsix","Tsix vs Linx","Chic vs Linx"),col=c(3,4,2),lty=c(1,1))
  
  dev.off()
}