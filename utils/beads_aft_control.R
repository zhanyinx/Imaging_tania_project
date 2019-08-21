setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
#input file containing the beads coordinates before/after manually/knime correction
before_correction_file = "/tungstenfs/scratch/ggiorget/_LIVECELL/Analysis_Data/PGK_G8_A11_B3_F6_SingleCells/Uncorrected_traj/beads_images_2h/20190704_uncorrected_beads.csv"
after_correction_file = "/tungstenfs/scratch/ggiorget/_LIVECELL/Analysis_Data/PGK_G8_A11_B3_F6_SingleCells/Corrected_traj/beads_images_2h/20190704_corrected_beads.csv"
manually_corrected_file ="/tungstenfs/scratch/ggiorget/_LIVECELL/Analysis_Data/PGK_G8_A11_B3_F6_SingleCells/manually_corrected_traj/beads_images_2h/20190704_uncorrected_beads_man_corrected.csv"
posColumns = 2:4

#reading beads coordinates
beads_before_correction = read.csv(before_correction_file,header=F)
beads_after_correction = read.csv(after_correction_file,header=F)
beads_after_manual_correction = read.csv(manually_corrected_file,header=F)

#removing the fucking factors
beads_before_correction$V2 = as.numeric(as.character(beads_before_correction$V2))
beads_before_correction$V3 = as.numeric(as.character(beads_before_correction$V3))
beads_before_correction$V4 = as.numeric(as.character(beads_before_correction$V4))

beads_after_correction$V2 = as.numeric(as.character(beads_after_correction$V2))
beads_after_correction$V3 = as.numeric(as.character(beads_after_correction$V3))
beads_after_correction$V4 = as.numeric(as.character(beads_after_correction$V4))

beads_after_manual_correction$V2 = as.numeric(as.character(beads_after_manual_correction$V2))
beads_after_manual_correction$V3 = as.numeric(as.character(beads_after_manual_correction$V3))
beads_after_manual_correction$V4 = as.numeric(as.character(beads_after_manual_correction$V4))


#same bead in three channel will have the same V5 column
beads_before_correction$V5 = gsub("w[0-9]_+", "", x=beads_before_correction$V5) 
beads_after_correction$V5 = gsub("w[0-9]_+", "", x=beads_after_correction$V5) 
beads_after_manual_correction$V5 = gsub("w[0-9]_+", "", x=beads_after_manual_correction$V5) 


#removing fake timepoint 1
beads_before_correction = beads_before_correction[!duplicated(paste0(beads_before_correction$V7,
                                                                     beads_before_correction$V6)),]
beads_after_correction = beads_after_correction[!duplicated(paste0(beads_after_correction$V7,
                                                                   beads_after_correction$V6)),]
beads_after_manual_correction = beads_after_manual_correction[!duplicated(paste0(beads_after_manual_correction$V7,
                                                                                 beads_after_manual_correction$V6)),]


#find_closest 
find_closest = function(beads_w1,beads_w2){
  order12 = rep(-1,nrow(beads_w1))
  mindist12 = rep(-1,nrow(beads_w1))
  for(i in 1:nrow(beads_w1)){
    min=99
    for(j in 1:nrow(beads_w2)){
      if(sqrt(sum((beads_w1[i,]-beads_w2[j,])**2))<min){
        min = sqrt(sum((beads_w1[i,]-beads_w2[j,])**2))
        mindist12[i] = min
        order12[i] = j
      }
    }
  }
  return(list(order12,mindist12))
}

#calculating distances
#extracting channels info (before correction and use it to find the coord)
beads_before_correction_w1 = beads_before_correction[beads_before_correction$V6=="w1_Tracks",posColumns]
beads_before_correction_w2 = beads_before_correction[beads_before_correction$V6=="w2_Tracks",posColumns]
beads_before_correction_w3 = beads_before_correction[beads_before_correction$V6=="w3_Tracks",posColumns]

#extract order (find the closest for each spot)
order_mindist_12 = find_closest(beads_before_correction_w1,beads_before_correction_w2)
order_mindist_23 = find_closest(beads_before_correction_w2,beads_before_correction_w3)
order_mindist_13 = find_closest(beads_before_correction_w1,beads_before_correction_w3)

#extracting channels info
beads_after_correction_w1 = beads_after_correction[beads_after_correction$V6=="w1_Tracks",posColumns]
beads_after_correction_w2 = beads_after_correction[beads_after_correction$V6=="w2_Tracks",posColumns]
beads_after_correction_w3 = beads_after_correction[beads_after_correction$V6=="w3_Tracks",posColumns]


beads_after_manual_correction_w1 = beads_after_manual_correction[beads_after_manual_correction$V6=="w1_Tracks",posColumns]
beads_after_manual_correction_w2 = beads_after_manual_correction[beads_after_manual_correction$V6=="w2_Tracks",posColumns]
beads_after_manual_correction_w3 = beads_after_manual_correction[beads_after_manual_correction$V6=="w3_Tracks",posColumns]



##calculating the distances
dist12_before = beads_before_correction_w1 - beads_before_correction_w2[order_mindist_12[[1]],]
dist23_before = beads_before_correction_w2 - beads_before_correction_w3[order_mindist_23[[1]],]
dist13_before = beads_before_correction_w1 - beads_before_correction_w3[order_mindist_13[[1]],]

dist12_after = beads_after_correction_w1 - beads_after_correction_w2[order_mindist_12[[1]],]
dist23_after = beads_after_correction_w2 - beads_after_correction_w3[order_mindist_23[[1]],]
dist13_after = beads_after_correction_w1 - beads_after_correction_w3[order_mindist_13[[1]],]

dist12_after_manual = beads_after_manual_correction_w1 - beads_after_manual_correction_w2[order_mindist_12[[1]],]
dist23_after_manual = beads_after_manual_correction_w2 - beads_after_manual_correction_w3[order_mindist_23[[1]],]
dist13_after_manual = beads_after_manual_correction_w1 - beads_after_manual_correction_w3[order_mindist_13[[1]],]

#plotting
making_df_for_ggplot = function(before,after,after_manual,name){
  before = as_tibble(before)
  before = gather(before)
  before = data.frame(gather(before))
  before$type = "before"
  
  after = as_tibble(after)
  after = gather(after)
  after = data.frame(gather(after))
  after$type = "after"
  
  after_manual = as_tibble(after_manual)
  after_manual = gather(after_manual)
  after_manual = data.frame(gather(after_manual))
  after_manual$type = "after_manual"
  
  df = rbind(before,after,after_manual)
  df[,1][df[,1]=="V2"]="x"
  df[,1][df[,1]=="V3"]="y"
  df[,1][df[,1]=="V4"]="z"
  
  names(df)=c("coord","distance","type")
  return(df)
  
}

ggdist12 = making_df_for_ggplot(dist12_before,dist12_after,dist12_after_manual,"12")
ggplot(ggdist12,aes(x=coord,y=distance,colour=type)) + geom_boxplot() +
  ylim(-0.5,0.5) + ggtitle("c1 vs c2")

ggdist23 = making_df_for_ggplot(dist23_before,dist23_after,dist23_after_manual,"23")
ggplot(ggdist23,aes(x=coord,y=distance,colour=type)) + geom_boxplot() +
  ylim(-0.5,0.5) + ggtitle("c2 vs c3")

ggdist13 = making_df_for_ggplot(dist13_before,dist13_after,dist13_after_manual,"13")
ggplot(ggdist13,aes(x=coord,y=distance,colour=type)) + geom_boxplot() +
  ylim(-0.5,0.5) + ggtitle("c1 vs c3")
