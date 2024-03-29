setwd(dirname(rstudioapi::getSourceEditorContext()$path))

####global options
####CHOOSE DATASET HERE BY CHOOSING FILEPATH:
dirpath = "/tungstenfs/scratch/ggiorget/_LIVECELL/Analysis_Data/PGK_G8_A11_B3_F6_SingleCells/Uncorrected_traj/2h_movies/" #set the filepath to files you are intested in
af_dirpath = "/tungstenfs/scratch/ggiorget/_LIVECELL/Analysis_Data/PGK_G8_A11_B3_F6_SingleCells/manually_corrected_traj/af_matrices_2h/"
##choose output directory for saving corrected fiels
outputdir = "/tungstenfs/scratch/ggiorget/_LIVECELL/Analysis_Data/PGK_G8_A11_B3_F6_SingleCells/manually_corrected_traj/2h_movies/" 
posColumns = 2:4


####functions
aff_trans = function(m, v, t){
  x = m[1,] %*% v
  y = m[2,] %*% v
  z = m[3,] %*% v
  v_r = c(x,y,z) - t
  return(as.vector(v_r))
}


###read in data and create output directories 
directories = paste(outputdir, list.files(dirpath), sep="")
lapply(directories, dir.create)

afs = sub("_AF.csv", "", list.files(af_dirpath))
filelist = list.files(path=dirpath, pattern="\\uncorrected.csv$",recursive = TRUE)

for (i in 1:length(afs)) {
  af_matrix = read.csv(paste(af_dirpath, afs[i], "_AF.csv", sep = ""))
  afm2 = af_matrix[af_matrix$channels=="W1-W2",]
  afm3 = af_matrix[af_matrix$channels=="W1-W3",]
  
  keep = grepl(afs[i], filelist)
  files = filelist[keep]

  for (q in 1:length(files)) {
    traj = read.csv(paste(dirpath,files[q], sep=""), sep = ",", header = TRUE)
    corrected = traj
    #clean up some of the dataset columns
    traj$RowID <- gsub("_w[0-9]_[0-9]+_[0-9]+", "", x=traj$RowID)
    traj$Original.RowID  = gsub("cell[0-9]_w", "C", x=as.character(traj$Original.RowID))
    names(traj)[names(traj) == "Original.RowID"] <- "channel"
    names(traj)[names(traj) == "particle"] <- "TRACK_ID"
    names(traj)[names(traj) == "RowID"] <- "Cell_no"
    
    c2 = corrected[traj$channel=="C2",]
    c3 = corrected[traj$channel=="C3",]
    
    for (p in 1:nrow(c2)) {
      c2[p,posColumns] = aff_trans(as.matrix(afm2[,1:3]), as.numeric(c2[p,posColumns]), as.matrix(afm2[,4]))
    }

    for (r in 1:nrow(c3)) {
      c3[r,posColumns] = aff_trans(as.matrix(afm3[,1:3]), as.numeric(c3[r,posColumns]), as.matrix(afm3[,4]))
    }
    
    corrected = rbind(corrected[traj$channel=="C1",], c2, c3)
    filename = sub(".csv", "", files[q])
    write.csv(x=corrected, file = paste(outputdir,filename, "_man_corrected.csv", sep = ""), row.names = FALSE)
  }

}


