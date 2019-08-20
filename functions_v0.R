setwd("/tungstenfs/scratch/ggiorget/zhan/2019/190627_Pia_MSD")
library(matrixStats) #for rowProds




#Make the auxiliary matrix from the chosen dataset using the sel_track as x and sel_ts as y and dataset=data
aux.matrix <- function(data,x,y){
  a=data[,trackColumn]
  b=data[,timeColumn]
  if(nrow(data[a==x & b==y,])==1){return(1)} else {return(0)}
}

#MSD for matrices where each row is a position in N dimension where N is the number of the columns of the matrices
msd <- function(mat1,mat2){
  return(rowSums((mat1-mat2)*(mat1-mat2)))
}

#centre of mass function
cm <- function(mat1){
  return(colMeans(mat1))
}


##distance between 2 channels matching time (time column is the first column)
dist_channels = function(channel1,channel2){
  ncol1 = ncol(channel1)
  ncol2 = ncol(channel2)
  if(ncol1 != ncol2) stop("Column dimensions are different in dist_channels")
  c12 = merge(channel1,channel2,by=1)
  dist=NULL
  for(i in 1:nrow(c12)){
    dist=c(dist,dist(t(matrix(c(c12[i,2:(ncol1)],c12[i,(ncol1+1):(2*ncol1-1)]),ncol=2))))
  }
  return(data.frame(t=c12[,1],dist=dist))
}




#Calculates the pairwise minimum distance between all tracks at a certain timepoint for two channels and return the tracks with that have the minimum distance
pair_dis <- function(data1, data2){
  pair=NULL
  for (k in 1:nrow(data1)) {
    dist=NULL
    for (p in 1:nrow(data2)) {
      temp <- dist(rbind(data1[k,posColumns], data2[p,posColumns]))
      dist <- rbind(dist, temp)
    }
    pair <- cbind(pair, dist)
  }
  colnames(pair) <- paste(data1$channel, "_", data1$TRACK_ID, sep = "")
  rownames(pair) <- paste(data2$channel, "_", data2$TRACK_ID, sep = "")
  keep <- which(pair == min(pair), arr.ind = TRUE)
  paste(data1[keep[2],6], data1[keep[2],7], sep = "_") -> data1
  paste(data2[keep[1],6], data2[keep[1],7], sep = "_") -> data2
  #return(pair)
  
  return(c(data1,data2))
}





findmaxtrack <- function(column){
  maxi <- max(table(column))
  pos <- which(table(column)==maxi, arr.ind = TRUE)
  
  return(c(maxi, rownames(pos)))
}





##given a vector containing the heights of a histogram, find the maximum area of the rectangle in the histogram
max_area_histogram = function(histogram,m_area,s,e,h){ 
  # This function calulates maximum  
  # rectangular area under given  
  # histogram with n bars 
  
  # Create an empty stack. The stack  
  # holds indexes of histogram[] list.  
  # The bars stored in the stack are 
  # always in increasing order of  
  # their heights. 
  stack = NULL
  if(is.null(m_area)){
    max_area = 0 # Initalize max area 
  }else{
    max_area = m_area
  }
  start = s
  end = e
  height = h
  # Run through all bars of 
  # given histogram 
  index = 1
  while(index <= length(histogram)){ 
  
  # If this bar is higher  
  # than the bar on top 
  # stack, push it to stack 
  
    if(length(stack)==0 || (histogram[stack[length(stack)]] <= histogram[index])){ 
      stack = c(stack,index) 
      index = index + 1
  
      # If this bar is lower than top of stack, 
      # then calculate area of rectangle with  
      # stack top as the smallest (or minimum 
      # height) bar.'i' is 'right index' for  
      # the top and element before top in stack 
      # is 'left index' 
    }else{ 
      # pop the top 
      top_of_stack = stack[length(stack)] 
      stack = stack[-length(stack)]
      # Calculate the area with  
      # histogram[top_of_stack] stack 
      # as smallest bar 
      if(length(stack)==0){
        area=histogram[top_of_stack] * (index-1)
      }else{
        area = histogram[top_of_stack] * ((index - stack[length(stack)] - 1 )) 
      }  
      # update max area, if needed 
      if(area>max_area){
        max_area = area
        end=index -1
        if(length(stack)==0){
          start=1
        }else{ 
          start=stack[length(stack)] +1
        }
        height = histogram[top_of_stack]
      }
    }
  }  

# Now pop the remaining bars from  
# stack and calculate area with  
# every popped bar as the smallest bar 
  while(length(stack)>0){ 
  
    # pop the top 
    top_of_stack = stack[length(stack)] 
    stack = stack[-length(stack)]

    # Calculate the area with  
    # histogram[top_of_stack]  
    # stack as smallest bar 
    if(length(stack)==0){
      area=histogram[top_of_stack] * (index-1)
    }else{
      area = histogram[top_of_stack] * ((index - stack[length(stack)] - 1)) 
    }

    # update max area, if needed 
    if(area>max_area){
      max_area = area
      end=index -1
      if(length(stack)==0){start=1}else{ start=stack[length(stack)] + 1 }
      height = histogram[top_of_stack]
    }
  }
# Return maximum area under  
# the given histogram 
  return(c(max_area,start,end,height))
}


##extract biggest submatrix of 1s in a boolean matrix (needs max_area_histogram function)
submatrix = function(m){
  v = m[1,]
  max_a1 = NULL
  start1 = NULL
  end1 = NULL
  height1 = NULL
  select_row = NULL
  for(i in 2:nrow(m)){
    v = m[i,]*(m[i,]+v)
    res = max_area_histogram(v,max_a1,start1,end1,height1)
    if(is.null(max_a1)){
      max_a1=res[1]
      select_row = i  
      start1 = res[2]
      end1 = res[3]
      height1= res[4]
    }else if(res[1]>max_a1){
      max_a1=res[1]
      select_row = i  
      start1 = res[2]
      end1 = res[3]
      height1= res[4]
    }
  }
  return(c(res,select_row))
}


##function to reorder matrix such that we have the biggest rectangular submatrix of 1 at the top left corner
reorder_matrix = function(m){
  position_zeros_r = apply(m,1,
                         function(x){
                           l = which(x==0); 
                           if(length(l)==0){
                             return(length(x)+1)
                           }else{
                             return(l[1])
                           }
                         }
  )
  m1 = m[order(position_zeros_r,decreasing=T),]
  position_zeros_c = apply(m1,2,
                           function(x){
                             l = which(x==0); 
                             if(length(l)==0){
                               return(length(x)+1)
                             }else{
                               return(l[1])
                             }
                           }
  )
  return(m1[,order(position_zeros_c,decreasing=T)])
  #return(m1)
  }


# 
#     a = sample(c(1,0),prob = c(0.8,0.2),replace = T,100)
#     m=matrix(a,ncol=10)
#     reorder_matrix(m)
#     res=submatrix(reorder_matrix(m))
#    m = (reorder_matrix(m))
#    m  
#    subset = m[((res[5]-res[4]+1):res[5]),res[2]:res[3]]
# # # subset
# # length(subset) == res[1]
# order = NULL
# for(i in 1:length(a)){
#   order= c(order,a[[i]][1])
# }

