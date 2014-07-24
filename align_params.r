setwd('~/deming_lab/psychrophilic_crude_genes/round_2')

cluster <- 'p450_cluster_0'

clusters <- c('p450_cluster_0', 'p450_cluster_1', 'FA_desaturase_cluster_0', 'FA_desaturase_cluster_1', 'Bac_luciferase_cluster_0', 'Pyr_redox_3_cluster_0')

cluster_points <- matrix(ncol = 3, nrow = 0)

for(cluster in clusters){

  aligned_params <- read.table(paste0(cluster, '_param_alignment.txt'), sep = ',')
  ap_cold <- as.matrix(aligned_params[which(aligned_params[,3] == 'cold'),6:length(colnames(aligned_params))])
  ap_control <- as.matrix(aligned_params[which(aligned_params[,3] == 'control'),6:length(colnames(aligned_params))])
  
  ap_cold_mean <- apply(ap_cold, 2, mean, na.rm = T)
  ap_control_mean <- apply(ap_control, 2, mean, na.rm = T)
  ap_delta_mean <- ap_cold_mean - ap_control_mean
  
  ap_cold_sd <- apply(ap_cold, 2, sd, na.rm = T)
  ap_control_sd <- apply(ap_control, 2, sd, na.rm = T)
  
  x <- (6:length(colnames(aligned_params)))
  
  pdf(paste0(cluster, '_param_align.pdf'),
      width = 6,
      height = 4)
  
  plot(ap_delta_mean ~ x,
       type = 'l',
       col = 'blue',
       ylab = 'Flexibility difference',
       xlab = 'Position',
       ylim = c(min(ap_delta_mean, na.rm = T), max(ap_control_sd + ap_cold_sd, na.rm = T)),
       xlim = c(min(which(is.na(ap_delta_mean) == F)), max(which(is.na(ap_delta_mean) == F))),
       main = cluster)
      
  
  points((ap_control_sd + ap_cold_sd) ~ x,
       type = 'l',
       col = 'black')
  
  lines(c(0,1000),
        c(0,0))
  
  ## draw lines where difference in means exceeds the sum of the standard deviations
  
  for(i in which(is.na(ap_delta_mean) == F)){
    try(
      if(abs(ap_delta_mean[i]) > (ap_cold_sd[i] + ap_control_sd[i])){
        print(ap_delta_mean[i])
        lines(c(i, i), c(-1, 1),
              lty = 2,
              col = 'orange')
        out <- c(cluster, i, ap_delta_mean[i])
        cluster_points <- rbind(cluster_points, out)
      }, silent = T)
  }
  
  legend('topleft',
         legend = c('Mean', 'SD sum'),
         col = c('blue', 'black'),
         lty = 2,
         bg = 'white')
  
  dev.off()
}


write.table(cluster_points, 'cluster_flex_consensus_points.txt', sep = '\t', quote = F, col.names = F, row.names = F)

#### raw points ####

for(i in 1:12){
  if(aligned_params[i, 3] == 'cold'){
    points(ap_matrix[i,] ~ x,
         type = 'l',
         col = 'blue')
  }else{points(ap_matrix[i,] ~ x,
             type = 'l',
             col = 'red')
  }
}
