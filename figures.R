library(cowplot)

calculate_place_cell_pct_by_ind_figures <- function(ind, 
                                                    spike_train, 
                                                    stim_trace){
  
  sliced_spike_train <- spike_train[,ind]
  sliced_stim_trace <- stim_trace[ind]
  
  firing_rate <- rowMeans(sliced_spike_train) / dt

  tmp <- compute_tuning(sliced_spike_train, 
                        sliced_stim_trace)
  stim_prob <- tmp[[1]]
  tuning_curve <- tmp[[2]]
  rm(tmp)
  
  sliced_SI <- compute_SI(stim_prob, 
                          tuning_curve, 
                          firing_rate)
  
  
  random_pval <- 
               lapply(c(1:nrow(sliced_spike_train)), 
                      function(i, var2) { 
                        compute_place_signif(i,
                                             sliced_spike_train,
                                             sliced_stim_trace,
                                             firing_rate,
                                             sliced_SI,
                                             shuffle_type=var2,
                                             num_shuffles = 500,
                                             return_pval = T,
                                             rate_to_use=2)},
                      var2="r")
  
  
  cyclic_pval <- 
               lapply(c(1:nrow(sliced_spike_train)), 
                      function(i, var2) { 
                        compute_place_signif(i,
                                             sliced_spike_train,
                                             sliced_stim_trace,
                                             firing_rate,
                                             sliced_SI,
                                             shuffle_type=var2,
                                             num_shuffles = 500,
                                             return_pval = T,
                                             rate_to_use=2,
                                             verbose=F)},
                      var2="c")
  
  res = list(SI=sliced_SI,
             random=random_pval,
             cyclic=cyclic_pval)
  
  return(res)
  
}

generate_sampled_ind_figure <- function(duration,
                                        length) {

    return(matrix(c(1:length)[1:(duration / dt)], nrow=1))
}

pct_by_sample_size_figures <- function(spike_train, stim_trace) {
  
  pct_by_duration <- c()
  dur <- seq(10,ncol(spike_train) * dt, length.out=20)
  
  result_rand <- list()
  result_cyclic <- list()
  for (i in 1:len(dur)) {
    duration <- dur[i]
    ind_mat <- generate_sampled_ind_figure(duration,  ncol(spike_train))
    random_mat <- c()
    cyclic_mat <- c()
    
    for (i in 1:nrow(ind_mat)) {
      print("############### DURATION, ind")
      print(duration)
      print(i)
      res = calculate_place_cell_pct_by_ind_figures(ind_mat[i,],
                                                    spike_train,
                                                    stim_trace)
      random_mat <- cbind(random_mat, unlist(res$random))
      cyclic_mat <- cbind(cyclic_mat, unlist(res$cyclic))
      
    }
    
    result_rand <- append(result_rand, list(random_mat))
    result_cyclic <- append(result_cyclic, list(cyclic_mat))
    
  }
  
  names(result_rand) <- dur
  names(result_cyclic) <- dur

  pval_rand <- 
    do.call(cbind, lapply(result_rand, function(mt) {apply(mt, 1, mean)}))
  
  pval_cyclic <- 
    do.call(cbind, lapply(result_cyclic, function(mt) {apply(mt, 1, mean)}))
 
  return(list(cyclic=pval_cyclic, rand=pval_rand, duration=dur))
   
}




plot_pval_by_sample_duration_figure <- 
  function(working_path_list, 
           ses_ind=c(5:8, 13:16), 
           ext="",
           fit_sim="",
           bwidth=0.75,
           verbose=F,
           type1=F) {
    df_list <- list()
    
    # Read all dataframes with computed place cell fraction by sample duration
    for (working_path in working_path_list) {
      
      session_paths <- 
        list.dirs(working_path, recursive=F)[grep("session", 
                                                  list.dirs(working_path, 
                                                            recursive=F))]
      # Run through all paths
      for (tpath in session_paths) {
        
        tpath3 <- str_replace(tpath, "/", "\\\\")
        tpath2 <- sprintf("%s\\%s%s%s.Rda",
                          tpath3,
                          fit_sim,
                          "cyclic_pval_by_dur",
                          ext)
        tpath4 <- sprintf("%s\\%s", tpath3, "duration_pval_matrices\\")
        
        if (len((grep(sprintf("cyclic_pval_by_dur%s.Rda", ext), list.files(tpath)))) == 0) {
          next
        }
        
        idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
        if (sum(idx == ses_ind) == 0){
          next
        }
        
        if (verbose) {
          print(sprintf("Loading path %s", tpath2))
        }
        
        load(tpath2)
        
        df_list <- append(df_list, list(cyclic_list))
      }
    }
    
    
    # # Properly bin all dataframes
    # df_list2 <- lapply(c(1:len(df_list)), 
    #                    function(i) {
    #                      df <- df_list[[i]]
    #                      st_size <-  df[[7]]
    #                      df <- df[-7]
    #                      
    #                      bins <- seq(10, st_size * dt, length.out=6)
    #                      bins <- bins / max(bins) * 20
    #                      
    #                      all_cells <- do.call(cbind, lapply(df, rowMeans))
    #                      all_cells_n <- apply(all_cells, 2, function(c) {(c + 10^-30) / (all_cells[,ncol(all_cells)] + 10^-30)})
    #                      all_cells_n_2 <- apply(all_cells_n, 2, function(c) {(c / rowSums(all_cells_n))})
    #                      f_vec <- ifelse(type1, list(all_cells_n), list(all_cells_n_2))[[1]]
    #                      
    #                      if (!type1){
    #                        f_vec  <- f_vec  - f_vec[,6]
    #                      }
    #                      
    #                      f_vec <- colMeans(f_vec)
    #                      #f_vec <- colMeans(abs(log2(all_cells + 10^-30) - log2(all_cells[,6] + 10^-30)))
    #                      
    #                      if (type1) {
    #                        f_vec <- f_vec / sum(f_vec)
    #                      }
    #                      
    #                      names(f_vec) <- bins
    #                      return(f_vec) 
    #                    })
    # 
    
    # Properly bin all dataframes
    df_list2 <- lapply(c(1:len(df_list)),
                       function(i) {
                         df <- df_list[[i]]
                         st_size <-  df[[7]]
                         df <- df[-7]

                         bins <- seq(10, st_size * dt, length.out=6)
                         bins <- bins / max(bins) * 20

                         all_cells <- do.call(cbind, lapply(df, rowMeans))

                         f_vec <- t(apply(all_cells, 1, function(r) {r / sum(r)}))
                         #f_vec <- all_cells
                         f_vec <- f_vec[!is.nan(f_vec[,6]),]
                         f_vec <- colMeans(f_vec)
                         #f_vec <- (f_vec - min(f_vec)) / (max(f_vec) - min(f_vec))

                         names(f_vec) <- bins
                        
                         return(f_vec)
                       })

    
    val_mat <- do.call(rbind, df_list2)
    bin_mat <- do.call(rbind, lapply(df_list2, function(r) {as.numeric(names(r))}))
    
    #val_mat <- val_mat[val_mat[,ncol(val_mat)] != 0,]
    plot_df <- data.frame(pvd=colMeans(val_mat),
                         pvdsd=apply(val_mat, 2, sd) * 1,
                         time=as.numeric(colnames(val_mat)))
    
    g <-
      ggplot(plot_df, aes(x=time, y=pvd)) + 
      geom_ribbon(aes(ymin=pvd-pvdsd, ymax=pvd+pvdsd), color="NA", alpha=0.2, fill="darkmagenta", size=5) +
      geom_line(size=1.2, color="darkmagenta") + 
      theme_light() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position="NA",
            panel.border = element_blank(),
            panel.background = element_blank()) +
      labs(x="Sample duration (minutes)",
           y = TeX(sprintf("%s fold change (AU)", r'($P_{value}$)'))) #y = TeX(r'(Normalized  \Delta $P_{value}$)'))
    
    
    g2 <-
      ggplot(plot_df, aes(x=time, y=pvd)) + 
      geom_errorbar(aes(ymin=pvd-pvdsd, ymax=pvd+pvdsd), color="darkmagenta", size=1, width=bwidth) +
      geom_line(size=1.2, color="darkmagenta") + 
      theme_light() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position="NA",
            panel.border = element_blank(),
            panel.background = element_blank()) +
      labs(x="Sample duration (minutes)",
           y = TeX(sprintf("%s fold change (AU)", r'($P_{value}$)'))) #y = TeX(r'(Normalized  \Delta $P_{value}$)'))
    #geom_point(size=2, color="gray50")
    
   
    
    return(list(g, plot_df, val_mat, bin_mat, g2))
  }



plot_place_cell_pct_by_sample_duration_figure <- 
                                          function(working_path_list, 
                                                   ses_ind=c(1:16), 
                                                   bin_size=2.5,
                                                   cyclic_only=F,
                                                   ext="",
                                                   prefix="",
                                                   fit_sim="",
                                                   verbose=F,
                                                   bwidth=0.75,
                                                   average_across_individuals=T) {
  df_list <- list()
  
  # Read all dataframes with computed place cell fraction by sample duration
  for (working_path in working_path_list) {
    
    session_paths <- 
      list.dirs(working_path, recursive=F)[grep("session", 
                                                list.dirs(working_path, 
                                                          recursive=F))]
    # Run through all paths
    for (tpath in session_paths) {
      
      
      
      tpath3 <- str_replace(tpath, "/", "\\\\")
      tpath_f <- sprintf("%s\\%s", tpath, fit_sim)
      tpath2 <- sprintf("%s%s%s%s.Rda",
                        tpath_f,
                        prefix,
                        "place_cells_by_sample_duration",
                        ext)
      tpath4 <- sprintf("%s\\%s", tpath3, "duration_pval_matrices\\")
      
      if (len((grep(sprintf("%splace_cells_by_sample_duration%s.Rda", prefix, ext), list.files(tpath_f)))) == 0) {
        print(sprintf("Missing: %s", tpath_f))
        next
      }
      
      idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
      if (sum(idx == ses_ind) == 0){
        next
      }

      if (verbose) {
        print(sprintf("Loading path %s", tpath2))
      }
      
      load(tpath2)

      df_list <- append(df_list, list(tmp))
    }
  }

  
  # # Properly bin all dataframes
  # df_list2 <- lapply(c(1:len(df_list)), 
  #                    function(i) {
  #                      df <- df_list[[i]]
  #                      dt = 0.05
  #                      frame_factor <- ((20 * 60 / dt) / (df)$total_frames)[[1]]
  #                      
  #                      a <- (as.numeric(rownames(df)) /  60) * frame_factor
  #                      
  #                      if (max(a) > 21) {
  #                        print("DT = 0.1!!!")
  #                        dt = 0.1
  #                        frame_factor <- ((20 * 60 / dt) / (df)$total_frames)[[1]]
  #                        
  #                        a <- (as.numeric(rownames(df)) /  60) * frame_factor
  #                        
  #                      }
  #                      
  #                      
  #                      new_df <- data.frame(bins=bin(a, bin_size),
  #                                           rand=df[,1],
  #                                           cyclic=df[,2])
  #                      new_df <- new_df[!duplicated(new_df$bin),]
  #                      
  #                      bins <- bin(a,bin_size)
  #                      bins <- bins[!duplicated(bins)]
  #                      
  #                      rownames(new_df) <- bins
  #                      return(new_df) 
  #                    })
  
  # Properly bin all dataframes
  df_list2 <- lapply(c(1:len(df_list)), 
                     function(i) {
                       
                       df <- df_list[[i]]
                       new_df <- data.frame(bins=seq(1, 20, length.out = 8),
                                            rand=df[,1],
                                            cyclic=df[,2])
                       
                       rownames(new_df) <- seq(1, 20, length.out = 8)
                       return(new_df) 
                     })

  
  merged_rand = df_list2[[1]][,c("bins", "rand")]
  merged_cyclic = df_list2[[1]][,c("bins", "cyclic")]
  
  
  
  # Merge dataframes
  for (i in 2:len(df_list2)) {
    merged_rand <- merge(merged_rand,
                         df_list2[[i]][,c("bins", "rand")],
                         by="bins",
                         all=T)
    
    merged_cyclic <- merge(merged_cyclic,
                           df_list2[[i]][,c("bins", "cyclic")],
                           by="bins",
                           all=T)        
  }
  
  
  
  if (average_across_individuals) {
    num_sessions <- len(ses_ind) * 2
    merged_cyclic_tmp <- as.matrix(merged_cyclic[,1], ncol=1)
    
    for (i in 1:(len(df_list2) / num_sessions)) {
        ind <- c(((i-1) * num_sessions + 1):(i * num_sessions)) + 1
        merged_cyclic_tmp <- cbind(merged_cyclic_tmp,
                                   rowMeans(merged_cyclic[,ind]))
    }
  
    colnames(merged_cyclic_tmp)  <- c("bins", 1:(len(df_list2) / num_sessions))  
    merged_cyclic <- merged_cyclic_tmp
  }
  
  
  
  
  
  means_rand <- rowMeans(merged_rand[,-1], na.rm=T)
  sem_rand <- apply(merged_rand[,-1], 1, sd, na.rm=T)
  #labels_rand <- sapply(merged_rand[,1], function(b) {as.numeric(str_split(str_split(b, "\\(")[[1]][2], ",")[[1]][1])})
  labels_rand <- merged_rand[,1]
  
  means_cyc <- rowMeans(merged_cyclic[,-1], na.rm=T)
  sem_cyc <- apply(merged_cyclic[,-1], 1, sd, na.rm=T)
  #labels_cyc <- sapply( merged_cyclic[,1], function(b) {as.numeric(str_split(str_split(b, "\\(")[[1]][2], ",")[[1]][1])})
  labels_cyc <- merged_cyclic[,1]
  
  
  df_r <- data.frame(fraction=means_rand,
                     time=as.numeric(labels_rand),
                     sd=sem_rand)
  
  df_r$Shuffle <- rep("Random shuffle", nrow(df_r))
  
  df_c <- data.frame(fraction=means_cyc,
                     time=as.numeric(labels_cyc),
                     sd=sem_cyc)
  
  df_c$Shuffle <- rep("Cyclic shuffle", nrow(df_c))
  
  if (cyclic_only) {
    plot_df <- df_c
  } else {
    plot_df <- rbind(df_r, df_c)  
  }
  
  g <- ggplot(plot_df, aes(x=time, y=fraction, group=Shuffle, color=Shuffle)) + 
    geom_linerange(aes(ymin=fraction-sd, ymax=fraction+sd, group=Shuffle, color=Shuffle), alpha=0.8, size=1) +
    geom_line(size=1.2) +
    scale_color_manual(values = c("darkmagenta", "lightsalmon2")) +
    scale_fill_manual(values = c("darkmagenta", "lightsalmon2")) +
    theme_light() +
    ylim(c(0, 0.9))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position="NA",
          panel.border = element_blank(),
          panel.background = element_blank()) +
    labs(x="Sample duration (minutes)", 
         y = "Fraction of place cells (%)")
  
  gn <- ggplot(plot_df, aes(x=time, y=fraction, group=Shuffle, color=Shuffle)) + 
    geom_ribbon(aes(ymin=fraction-sd, ymax=fraction+sd, group=Shuffle, fill=Shuffle), color="NA", alpha=0.1, size=1, width=bwidth) +
    geom_line(size=1.2) +
    scale_color_manual(values = c("darkmagenta", "lightsalmon2")) +
    scale_fill_manual(values = c("darkmagenta", "lightsalmon2")) +
    theme_light() +
    ylim(c(0, 0.9))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position="NA",
          panel.border = element_blank(),
          panel.background = element_blank()) +
    labs(x="Sample duration (minutes)", 
         y = "Fraction of place cells (%)")
  
  gn2 <- ggplot(plot_df, aes(x=time, y=fraction, group=Shuffle, color=Shuffle)) + 
    geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd, group=Shuffle, fill=Shuffle), width=bwidth, size=1) +
    geom_line(size=1.2) +
    scale_color_manual(values = c("darkmagenta", "lightsalmon2")) +
    scale_fill_manual(values = c("darkmagenta", "lightsalmon2")) +
    theme_light() +
    ylim(c(0, 0.9))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position="NA",
          panel.border = element_blank(),
          panel.background = element_blank()) +
    labs(x="Sample duration (minutes)", 
         y = "Fraction of place cells (%)")
    #geom_point(size=2, color="gray50")
  
  
    g2 <- ggplot(data.frame(Fraction=means_rand, Time=labels_rand)) +
          geom_line(aes(y=Fraction, x=Time), col="lightsalmon2", size=1.6)
    
    
    for (i in 2:ncol(merged_rand)) {
      sub_ind <- !is.na(merged_rand[,i])
      tmp_df <- data.frame(Fraction=merged_rand[sub_ind,i],  Time=labels_rand[sub_ind])
      
      g2 <- g2 + 
            geom_line(data=tmp_df,aes(x=Time, y=Fraction), col="lightsalmon2", alpha=0.2, size=0.6)
    }
    
    cyc_df <- data.frame(Fraction=means_cyc, Time=labels_cyc)
    g2 <- g2 + geom_line(data=cyc_df, aes(y=Fraction, x=Time), col="darkmagenta", size=1.6)
    
    
    for (i in 2:ncol(merged_cyclic)) {
      sub_ind <- !is.na(merged_cyclic[,i])
      tmp_df <- data.frame(Fraction=merged_cyclic[sub_ind,i], Time=labels_cyc[sub_ind])
      
      g2 <- g2 + geom_line(data=tmp_df,
                           aes(x=Time, y=Fraction), col="darkmagenta", alpha=0.2, size=0.6)
    }
    
    g3 <- g2 +
      labs(x="Sample duration (minutes)", 
           y = "Fraction of place cells (%)") + 
      theme_light() +
      ylim(c(0, 0.9)) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position="NA",
            panel.border = element_blank(),
            panel.background = element_blank()) 
    
    
    sub_ind <- !is.na(merged_cyclic[,2])
    tmp_df <- data.frame(Fraction=merged_cyclic[sub_ind,2], Time=labels_cyc[sub_ind])
    
    g4 <- ggplot(tmp_df, aes(x=Time, y=Fraction)) + geom_line(col="darkmagenta", alpha=0.2, size=0.6)
    
    
    for (i in 3:ncol(merged_cyclic)) {
      sub_ind <- !is.na(merged_cyclic[,i])
      tmp_df <- data.frame(Fraction=merged_cyclic[sub_ind,i], Time=labels_cyc[sub_ind])
      
      g4 <- g4 + geom_line(data=tmp_df,
                           aes(x=Time, y=Fraction), col="darkmagenta", alpha=0.2, size=0.6)
    }
    
    if (!cyclic_only) {
      sub_ind <- !is.na(merged_rand[,2])
      tmp_df <- data.frame(Fraction=merged_rand[sub_ind,2],  Time=labels_rand[sub_ind])
      
      g4 <- g4 + 
        geom_line(data=tmp_df,aes(x=Time, y=Fraction), col="lightsalmon2", alpha=0.2, size=0.6)
      
      for (i in 3:ncol(merged_rand)) {
        sub_ind <- !is.na(merged_rand[,i])
        tmp_df <- data.frame(Fraction=merged_rand[sub_ind,i],  Time=labels_rand[sub_ind])
        
        g4 <- g4 + 
          geom_line(data=tmp_df,aes(x=Time, y=Fraction), col="lightsalmon2", alpha=0.2, size=0.6)
      }
    }
    
    g5 <-  g4 +
      geom_line(data=plot_df, aes(x=time, y=fraction, group=Shuffle, color=Shuffle), size=1.6) + 
      geom_linerange(data=plot_df, aes(x=time, y=fraction, ymin=fraction-sd, ymax=fraction+sd, group=Shuffle, color=Shuffle), alpha=0.8, size=1.6) +
      scale_color_manual(values = c("darkmagenta", "lightsalmon2")) +
      scale_fill_manual(values = c("darkmagenta", "lightsalmon2")) +
      theme_light() +
      ylim(c(0, 0.9))+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position="NA",
            panel.border = element_blank(),
            panel.background = element_blank())  +
      labs(x="Sample duration (minutes)", 
           y = "Fraction of place cells (%)")
    

    
    
    return(list(gn, g3, g5, plot_df, merged_cyclic, gn2))
}


create_pct_by_activity_df <- 
  function(path_list,indices=c(5:8, 13:16), activity_bin_size=15, critical_pval=0.05, fit_sim="") {
    cells_pvalue_df <- list()
    
    subfields <- c()
    for (working_path in path_list) {
      
      # extract all paths
      session_paths <- 
        list.dirs(working_path, recursive=F)[grep("session", list.dirs(working_path, recursive=F))]
      
      # Run through all paths
      for (tpath in session_paths) {
        
        idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
        if (sum(idx == indices) == 0){
          next
        }
        
        if (len((grep("tuning_df", list.files(tpath)))) == 0) {
          next
        }
        
        
        
        load(sprintf("%s/%s%s", tpath, fit_sim, "activity_tuning_df_cor.Rda"))
        cells_pvalue_df <- append(cells_pvalue_df, list(tuning_df))
        subfields <- c(subfields, 
                       ifelse(grepl("CA1", tpath), "CA1", "CA3"))
      }
      
    }
    
    pct_df <- 
      lapply(cells_pvalue_df,
             function(sub_df) {
               #binned <- bin(sub_df$active_bins, activity_bin_size)
               binned <- cut(sub_df$act, seq(0,150, by=10))
               random_pct <- lapply(1:len(levels(binned)),
                                    function(b) {
                                      ind <- which(as.numeric(binned) == b)
                                      sum(sub_df[ind,]$random_significance < critical_pval) / len(ind)
                                    })
               
               cyclic_pct <- lapply(1:len(levels(binned)),
                                    function(b) {
                                      ind <- which(as.numeric(binned) == b)
                                      sum(sub_df[ind,]$cyclic_significance < critical_pval) / len(ind)
                                    })
               
               return(data.frame(bins=levels(binned)[!is.nan(unlist(random_pct))],
                                 random=unlist(random_pct)[!is.nan(unlist(random_pct))],
                                 cyclic=unlist(cyclic_pct)[!is.nan(unlist(random_pct))]))
             })
    
    
    
    
    random_df_merged = pct_df[[1]][,c(1,2)]
    cyclic_df_merged = pct_df[[1]][,c(1,3)]
    
    cor_df_merged = cells_pvalue_df[[1]][,1:2]
    
    for(i in 2:len(cells_pvalue_df)) {
      # random_df_merged <- merge(random_df_merged, 
      #                           pct_df[[i]][,c(1,2)],
      #                           by="bins",
      #                           all=T)
      # 
      # cyclic_df_merged <- merge(cyclic_df_merged, 
      #                           pct_df[[i]][,c(1,3)],
      #                           by="bins",
      #                           all=T)  
      
      cor_df_merged <- merge(cor_df_merged,
                             cells_pvalue_df[[i]][,1:2],
                             by="bins",
                             all=T)
    }
    
    rownames(cor_df_merged) = cor_df_merged[,1]
    cor_df_merged = cor_df_merged[,-1]
    rownames(random_df_merged) <- random_df_merged[,1]
    random_df_merged <- random_df_merged[,-1]
    rownames(random_df_merged) <- 
      unlist(lapply(strsplit(rownames(random_df_merged), ","), function(sp) {strsplit(sp[2], "]")[[1]][1]}))
    
    rownames(cyclic_df_merged) <- cyclic_df_merged[,1]
    cyclic_df_merged <- cyclic_df_merged[,-1]
    rownames(cyclic_df_merged) <- 
      unlist(lapply(strsplit(rownames(cyclic_df_merged), ","), function(sp) {strsplit(sp[2], "]")[[1]][1]}))
    
    return(list(random_df=random_df_merged, 
                cyclic_df=cyclic_df_merged, 
                all=pct_df,
                subfields=subfields))
  }

generate_spike_train_figures <- 
  function(average_width,
           sd_width,
           noise_level,
           double_peak_percent,
           true_firing_rate,
           true_time_bins_per_cells,
           n_neur,
           stim_trace,
           place_cell_percentage,
           return_TC=F) {
    
    tuning_curve <-
      generate_tuning_curves_cost(n = n_neur,
                                  percentage = place_cell_percentage,
                                  average_width = average_width, 
                                  sd_width = sd_width,
                                  fixed_fr=true_firing_rate,
                                  noise=noise_level,
                                  double_peak_pct = double_peak_percent,
                                  plot=F)
    
    
    pois_factor <- currate_spike_train_cost(tuning_curve, true_time_bins_per_cells, stim_trace)
    gst <- generate_spike_trains_cost(tuning_curves = tuning_curve, stim_trace = stim_trace,
                                      factor=pois_factor, fs=1)
      
    if (return_TC) {
      return(list(gst, tuning_curve))
    }
    
    return(gst)
}


smooth_spline_vector <- function(vec, smooth_lambda) {
  barp <- barplot(vec, plot=F)
  smt <- smooth.spline(barp, vec, all.knots = T, lambda=smooth_lambda)
  smt$y[smt$y < 0] <- 0
  return(smt$y);
}

plot_gen_place_cell_histograms <- 
  function(average_width,
           sd_width,
           noise_level,
           double_peak_percent,
           true_spike_train,
           stim_trace,
           place_cell_percentage,
           pval_func=wilcox.test,
           peaks_significance_threshold=0.3,
           smooth_lambda=-1,
           hist_density=T,
           breaks_SI=31,
           breaks_tb=31,
           breaks_sb=15,
           breaks_peaks=8)  {
    
    firing_rate <-rowMeans(true_spike_train) / dt
    processed_real <- preprocess_spike_train(true_spike_train, stim_trace=stim_trace, verbose=F)
    true_cells_spike_train <- processed_real$working_cells_spike_train
    true_firing_rate <- processed_real$working_firing_rate
    true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
    
    tmp <- compute_tuning(true_cells_spike_train, stim_trace)
    stim_prob <- tmp[[1]]
    true_tuning_curve <- tmp[[2]]
    rm(tmp)
    
    
    if (smooth_lambda == - 1){
      tuning_for_peaks <- true_tuning_curve
    } else {
      tuning_for_peaks <- 
        t(apply(true_tuning_curve, 1, function(trc) {
          barp <- barplot(trc, plot=F)
          smt <- smooth.spline(barp, trc, all.knots = T, lambda=smooth_lambda)
          smt$y[smt$y < 0] <- 0
          return(smt$y);
        }))
    }
    
    true_peaks <- unlist(apply(tuning_for_peaks, 1, 
                               function(n){return(len(get_peaks(n, threshold_from_max = peaks_significance_threshold)))}))
    true_active_spatial_bins_per_cell <- apply(true_tuning_curve, 1, function(n){sum(n>0)})
    #print(true_active_spatial_bins_per_cell)
    true_SI <- compute_SI(stim_prob, true_tuning_curve, true_firing_rate)
    

    tuning_curve <-
      generate_tuning_curves_cost(n = nrow(true_spike_train),
                                  percentage = place_cell_percentage,
                                  average_width = average_width, 
                                  sd_width = sd_width,
                                  fixed_fr=firing_rate,
                                  noise=noise_level,
                                  double_peak_pct = double_peak_percent,
                                  plot=F)
    
    
      pois_factor <- currate_spike_train_cost(tuning_curve, 
                                              true_time_bins_per_cells,
                                              stim_trace, verbose=F)
      
      
      
    generated_spike_train <- 
      generate_spike_trains_cost(tuning_curves = tuning_curve,
                                 stim_trace = stim_trace,
                                 factor=pois_factor,
                                 fs=1)
    
    processed_generated <- preprocess_spike_train(generated_spike_train, stim_trace, verbose=F)
    generated_cells_active_st <- processed_generated$working_cells_spike_train
    generated_firing_rate <- processed_generated$working_firing_rate
    generated_time_bins <- processed_generated$working_time_bins_per_cells
    
    tmp <- compute_tuning(generated_cells_active_st, stim_trace)
    gen_stim_prob <- tmp[[1]]
    gen_tuning_curve <- tmp[[2]]
    rm(tmp)
    

    
    if (smooth_lambda == - 1){
      tuning_for_peaks <- gen_tuning_curve
    } else {
      tuning_for_peaks <- 
        t(apply(gen_tuning_curve, 1, function(trc) {
          barp <- barplot(trc, plot=F)
          smt <- smooth.spline(barp, trc, all.knots = T, lambda=smooth_lambda)
          smt$y[smt$y < 0] <- 0
          return(smt$y);
        }))
    }
    
    generated_peaks <- unlist(apply(tuning_for_peaks, 1, 
                                    function(n){return(len(get_peaks(n, threshold_from_max = peaks_significance_threshold)))}))
    generated_active_spatial_bins_per_cell <- apply(gen_tuning_curve, 1, function(n){sum(n>0)})
    generated_SI <- compute_SI(gen_stim_prob, gen_tuning_curve, generated_firing_rate)
    ks.test(generated_SI[[1]], true_SI[[1]])
    
    max_SI <- max(c(generated_SI[[1]], true_SI[[1]]))
    SI_breaks =seq(0, max_SI, length.out=breaks_SI)
    hgen <- hist(generated_SI[[1]], breaks=SI_breaks, plot=F)$counts
    htrue <- hist(true_SI[[1]], breaks=SI_breaks, plot=F)$counts

    if (hist_density) {
      SI_df <- data.frame(counts_gen=smooth_spline_vector(hgen, 5e-5),
                          counts_true=smooth_spline_vector(htrue, 5e-5),
                          SI=SI_breaks[2:breaks_SI])
      
      
        gSI <- 
          ggplot(SI_df, aes(y=counts_true, x=SI)) + 
          geom_bar(stat="identity", 
                   color="black", 
                   position="identity", 
                   fill="cornflowerblue",
                   alpha=0.8) + 
          theme_light() +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                legend.position="NA",
                panel.border = element_blank(),
                panel.background = element_blank()) +
          xlab("SI (bit/spikes)") +
          ylab("Frequency") +
          ggtitle(sprintf("Statistic=%f; Pval=%f",
                          pval_func(generated_SI[[1]], true_SI[[1]])$statistic,
                          pval_func(generated_SI[[1]], true_SI[[1]])$p.value)) +
          geom_line(aes(y=counts_gen), size=1, linetype="longdash")
        
    } else {
      SI_df <- 
        data.frame(Counts=c(htrue,hgen),
                   SI=c(SI_breaks[2:breaks_SI], SI_breaks[2:breaks_SI]),
                   Type=c(rep("Generated", times=breaks_SI - 1), rep("True", times=breaks_SI - 1)))
      
      
      gSI <- 
        ggplot(SI_df, aes(y=Counts, x=SI)) + 
        geom_bar(aes(fill=Type), 
                 stat="identity", 
                 color="black", 
                 position="identity", 
                 alpha=0.8) + 
        theme_light() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.position="NA",
              panel.border = element_blank(),
              panel.background = element_blank()) +
        xlab("Spatial information (bit/spikes)") +
        ylab("Frequency") +
        scale_fill_manual(values=c("cornflowerblue", "cornflowerblue")) + 
        ggtitle(sprintf("Statistic=%f; Pval=%f",
                              pval_func(generated_SI[[1]], true_SI[[1]])$statistic,
                              pval_func(generated_SI[[1]], true_SI[[1]])$p.value))
   
      
    }
    # Compare time bins distribution distances
    max_time_bins <- max(generated_time_bins, true_time_bins_per_cells)
    time_bins_breaks = seq(0, max_time_bins, length.out=breaks_tb)
    hgen <- hist(generated_time_bins, breaks=time_bins_breaks, plot=F)$counts
    htrue <- hist(true_time_bins_per_cells, breaks=time_bins_breaks, plot=F)$counts

    
    if (hist_density) {
      time_bins_df <- data.frame(counts_gen=smooth_spline_vector(hgen, 5e-5),
                          counts_true=smooth_spline_vector(htrue, 5e-5),
                          Tbins=time_bins_breaks[2:breaks_tb])
      
      
      gtimebins <- 
        ggplot(time_bins_df, aes(y=counts_true, x=Tbins)) + 
        geom_bar(stat="identity", 
                 color="black", 
                 position="identity", 
                 fill="cornflowerblue",
                 alpha=0.8) + 
        theme_light() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.position="NA",
              panel.border = element_blank(),
              panel.background = element_blank()) +
        xlab("Time bins (#)") +
        ylab("Frequency") +
        ggtitle(sprintf("Statistic=%f; Pval=%f",
                        pval_func(generated_time_bins, true_time_bins_per_cells)$statistic,
                        pval_func(generated_time_bins, true_time_bins_per_cells)$p.value)) +
        geom_line(aes(y=counts_gen), size=1, linetype="longdash")
      
    } else {
      time_bins_df <- 
        data.frame(Counts=c(hgen, htrue),
                   Tbins=c(time_bins_breaks[2:breaks_tb], time_bins_breaks[2:breaks_tb]),
                   Type=c(rep("Generated", times=breaks_tb - 1), rep("True", times=breaks_tb - 1)))
      gtimebins <- 
        ggplot(time_bins_df, aes(y=Counts, x=Tbins)) + 
        geom_bar(aes(fill=Type), 
                 stat="identity", 
                 color="black", 
                 position="identity", 
                 alpha=0.8) + 
        theme_light() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.position="NA",
              panel.border = element_blank(),
              panel.background = element_blank()) +
        xlab("Active time bins (#)") +
        ylab("Frequency") +
        scale_fill_manual(values=c("cornflowerblue", "cornflowerblue")) + 
        ggtitle(sprintf("Statistic=%f; Pval=%f",
                        pval_func(generated_time_bins, true_time_bins_per_cells)$statistic,
                        pval_func(generated_time_bins, true_time_bins_per_cells)$p.value))
    }
    
    # Compare active spatial bins distribution distances
    max_spatial_bins <- max(generated_active_spatial_bins_per_cell, true_active_spatial_bins_per_cell)
    spat_bins_breaks =seq(0, max_spatial_bins, length.out=breaks_sb)
    hgen <- hist(generated_active_spatial_bins_per_cell, breaks=spat_bins_breaks, plot=F)$counts
    htrue <- hist(true_active_spatial_bins_per_cell, breaks=spat_bins_breaks, plot=F)$counts
    
    if (hist_density) {
      spat_bins_df <- data.frame(counts_gen=smooth_spline_vector(hgen, 5e-5),
                                 counts_true=smooth_spline_vector(htrue, 5e-5),
                                 Sbins=spat_bins_breaks[2:breaks_sb])
      
      
      gspatbins <- 
        ggplot(spat_bins_df, aes(y=counts_true, x=Sbins)) + 
        geom_bar(stat="identity", 
                 color="black", 
                 position="identity", 
                 fill="cornflowerblue",
                 alpha=0.8) + 
        theme_light() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.position="NA",
              panel.border = element_blank(),
              panel.background = element_blank()) +
        xlab("Spatial bins (#)") +
        ylab("Frequency") + 
        ggtitle(sprintf("Statistic=%f; Pval=%f",
                        pval_func(generated_active_spatial_bins_per_cell, true_active_spatial_bins_per_cell)$statistic,
                        pval_func(generated_active_spatial_bins_per_cell, true_active_spatial_bins_per_cell)$p.value)) + 
        geom_line(aes(y=counts_gen), size=1, linetype="longdash")
      
    } else {
    
    spat_bins_df <- 
      data.frame(Counts=c(hgen, htrue),
                 Sbins=c(spat_bins_breaks[2:breaks_sb], spat_bins_breaks[2:breaks_sb]),
                 Type=c(rep("Generated", times=breaks_sb - 1), rep("True", times=breaks_sb - 1)))
    
    gspatbins <- 
      ggplot(spat_bins_df, aes(y=Counts, x=Sbins)) + 
      geom_bar(aes(fill=Type), 
               stat="identity", 
               color="black", 
               position="identity", 
               alpha=0.8) + 
      theme_light() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position="NA",
            panel.border = element_blank(),
            panel.background = element_blank()) +
      xlab("Active spatial bins (#)") +
      ylab("Frequency") + 
      scale_fill_manual(values=c("cornflowerblue", "cornflowerblue")) + 
      ggtitle(sprintf("Statistic=%f; Pval=%f",
                      pval_func(generated_active_spatial_bins_per_cell, true_active_spatial_bins_per_cell)$statistic,
                      pval_func(generated_active_spatial_bins_per_cell, true_active_spatial_bins_per_cell)$p.value))
    
    }
    
    # Compare peaks bins distribution distances
    max_peaks <- max(generated_peaks, true_peaks)
    peaks_breaks =seq(0, max_peaks, length.out=breaks_peaks)
    hgen <- hist(generated_peaks, breaks=peaks_breaks, plot=F)$counts
    htrue <- hist(true_peaks,      breaks=peaks_breaks, plot=F)$counts
    
    if (hist_density) {
      peaks_df <- data.frame(counts_gen=smooth_spline_vector(hgen, 5e-5),
                             counts_true=smooth_spline_vector(htrue, 5e-5),
                             Peaks=peaks_breaks[2:breaks_peaks])
      
      
      gpeaks <- 
        ggplot(peaks_df, aes(y=counts_true, x=Peaks)) + 
        geom_bar(stat="identity", 
                 color="black", 
                 position="identity", 
                 fill="cornflowerblue",
                 alpha=0.8) + 
        theme_light() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.position="NA",
              panel.border = element_blank(),
              panel.background = element_blank()) +
        xlab("Peaks (#)") +
        ylab("Frequency") + 
        ggtitle(sprintf("Statistic=%f; Pval=%f",
                        ks.test(generated_peaks, true_peaks)$statistic,
                        pval_func(generated_peaks, true_peaks)$p.value)) +  
        geom_line(aes(y=counts_gen), size=1, linetype="longdash")
      
    } else {
    
    peaks_df <- data.frame(Counts=c(hgen, htrue),
                           Peaks=c(peaks_breaks[2:breaks_peaks], peaks_breaks[2:breaks_peaks]),
                           Type=c(rep("Generated", times=breaks_peaks - 1), rep("True", times=breaks_peaks - 1)))
    
    gpeaks <- 
      ggplot(peaks_df, aes(y=Counts, x=Peaks)) + 
      geom_bar(aes(fill=Type), 
               stat="identity", 
               color="black", 
               position="identity", 
               alpha=0.8) + 
      theme_light() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position="NA",
            panel.border = element_blank(),
            panel.background = element_blank()) +
      xlab("Peaks (#)") +
      ylab("Frequency") + 
      scale_fill_manual(values=c("cornflowerblue", "cornflowerblue")) + 
      ggtitle(sprintf("Statistic=%f; Pval=%f",
                      pval_func(generated_peaks, true_peaks)$statistic,
                      pval_func(generated_peaks, true_peaks)$p.value))
    
    }
    
    
    
    
    return(list(gSI, gtimebins, gspatbins, gpeaks))
    
  }


plot_gen_place_cell_qqplots <- 
  function(average_width,
           sd_width,
           noise_level,
           double_peak_percent,
           true_spike_train,
           stim_trace,
           place_cell_percentage,
           pval_func=wilcox.test,
           peaks_significance_threshold=0.3,
           smooth_lambda=-1,
           hist_density=T,
           breaks_SI=31,
           breaks_tb=31,
           breaks_sb=15,
           breaks_peaks=8)  {
    
    firing_rate <-rowMeans(true_spike_train) / dt
    processed_real <- preprocess_spike_train(true_spike_train, stim_trace=stim_trace, verbose=F)
    true_cells_spike_train <- processed_real$working_cells_spike_train
    true_firing_rate <- processed_real$working_firing_rate
    true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
    
    tmp <- compute_tuning(true_cells_spike_train, stim_trace)
    stim_prob <- tmp[[1]]
    true_tuning_curve <- tmp[[2]]
    rm(tmp)
    
    
    if (smooth_lambda == - 1){
      tuning_for_peaks <- true_tuning_curve
    } else {
      tuning_for_peaks <- 
        t(apply(true_tuning_curve, 1, function(trc) {
          barp <- barplot(trc, plot=F)
          smt <- smooth.spline(barp, trc, all.knots = T, lambda=smooth_lambda)
          smt$y[smt$y < 0] <- 0
          return(smt$y);
        }))
    }
    
    true_peaks <- unlist(apply(tuning_for_peaks, 1, 
                               function(n){return(len(get_peaks(n, threshold_from_max = peaks_significance_threshold)))}))
    true_active_spatial_bins_per_cell <- apply(true_tuning_curve, 1, function(n){sum(n>0)})
    #print(true_active_spatial_bins_per_cell)
    true_SI <- compute_SI(stim_prob, true_tuning_curve, true_firing_rate)
    
    
    tuning_curve <-
      generate_tuning_curves_cost(n = nrow(true_spike_train),
                                  percentage = place_cell_percentage,
                                  average_width = average_width, 
                                  sd_width = sd_width,
                                  fixed_fr=firing_rate,
                                  noise=noise_level,
                                  double_peak_pct = double_peak_percent,
                                  plot=F)
    
    
    pois_factor <- currate_spike_train_cost(tuning_curve, 
                                            true_time_bins_per_cells,
                                            stim_trace, verbose=F)
    
    
    
    generated_spike_train <- 
      generate_spike_trains_cost(tuning_curves = tuning_curve,
                                 stim_trace = stim_trace,
                                 factor=pois_factor,
                                 fs=1)
    
    processed_generated <- preprocess_spike_train(generated_spike_train, stim_trace, verbose=F)
    generated_cells_active_st <- processed_generated$working_cells_spike_train
    generated_firing_rate <- processed_generated$working_firing_rate
    generated_time_bins <- processed_generated$working_time_bins_per_cells
    
    tmp <- compute_tuning(generated_cells_active_st, stim_trace)
    gen_stim_prob <- tmp[[1]]
    gen_tuning_curve <- tmp[[2]]
    rm(tmp)
    
    
    
    if (smooth_lambda == - 1){
      tuning_for_peaks <- gen_tuning_curve
    } else {
      tuning_for_peaks <- 
        t(apply(gen_tuning_curve, 1, function(trc) {
          barp <- barplot(trc, plot=F)
          smt <- smooth.spline(barp, trc, all.knots = T, lambda=smooth_lambda)
          smt$y[smt$y < 0] <- 0
          return(smt$y);
        }))
    }
    
    generated_peaks <- unlist(apply(tuning_for_peaks, 1, 
                                    function(n){return(len(get_peaks(n, threshold_from_max = peaks_significance_threshold)))}))
    generated_active_spatial_bins_per_cell <- apply(gen_tuning_curve, 1, function(n){sum(n>0)})
    generated_SI <- compute_SI(gen_stim_prob, gen_tuning_curve, generated_firing_rate)
  
    get_quantile_df <- function(a,b, quantiles=seq(0.01, 1, length.out=100)) {sapply(quantiles, function(q) {c(quantile(a, q), quantile(b,q))})}
    
    
    qq_SI <- t(get_quantile_df(true_SI[[1]], generated_SI[[1]]));
    qq_timebins  <- t(get_quantile_df(true_time_bins_per_cells, generated_time_bins));
    qq_sbins  <- t(get_quantile_df(true_active_spatial_bins_per_cell, generated_active_spatial_bins_per_cell));
    qq_peaks <- t(get_quantile_df(true_peaks, generated_peaks))
    
    print(dim(qq_SI))
    
    colnames(qq_SI) <- c("True", "Generated")
    colnames(qq_sbins) <- c("True", "Generated")
    colnames(qq_timebins) <- c("True", "Generated")
    colnames(qq_peaks) <- c("True", "Generated")
    
    print(qq_SI)
    gSI <- 
      ggplot(as.data.frame(qq_SI), aes(x=True, y=Generated)) + 
      geom_line() + 
      geom_line(data=data.frame(x=c(min(qq_SI), max(qq_SI)),
                                y=c(min(qq_SI), max(qq_SI))),
                aes(x=x,y=y))

    gspatbins <- 
      ggplot(as.data.frame(qq_sbins), aes(x=True, y=Generated)) + 
      geom_line() +
      geom_line(data=data.frame(x=c(min(qq_sbins), max(qq_sbins)),
                                y=c(min(qq_sbins), max(qq_sbins))),
                aes(x=x,y=y))
    
    gtimebins <- 
      ggplot(as.data.frame(qq_timebins), aes(x=True, y=Generated)) + 
      geom_line() + 
      geom_line(data=data.frame(x=c(min(qq_timebins), max(qq_timebins)),
                                y=c(min(qq_timebins), max(qq_timebins))),
                aes(x=x,y=y))
    
    gpeaks <- 
      ggplot(as.data.frame(qq_peaks), aes(x=True, y=Generated)) + 
      geom_line() + 
      geom_line(data=data.frame(x=c(min(qq_peaks), max(qq_peaks)),
                                y=c(min(qq_peaks), max(qq_peaks))),
                aes(x=x,y=y))
    
    return(list(gSI, gtimebins, gspatbins, gpeaks))
    
  }


figure_1A_1B <- function(stim_trace, 
                         spike_train,
                         simulated_tuning_curve=list(),
                         n_to_use = 10,
                         session_length=20, 
                         dt=0.05,
                         out_path,
                         used_pois_factor=1.5,
                         ext="") {
  
  figs_path <- sprintf("%s\\%s%s", out_path, "figure_1AB_pdf", ext)
  dir.create(figs_path)
  pval_by_time = pct_by_sample_size_figures(stim_trace = stim_trace,
                                            spike_train=spike_train)
  
  decreasing_rand <- apply(t(apply(pval_by_time$rand, 1, diff)), 1, function(r) {sum(r < 0)})
  
  indices <- order(decreasing_rand, decreasing=T)[1:50]
  run_df <- data.frame(Time=((1:len(stim_trace) / len(stim_trace)) * session_length), 
                       Position=stim_trace * 4) 
  run_df$Position[which(run_df$Position == 4)] <- 0
  actual_times <- (pval_by_time$duration / (ncol(spike_train) * dt)) * session_length
  
  
  for (idx in 1:nrow(spike_train)) {
  

  
  firing_ind <- which(spike_train[idx,] > 0)
  
  cell_df <- data.frame(Positions=stim_trace[firing_ind] * 4,
                        Times=run_df$Time[firing_ind])
  
  pval_df_all <- data.frame(Pvalue=c(pval_by_time$rand[idx,], pval_by_time$cyclic[idx,]),
                        Time=c(actual_times,actual_times),
                        `Shuffle type`=c(rep("Random", times=length(actual_times)),
                                  rep("Cyclic", times=length(actual_times))))
  
  pval_df <- data.frame(Pvalue=c(pval_by_time$cyclic[idx,]),
                        Time=c(actual_times),
                        `Shuffle type`=c(rep("Cyclic", times=length(actual_times))))
  
  tuning_curve <- (unlist(lapply(sort(unique(stim_trace)), function (sb) {mean(spike_train[idx,which(stim_trace == sb)]) / dt}))) 
  barp <- barplot(tuning_curve)
  smoothed_tuning <- smooth.spline(barp, tuning_curve, all.knots = T, lambda=1e-4)
  smoothed_tuning$y[smoothed_tuning$y < 0] <- 0
  pos <- sort(unique(stim_trace)) * 4; pos[pos == 4] <- 0
  smoothed_df <- data.frame(FR=smoothed_tuning$y,
                            Positions=pos)
  
  if (len(simulated_tuning_curve) != 0) {
    sim_barp <- barplot(simulated_tuning_curve[idx,])
    sim_smoothed_tuning <- smooth.spline(sim_barp, simulated_tuning_curve[idx,], all.knots = T, lambda=1e-4)
    sim_smoothed_tuning$y[sim_smoothed_tuning$y < 0] <- 0
    sim_smoothed_df <- data.frame(FR=sim_smoothed_tuning$y / (dt * pois_factor),
                                  Positions=pos)
  }
  
  polygon_bin <- which(pval_df$Pvalue[2:30] < 0.05)[1]
  if(len(polygon_bin) > 0) {
  polygon_df <- polygon_df <- data.frame(x=c(pval_df$Time[polygon_bin],# - (pval_df$Time[polygon_bin] - pval_df$Time[polygon_bin-1]) * 0.5,
                                             20, 20,
                                             pval_df$Time[polygon_bin]),# -(pval_df$Time[polygon_bin] - pval_df$Time[polygon_bin-1]) * 0.5),
                                         y=c(0,0,96,96))
  signif_text_pos = ifelse(polygon_bin - 5 > 0, polygon_bin -5, polygon_bin +5)
  } else {
    polygon_df(x=c(20,20,20,20), y=c(0,0,96,96))
    signif_text_pos = 4
  }
  
  grundf <-
          ggplot(polygon_df, aes(x=x,y=y)) + 
          geom_polygon(fill=adjustcolor("deepskyblue1", alpha=0.2)) +
          #ggplot(run_df) +
          geom_line(data=run_df, aes(x=Time, y=Position),
                    col="gray30") + 
          theme_light() +  
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                legend.position="NA",
                #panel.border = element_blank(),
                panel.background = element_blank()) +
          xlim(20,0) +#plot.margin = margin(t=0,r=0,b=0,l=0,unit="cm")) +
          ylim(0,96) +
          ylab("Position (cm)") +
          xlab("Time (minutes)") +
          geom_vline(xintercept = polygon_df$x[1],
                     linetype="dashed",
                      color="black") +
          geom_point(data=cell_df, aes(x=Times,y=Positions),
                     fill="red",
                     color="red",
                     size=2) + 
          coord_flip()
          

  
  
  gtuning_df <- ggplot(smoothed_df) + 
                geom_line(aes(x=Positions, y=FR), col="red", size=2) +
                theme_light() +
                ylab("Firing rate (spikes/sec)") + 
                xlab("Position (cm)") +
                xlim(0,96) +
                ylim(0,ifelse(max(smoothed_df$FR) > 1.5, max(smoothed_df$FR), 1.5)) +
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "black"),
                    legend.position="NA",
                    panel.border = element_blank(),
                    panel.background = element_blank()) 
  
  if (len(simulated_tuning_curve) != 0) { 
    
    scale_factor = 0.5
    combined_tuning_df <- rbind(sim_smoothed_df, smoothed_df)
    combined_tuning_df$Type <- rep(c("Defined tuning curve", "Measured"), each=nrow(smoothed_df))
    
    gtuning_df <- ggplot(combined_tuning_df) + 
      geom_line(aes(x=Positions, y=FR, group=Type, color=Type), size=2, alpha=0.8) +
      theme_light() +
      ylab("Firing rate (spikes/sec)") + 
      xlab("Position (cm)") +
      scale_color_manual(breaks=c("Defined tuning curve", "Measured"), values=c("royalblue4", "red")) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position="NA",
            panel.border = element_blank(),
            panel.background = element_blank()) +
      ylim(c(0,ifelse(max(smoothed_df$FR) > 1.5, max(smoothed_df$FR), 1.5))) +
      xlim(0,96)

            
    
  }
  #plot(grundf)
  #dev.off()
  gpvaldf <- ggplot(pval_df) + 
             geom_line(aes(x=Time, y=Pvalue, color=`Shuffle.type`)) + 
             geom_point(aes(x=Time, y=Pvalue, color=`Shuffle.type`)) + 
             theme_light() +  
             theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(colour = "black"),
                  legend.position="NA",
                  panel.border = element_blank(),
                  panel.background = element_blank()) + 
             ylab(TeX(r'($P_{value}$)')) +
             ylim(0,1) + 
             xlim(0, 20) +
             geom_hline(yintercept = 0.05, linetype="dashed", col="gray60") +
             # geom_text(y=0.08, x=17.3, label="Significance threshold",
             #           col="gray60",
             #           size=2.5) +
             scale_color_manual(name = "Shuffle type", 
                                  breaks = c("Cyclic", "Random"),
                                  values = c("mediumorchid2", "lightsalmon2"))
     
  #   
  #             
  # tiff(filename=sprintf("%s\\%d_pval_df.tiff", 
  #                       out_path, 
  #                       idx),
  #      units="in",
  #      width=4,
  #      height=5,
  #      res=500)
  res <- align_plots(grundf, gtuning_df, gpvaldf, align="v")
  g <- grid.arrange(res[[1]], res[[2]], res[[3]], heights=c(1.8,1,1))
  pdf(file=sprintf("%s\\%d_1ab.pdf", 
                       figs_path, 
                        idx),
       width=4,
       height=3.8 * 4)
  plot(g)
  dev.off()       
  
  png(file=sprintf("%s\\%d_1ab.png", 
                   figs_path, 
                   idx),
      units="in",
      res=500,
      width=4,
      height=10)
  plot(g)
  dev.off()      
  # if (len(simulated_tuning_curve) != 0) {
  #   g <- grid.arrange(grundf, gsimtuning_df, gtuning_df, gpvaldf, nrow=4, heights=c(1.8,1,1, 1))
  #   pdf(file=sprintf("%s\\%d_1ab_sim_tuning.pdf", 
  #                        figs_path, 
  #                        idx),
  #       width=4,
  #       height=12.5)
  #   plot(g)
  #   dev.off()       
  # 
  #   png(file=sprintf("%s\\%d_1ab_sim_tuning.png", 
  #                    figs_path, 
  #                    idx),
  #       units="in",
  #       res=500,
  #       width=4,
  #       height=12.5)
  #   plot(g)
  #   dev.off()           
  # }
  # 
  }
 }

plot_place_cell_pct_by_activity_plot <- function(paths, activity_bin_size =15, cyclic_only=F, fit_sim="") {
  pct_by_activity <- create_pct_by_activity_df(paths, activity_bin_size = activity_bin_size, fit_sim = fit_sim)
  
  random_df <- pct_by_activity$random_df
  cyclic_df <- pct_by_activity$cyclic_df
  random_bins <- as.numeric(rownames(random_df))
  cyclic_bins <- as.numeric(rownames(cyclic_df))
  bins_order_cyc <- order(cyclic_bins)
  bins_order_rand <- order(random_bins)
  
  random_df_means <- rowMeans(random_df, na.rm=T)[bins_order_rand]
  cyclic_df_means <- rowMeans(cyclic_df, na.rm=T)[bins_order_cyc]
  random_df_sem <- apply(random_df, 1, sd, na.rm=T)[bins_order_rand]
  cyclic_df_sem <- apply(cyclic_df, 1, sd, na.rm=T)[bins_order_cyc]
  random_bins <- sort(random_bins) # same as order
  cyclic_bins <- sort(cyclic_bins) # same as order
  
  
  if (cyclic_only) {
    plot_df <- data.frame(fraction=c(cyclic_df_means),
                          sem=c(cyclic_df_sem),
                          activity=(c(cyclic_bins)))
    
    plot_df$Shuffle <- c(rep("Cyclic shuffle", times=len(cyclic_bins)))    
  } else {
  plot_df <- data.frame(fraction=c(random_df_means, cyclic_df_means),
                        sem=c(random_df_sem, cyclic_df_sem),
                        activity=c(random_bins, cyclic_bins))
  
  plot_df$Shuffle <- c(rep("Random shuffle", times=len(random_bins)),
                       rep("Cyclic shuffle", times=len(cyclic_bins)))
  }
  
  plot_df$sem[is.na(plot_df$sem)] <- 0
  #plot_df$sem <- plot_df$sem + 0.1
  g <- ggplot(plot_df, aes(x=activity, y=fraction, color=Shuffle)) + 
    geom_line(size=1.2, aes(color=Shuffle, group=Shuffle)) +
    geom_linerange(aes(ymin=fraction-sem, ymax=fraction+sem, group=Shuffle, color=Shuffle),size=1) +    
    scale_color_manual(values = c("darkmagenta", "lightsalmon2")) +
    theme_light() +
    ylim(c(0, 1.1))+
    xlim((c(min(cyclic_bins),max(cyclic_bins)))) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position="NA",
          panel.border = element_blank(),
          panel.background = element_blank()) + 
    labs(x="Active time bins (frames)", y = "Fraction of place cells (%)")
  
  gn <- ggplot(plot_df, aes(x=activity, y=fraction, color=Shuffle)) + 
    geom_line(size=1.2, aes(color=Shuffle, group=Shuffle)) +
    geom_ribbon(aes(ymin=fraction-sem, ymax=fraction+sem, group=Shuffle, fill=Shuffle), color="NA", alpha=0.1, size=1) +
    scale_color_manual(values = c("darkmagenta", "lightsalmon2")) +
    scale_fill_manual(values = c("darkmagenta", "lightsalmon2")) +
    theme_light() +
    ylim(c(0, 1.1))+
    xlim((c(min(cyclic_bins),max(cyclic_bins)))) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position="NA",
          panel.border = element_blank(),
          panel.background = element_blank()) + 
    labs(x="Active time bins (frames)", y = "Fraction of place cells (%)") +
    geom_hline(yintercept=1, linetype="dashed")
  
  gn2 <- ggplot(plot_df, aes(x=activity, y=fraction, color=Shuffle)) + 
    geom_line(size=1.2, aes(color=Shuffle, group=Shuffle)) +
    geom_errorbar(aes(ymin=fraction-sem, ymax=fraction+sem, group=Shuffle, fill=Shuffle), size=1, width=0.5) +
    scale_color_manual(values = c("darkmagenta", "lightsalmon2")) +
    scale_fill_manual(values = c("darkmagenta", "lightsalmon2")) +
    theme_light() +
    ylim(c(0, 1.1))+
    xlim(c(0,135)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position="NA",
          panel.border = element_blank(),
          panel.background = element_blank()) + 
    labs(x="Active time bins (frames)", y = "Fraction of place cells (%)") 
  
  
  sub_ind <- !is.na(cyclic_df[,1][bins_order_cyc])
  
  tmp_df <- data.frame(fraction=cyclic_df[,1][bins_order_cyc][sub_ind],  activity=cyclic_bins[sub_ind])
  
  g2 <- g + 
    geom_line(data=tmp_df,aes(x=activity, y=fraction), col="darkmagenta", alpha=0.2, size=0.6)
  
  for (i in 2:ncol(cyclic_df)) {
    sub_ind <- !is.na(cyclic_df[,i][bins_order_cyc])
    
    tmp_df <- data.frame(fraction=cyclic_df[,i][bins_order_cyc][sub_ind],  activity=cyclic_bins[sub_ind])
    
    g2 <- g2 + 
      geom_line(data=tmp_df,aes(x=activity, y=fraction), col="darkmagenta", alpha=0.2, size=0.6)
  }
  
  return(list(gn, g2, plot_df, cyclic_df, gn2))
}


figure_1C <- 
  function(paths,
          out_path,
          ext="") { 
  
  figs_path <- sprintf("%s\\%s", out_path, "figure_1C")
  dir.create(figs_path)
  print(figs_path)
  true_data_paths <- sprintf("%s\\%s", paths, "equalized")
  
  
  
  tplots <- plot_place_cell_pct_by_activity_plot(true_data_paths, cyclic_only=T)
  fplots_1c <- plot_place_cell_pct_by_activity_plot(true_data_paths, cyclic_only=T, fit_sim = "fit_simulation_new_v_fit\\")
  fplots_1 <- plot_place_cell_pct_by_activity_plot(true_data_paths, cyclic_only=T, fit_sim = "fit_simulation\\")
  #simulated_data_plots <- plot_place_cell_pct_by_activity_plot(simulated_data_paths)
  
  pdf(file=sprintf("%s\\%s1c.pdf", figs_path,ext),
       width=5,
       height=5)
  
  # g <- grid.arrange(true_data_plots, 
  #                   simulated_data_plots, nrow=2)
  plot(true_data_plots[[2]])
  dev.off()         
}
figure_1D_1E <- function(paths, 
                         out_path,
                         ext="") {
  
  figs_path <- sprintf("%s\\%s", out_path, "figure_1DE")
  dir.create(figs_path)
  print(figs_path)
  true_data_paths <- sprintf("%s\\%s", paths, "equalized")
  simulated_data_paths <- sprintf("%s\\%s", paths, "simulated\\equalized")
  fit_sim_paths <- sprintf("%s\\%s", paths, "equalized\\fit_simulation")
  
  tdp <- plot_place_cell_pct_by_sample_duration_figure(true_data_paths, bin_size = 2.5, ext = "_2", cyclic_only=T, fit_sim = "", verbose=T)
  fit_c <- plot_place_cell_pct_by_sample_duration_figure(true_data_paths, bin_size = 2.5, ext = "_2", cyclic_only=T, fit_sim = "fit_simulation_new_v_fit\\", verbose=T)
  fit_w <- plot_place_cell_pct_by_sample_duration_figure(true_data_paths, bin_size = 2.5, ext = "_2", cyclic_only=T, fit_sim = "fit_simulation_new_v\\", verbose=F)
  
  
  pdf(file=sprintf("%s\\true_%s1d.pdf",figs_path,ext),
       width=5,
       height=5)
  plot(true_data_plots[[3]])
  dev.off()       
  
  pdf(file=sprintf("%s\\sim_%s1d.pdf",figs_path,ext),
      width=5,
      height=5)
  plot(simulated_data_plots[[3]])
  dev.off()       
  
  # pdf(file=sprintf("%s\\%s1e.pdf", figs_path,ext),
  #      width=5,
  #      height=5)
  # 
  # plot(simulated_data_plots[[3]])
  # dev.off()         

}


figure_1F <- function(paths, 
                         out_path,
                         ext="") {
  
  figs_path <- sprintf("%s\\%s", out_path, "figure_1F")
  dir.create(figs_path)
  print(figs_path)
  true_data_paths <- sprintf("%s\\%s", paths, "equalized")
  
  true_p1f <- plot_pval_by_sample_duration_figure(true_data_paths, ext = "", fit_sim = "", verbose=F, type1 = F)
  fit_c1f <- plot_pval_by_sample_duration_figure(true_data_paths, ext = "", fit_sim = "fit_simulation_new_v_fit\\", verbose=F, type1=F)
  
  true_p1f_t2 <- plot_pval_by_sample_duration_figure(true_data_paths, ext = "", fit_sim = "", verbose=F, type1 = T)
  fit_c1f_t2 <- plot_pval_by_sample_duration_figure(true_data_paths, ext = "", fit_sim = "fit_simulation_new_v_fit\\", verbose=F, type1=T)
  
  pdf(file=sprintf("%s\\true_%s1d.pdf",figs_path,ext),
      width=5,
      height=5)
  plot(true_data_plots[[3]])
  dev.off()       
  
  pdf(file=sprintf("%s\\sim_%s1d.pdf",figs_path,ext),
      width=5,
      height=5)
  plot(simulated_data_plots[[3]])
  dev.off()       
  
    
  
}


figure_2D_2E_2F_2G <- 
         function(average_width,
                  sd_width,
                  noise_level,
                  double_peak_percent,
                  true_spike_train,
                  stim_trace,
                  out_path) 
{
          out_path <- sprintf("%s\\%s", out_path, "figure_2DEFG")
          dir.create(out_path)
  
           graphs <- plot_gen_place_cell_histograms(average_width = average_width,
                                                    sd_width = sd_width,
                                                    noise_level = noise_level,
                                                    double_peak_percent = double_peak_percent,
                                                    true_spike_train = true_spike_train,
                                                    stim_trace = stim_trace,
                                                    place_cell_percentage = 0.8)
           
           
            names(graphs) <- c("SI", "timebins", "spatbins", "peaks")
           
            for (gn in names(graphs)) {
              
              g <- graphs[[gn]]
              
              tiff(filename=sprintf("%s\\%s_compare_hist.tiff", 
                                    out_path, 
                                    gn),
                   units="in",
                   width=4,
                   height=4,
                   res=500)
              
              plot(g)
              dev.off()
            }

}                    
    
figure_2A_2B_2C <- 
  function(true_spike_train,
           stim_trace,
           out_path,
           ext="")
  {
    out_path <- sprintf("%s\\%s", out_path, "figure_2ABC")
    dir.create(out_path)
    n_cells = 200
    session_length = 20
    
    true_firing_rate <-rowMeans(true_spike_train) / 0.05
    processed_real <- preprocess_spike_train(true_spike_train, stim_trace=stim_trace, verbose=F)
    
    ### Generate a spike train
    tmp <- 
      generate_spike_train_figures(average_width=0.13,
                                   sd_width=0.007,
                                   noise_level=0.01,
                                   double_peak_percent=0.4,
                                   true_firing_rate,
                                   true_time_bins_per_cells = processed_real$working_time_bins_per_cells,
                                   nrow(true_spike_train),
                                   1, 
                                   return_TC = T)
    
    generated_spike_train <- tmp[[1]]
    tuning_curve <- tmp[[2]]
    rm(tmp)
    
    run_df <- data.frame(Time=((1:len(stim_trace) / len(stim_trace)) * session_length), 
                         Position=stim_trace * 4) 
    run_df$Position[which(run_df$Position == 4)] <- 0
    
    
    neur_ind <- sample(1:nrow(generated_spike_train), n_cells)
    measured_tuning <- compute_tuning(generated_spike_train, stim_trace)[[2]]
    
    for (idx in neur_ind) {
      print(idx)
      firing_ind <- which(generated_spike_train[idx,] > 0)
      cell_df <- data.frame(Positions=stim_trace[firing_ind] * 4, Times=run_df$Time[firing_ind])
      tc_df <- data.frame(created= c(0,tuning_curve[idx,]), 
                          measured=c(0,measured_tuning[idx,]),
                          position=c(0,sort(unique(stim_trace * 4))))
      
      grundf <-
        ggplot(run_df) +
        geom_line(data=run_df, aes(x=Time, y=Position),
                  col="gray50") + 
        theme_light() +  
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank()) + #plot.margin = margin(t=0,r=0,b=0,l=0,unit="cm")) +
        ylim(0,96) +
        xlim(0,20) +
        ylab("Position (cm)") +
        xlab("Time (minutes)") +
        geom_point(data=cell_df, aes(x=Times,y=Positions),
                   fill="red",
                   color="red",
                   size=2)
      
      gtuning_df <- ggplot(tc_df)
      
      gcreated <-  gtuning_df + 
                   geom_line(aes(x=position, y=created), col="red", size=2) +
                   theme_light() +
                   ylab("Firing rate (spikes/sec)") + 
                   xlab("Position (cm)") +
                   xlim(0,100)+
                   theme(panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank())
      
      gmeasured <- gtuning_df + 
        geom_bar(aes(x=position, y=measured),stat="identity" , col="black", size=1, fill="gray60") +
        theme_light() +
        ylab("Firing rate (spikes/sec)") + 
        xlab("Position (cm)") +
        xlim(0,100)+
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())
        
      tiff(filename=sprintf("%s\\%s%d_2abc.tiff", 
                            out_path, 
                            ext,
                            idx),
           units="in",
           width=10,
           height=8,
           res=500)
      
      g <- grid.arrange(gcreated, 
                        gmeasured,
                        grundf,
                        nrow=3)
      plot(g)
      dev.off()    
      
    }
  }


figure_4CD <- function(paths, out_path, indics=1:16) { 
  
  figs_path <- sprintf("%s\\%s", out_path, "figure_4")
  dir.create(figs_path)
  print(figs_path)
  true_data_paths <- sprintf("%s\\%s", paths, "equalized")
  
  estimation_list <- list()
  likelihood_list <- list()
  session_idx <- c()
  
  for (working_path in true_data_paths) {
    
    # extract all paths
    session_paths <- 
      list.dirs(working_path, recursive=F)[grep("session", list.dirs(working_path, recursive=F))]
    
    # Run through all paths
    for (tpath in session_paths) {
      
      idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
      if (sum(idx == indices) == 0){
        next
      }
      
      if (len((grep("pct_estimate_df.R", list.files(sprintf("%s\\simulations_2", tpath))))) == 0) {
        next
      }
      
      
      
      load(sprintf("%s\\simulations_2\\%s", tpath,"pct_estimate_df.R"))

      
      likelihood <- final_df[,ncol(final_df)] # Put last column on different vector
      final_df <- final_df[,-ncol(final_df)] # remove last column
      
      percents_measured <- c()
      
      for (i in 1:nrow(final_df)) {
        
        kdf_func <- approxfun(density(final_df[i,], width=0.10, from=0, to=1))
        measured_p1 <- sapply(seq(0,0.999, by=0.001), function(v) {euclidean(integrate(kdf_func, v, 1, stop.on.error = F)$value, likelihood[i])})
        measured_p2 <- sapply(seq(0.001,1,  by=0.001), function(v) {euclidean(integrate(kdf_func, 0, v, stop.on.error = F)$value, likelihood[i])})
        
        measured_p1 <- seq(0,0.999, by=0.001)[which.min(measured_p1)]
        measured_p2 <- seq(0.001,1,  by=0.001)[which.min(measured_p2)]
        
        percents_measured <- rbind(percents_measured, c(measured_p1, measured_p2))
      }
      
      session_idx <- c(session_idx, idx)
      estimation_list <- append(estimation_list, list(as.numeric(names(which.max(table(percents_measured))))))
      likelihood_list <- append(likelihood_list, list(as.numeric(names(which.max(likelihood)))))

    }
  }
  
  percentage_df <- cbind(unlist(estimation_list), unlist(likelihood_list), session_idx)
  colnames(percentage_df) <- c("Measured", "Estimated", "Session")
  percentage_df[percentage_df[,"Session"] > 8, "Session"] <- percentage_df[percentage_df[,"Session"] > 8, "Session"] - 8
  
  percentage_df <- percentage_df[percentage_df[,"Session"] %in% c(1,2,7,8),]
  percentage_df <- data.frame(percentage_df)
  
  summary_df <- 
    ddply(percentage_df, .(Session), function(sub_df) {c(nrow(sub_df), apply(sub_df[,1:2], 2, sd), colMeans(sub_df))})
  
  percentage_df$Session <- factor(as.character(percentage_df[,"Session"]), levels=c("1", "2", "7", "8"))
  summary_df$Session <- factor(as.character(summary_df[,"Session"]), levels=c("1", "2", "7", "8"))
  
  
  
  
  colnames(summary_df) <- c("N", "sd_measured", "sd_est", "Measured", "Estimated", "Session")
  g_est <- ggplot(summary_df, aes(y=Estimated, x=Session)) + 
       geom_bar(stat = "identity", fill="darkgoldenrod2", color="black") + 
       geom_jitter(data=percentage_df, position = position_jitter(0.1),color = "black", size=1.5) +
       geom_errorbar(aes(ymin = Estimated-sd_est, ymax = Estimated+sd_est), width=0.5, size=1) +
       theme_light() +
       theme(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank())
  
  g_measured <- ggplot(summary_df, aes(y=Measured, x=Session)) + 
    geom_bar(stat = "identity", fill="darkgoldenrod2", color="black") + 
    geom_jitter(data=percentage_df, position = position_jitter(0.1),color = "black", size=1.5) +
    geom_errorbar(aes(ymin = Measured-sd_measured, ymax = Measured+sd_measured), width=0.5, size=1) +
    theme_light() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  
  pdf(file=sprintf("%s\\%d_4c.pdf", 
                   figs_path, 
                   idx),
      width=5,
      height=5)
  plot(g_est)
  dev.off()
  
  pdf(file=sprintf("%s\\%d_4d.pdf", 
                   figs_path, 
                   idx),
      width=5,
      height=5)
  plot(g_measured)
  dev.off()

}


figure_4AB <- function(paths, out_path, indics=1:16) { 
  
  figs_path <- sprintf("%s\\%s", out_path, "figure_4")
  dir.create(figs_path)
  print(figs_path)
  true_data_paths <- sprintf("%s\\%s", paths, "equalized")
  

  
  for (working_path in true_data_paths) {
    estimation_list <- list()
    likelihood_list <- list()
    session_idx <- c()
    # extract all paths
    session_paths <- 
      list.dirs(working_path, recursive=F)[grep("session", list.dirs(working_path, recursive=F))]
    s = 0
    # Run through all paths
    for (tpath in session_paths) {
      
      idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
      if (sum(idx == indices) == 0){
        next
      } else {
        s <- s + 1
      }
      
      if (len((grep("pct_estimate_df.R", list.files(sprintf("%s\\simulations_2", tpath))))) == 0) {
        next
      }
      
      
      
      load(sprintf("%s\\simulations_2\\%s", tpath,"pct_estimate_df.R"))
      likelihood <- final_df[,ncol(final_df)] # Put last column on different vector
      session_idx <- c(session_idx, idx)
      likelihood_list <- append(likelihood_list, list(likelihood))
      
    }
    
    if (s == 16) {
      break
    }
  }
  
  mdf <- data.frame(session=rep(session_idx, each=7), 
                    percentage=rep(as.numeric(names(likelihood_list[[1]])), times=16),
                    likelihood=unlist(likelihood_list))
  
  # gdensity <- 
  #   ggplot(mdf,  aes(y = Var1, x = value)) + 
  #   geom_density_ridges(bandwidth=0.02, aes(fill=Var1), size=1, alpha=0.8, scale=0.9) + 
  #   coord_flip() + 
  #   scale_fill_brewer(palette="Spectral") +
  #   theme_light() +
  #   ylab("True simulation percentage (%)") +
  #   xlab("Measured simulation percentage (%)") +
  #   theme(panel.grid.major = element_blank(), panel.grid.minor= element_blank(), legend.position =  "NA")  +
  #   geom_vline(xintercept=place_cell_percentage, linetype="dashed", col="gray60", size=2) +
  #   geom_text(y=place_cell_percentage + 0.01, x=2.5, label="Measure place cell percentage (actual data)",col="gray60") 
}



plot_neur_raster <- function(stim_trace, spike_train, idx, coltoc="red", cx) {
  plot(stim_trace, type="l", col="gray60")
  points(which(spike_train[idx,] != 0 ), stim_trace[which(spike_train[idx,] != 0 )], col=coltoc,
         pch=19, cex=cx)
}


compute_MI <- function(stim_trace, spike_train) {
  eps <- 10^-30
  spatial_bins <- spatial_bins_t
  rounded_spike_train <- round(spike_train)
  
  MI <- lapply(1:nrow(rounded_spike_train),
                       
                       function (n_indx) {
                        neur <- rounded_spike_train[n_indx,]
                        spikes <- c(0,1)
                        #spikes <- sort(unique(neur))
                        
                        joint_prob_mat <- matrix(rep(0, times=len(spikes) * len(spatial_bins)), 
                                            nrow=len(spatial_bins))
                        rownames(joint_prob_mat) <- spatial_bins
                        colnames(joint_prob_mat) <- spikes
                        
                        for (sbin in spatial_bins) {
                          ind <- which(stim_trace == sbin)
                          # for (s_count in spikes) {
                          #   s_count_idx <- which(s_count == spikes)
                          #   joint_prob_mat[sbin, s_count_idx] <- (sum(neur[ind] == s_count))
                          # }
                          
                          joint_prob_mat[sbin, 1] <- (sum(neur[ind] == 0))
                          joint_prob_mat[sbin, 2] <-  len(neur[ind]) - joint_prob_mat[sbin, 1]
                        }
                        
                        sm <- sum(joint_prob_mat)
                        spike_prob <- colSums(joint_prob_mat) / sm
                        bin_prob <- rowSums(joint_prob_mat) / sm 
                        
                        MI <- c()
                        
                        
                        for (sbin in spatial_bins) {
                          for (s_count in spikes) {
                            s_count_idx <- which(s_count == spikes)
                            pso <- joint_prob_mat[sbin, s_count_idx] / sm
                            MI <- c(MI,
                                  pso * log((pso/(spike_prob[s_count_idx] * bin_prob[sbin])) + eps))
                            
                          }
                        }
                        
                        return(sum(MI))
                      })
  
  return(unlist(MI))

}


# compute_pois_product <- function(tuning_curves, firing_rate) {
#   
#   
#   pois_product <- 
#     unlist(lapply(1:nrow(tuning_curves),
#            function(n_indx) {
#              prod(ppois(tuning_curves[n_indx,], lambda=firing_rate[n_indx]))
#            }))
#   
#   return(pois_product)
# }


compute_pois_product <- function(spike_train, stim_trace) {

rm <- rowMeans(spike_train)
rounded <- round(spike_train)
for (sb in sort(unique(stim_trace))) {
  
  ind <- which(stim_trace == sb)
  rs <- rounded[,ind]
  res <- cbind(res, 
               unlist(lapply(1:nrow(spike_train), function(idx) {sum(ppois(rs[idx,], lambda=mean(rs[idx,])))})))
}

  return(apply(res, 1, sum))
}

out_path_f <- "C:\\Users\\itayta.WISMAIN\\Desktop\\allcells\\paper_figures"

figure_3 <- function(path, out_path, spike_train, stim_trace) {
  
  figs_path <- sprintf("%s\\%s", out_path, "figure_3")
  dir.create(figs_path)
  print(figs_path)
  true_path <- sprintf("%s\\%s", path, "equalized")
  
  
  session = 7
  likelihood_path_to_use <- "simulations_3"
  param_estimation_path_to_use <- "simulations_2"
  comparsion_path_to_use <- "sim_vs_real"
  
  fit_params_df <- extract_simulation_params(sprintf("%s\\equalized\\session_%d", path, session), plot=F,  absolute_min = T)
  params <- fit_params_df["1",]
  graphs <- plot_gen_place_cell_histograms(average_width = params["average_width"],
                                           sd_width = params["sd_width"],
                                           noise_level = params["noise"],
                                           double_peak_percent = params["double_peak_pct"],
                                           true_spike_train = spike_train,
                                           stim_trace = stim_trace,
                                           place_cell_percentage = 1)
  
  graphs$nrow = 1
  first_row <- do.call(grid.arrange, graphs)
  load(file=sprintf("%s\\session_%d\\%s\\pct_estimate_df.R", true_path,session, likelihood_path_to_use))
  
  likelihood_vec <- final_df[,"likelihood"]
  likelihood_df <- final_df[,-which(colnames(final_df) == "likelihood")] # Remove last
  
  percents_measured <- c()
  
  for (i in 1:nrow(final_df)) {
    
    kdf_func <- approxfun(density(final_df[i,], width=0.10, from=0, to=1))
    measured_p1 <- sapply(seq(0,0.999, by=0.001), function(v) {euclidean(integrate(kdf_func, v, 1, stop.on.error = F)$value, likelihood_vec[i])})
    measured_p2 <- sapply(seq(0.001,1,  by=0.001), function(v) {euclidean(integrate(kdf_func, 0, v, stop.on.error = F)$value, likelihood_vec[i])})
    
    measured_p1 <- seq(0,0.999, by=0.001)[which.min(measured_p1)]
    measured_p2 <- seq(0.001,1,  by=0.001)[which.min(measured_p2)]
    
    percents_measured <- rbind(percents_measured, c(measured_p1, measured_p2))
    
    
  }
  
  percents_tab <- table(percents_measured)
  percents_tab <- percents_tab[as.numeric(names(percents_tab)) > 0.01]
  measured_percent <- as.numeric(names(which.max(percents_tab)))
  
  load(file=sprintf("%s\\session_%d\\%s\\sim_vs_real_by_sample.R", true_path, session, comparsion_path_to_use))
  simulation_vs_real_df <- final_result
  
  
  mdf <- melt(likelihood_df[-1,])
  mdf$Var1 <- factor(mdf$Var1)
  
  bar_df <- data.frame(Likelihood=likelihood_vec[-1], 
                       Percentage=factor(sprintf("%d%%",
                                                 as.numeric(names(likelihood_vec[-1])) * 100),
                                         levels=sprintf("%d%%",
                                                        seq(40,100,by=10))))
  
  limits_of_plots <- c(min(mdf$value) - 0.05, max(mdf$value) + 0.1)
  gbox <- 
    ggplot(mdf) +
    geom_boxplot(aes(y=value, group=Var1, fill=Var1, x=Var1), size=1, alpha=0.9) +
    
    scale_fill_brewer(palette="Spectral") +
    theme_light() +
    xlab("True simulation percentage (%)") +
    ylab("Measured simulation percentage (%)") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position="NA",
          panel.border = element_blank(),
          panel.background = element_blank()) +
    geom_hline(yintercept=measured_percent, linetype="dashed", col="gray60", size=2) +
    geom_text(y=measured_percent + 0.01, x=2.5, label="Measured place cell percentage (actual data)",col="gray60",size=3) 
  
  gbar <- ggplot(bar_df, aes(x=Percentage, y=Likelihood)) +
    geom_bar(stat="identity", fill="gray50", color="black") + 
    theme_light() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position="NA",
          panel.border = element_blank(),
          panel.background = element_blank())
  
  gdensity <- 
    ggplot(mdf,  aes(y = Var1, x = value)) + 
    geom_density_ridges(bandwidth=0.02, aes(fill=Var1), size=1, alpha=0.8, scale=0.9) + 
    coord_flip() + 
    scale_fill_brewer(palette="Spectral") +
    theme_light() +
    ylab("True simulation percentage (%)") +
    xlab("Measured simulation percentage (%)") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position="NA",
          panel.border = element_blank(),
          panel.background = element_blank()) +
    geom_vline(xintercept=measured_percent, linetype="dashed", col="gray60", size=2) +
    geom_text(y=measured_percent + 0.01, x=2.5, label="Measure place cell percentage (actual data)",col="gray60") 
  
  
  second_row <- grid.arrange(gbox, gdensity, gbar, nrow=1)
  
  
  g_real_cyc <-  ggplot(simulation_vs_real_df$real, aes(x=Time, y=Cyclic)) + geom_line(size=2) + 
    theme_light() +
    ylim(c(0, 1)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position="NA",
          panel.border = element_blank(),
          panel.background = element_blank()) +
    xlim(0, max(simulation_vs_real_df$real$Time)) + 
    labs(x="Sample duration (minutes)", y = "Fraction of place cells (%)")
  
  
  sim_pct_colors = brewer.pal(n = len(levels(mdf$Var1)),  name = "Spectral")
  names(sim_pct_colors) <- as.numeric(levels(mdf$Var1))
  pval_df <- c()
  plots_list <- list()
  
  for (sim_pct in as.numeric(levels(mdf$Var1))) {
    print(sim_pct)
    simulation_df <- final_result[[as.character(sim_pct)]]
    g_sim_cyc <- g_real_cyc + ggtitle(sprintf("%d%% place cells simulation", as.numeric(sim_pct) * 100))
    pval_vec <- c()
    
    for (i in 1:(ncol(simulation_df) / 2)) {
      tmp_df_cyc <- data.frame(Fraction=simulation_df[, (i * 2)],
                               Time=simulation_vs_real_df$real$Time)
     
      g_sim_cyc <- 
        g_sim_cyc + geom_line(data=tmp_df_cyc,
                              aes(x=Time, y=Fraction),
                              color=sim_pct_colors[as.character(sim_pct)],
                              size=2,
                              alpha=0.5)
      
      
      pval_vec <- c(pval_vec,
                    wilcox.test(simulation_df[, (i * 2)],
                                simulation_vs_real_df$real$Cyclic)$p.value)
      
    }
    
    pval_df <- rbind(pval_df, pval_vec)
    plots_list <- append(plots_list, list(g_sim_cyc))
  }
  
  plots_list$nrow = 1
  third_row <- do.call(grid.arrange, plots_list)
  
  aligned_rows <- align_plots(first_row, second_row, third_row)
  final_aligned = grid.arrange(aligned_rows[[1]],
                               aligned_rows[[2]],
                               aligned_rows[[3]])
  
  plot_w = 15
  pdf(file=sprintf("%s\\figure_3.pdf", figs_path),
      height=(plot_w/4 + plot_w/3 + plot_w/6),
      width=plot_w)
  plot(final_aligned)
  dev.off()
  
 
  pdf(file=sprintf("%s\\figure_3_first_row.pdf", figs_path),
      height=(plot_w/4),
      width=plot_w)
  plot(aligned_rows[[1]])
  dev.off() 
  
  
  pdf(file=sprintf("%s\\figure_3_second_row.pdf", figs_path),
      height=(plot_w/3),
      width=plot_w)
  plot(aligned_rows[[2]])
  dev.off() 
  
  pdf(file=sprintf("%s\\figure_3_third_row.pdf", figs_path),
      height=(plot_w/6),
      width=plot_w)
  plot(aligned_rows[[3]])
  dev.off() 
}

figure_4_6 <- function(paths, out_path,  fig_ext = "4", path_region_str_split = "CA1\\\\", ext="") { 
  
  figs_path <- sprintf("%s\\%s_%s", out_path, "figure", fig_ext)
  dir.create(figs_path)
  print(figs_path)
  true_data_paths <- sprintf("%s\\%s", paths, "equalized")
  
  sessions_to_use=c(5:8,13:16)
  
  unique_paths_both <- 
    unique(unlist(lapply(str_split(true_data_paths, "Matrices"), function(sr) {return(sr[[1]])})))
  
  percent_df_list <- lapply(unique_paths_both, 
                             function(p) {
                               load(sprintf("%s\\percent_df.R", p))
                               return(final_df_list)
                             })
    
  all_df <- c()
  
  for(i in 1:len(percent_df_list)) {
    measured_percent_df <- percent_df_list[[i]] 
    measurements_df <- c()
    
    for (col_idx in 1:ncol(measured_percent_df[[1]])) {
      
      
      col_sessions <- 
      do.call(rbind,
              lapply(measured_percent_df,
                     function(df_rep) {
                       return(df_rep[sessions_to_use, col_idx])
                     }))
      
      measurements_df <- cbind(measurements_df,
                              colMeans(col_sessions))
    }
    
    measurements_df <- cbind(measurements_df,
                             sessions_to_use)
    
    measurements_df <- as.data.frame(measurements_df)
    
    measurements_df <- cbind(measurements_df,
                             rep(str_split(str_split(unique_paths_both[[i]], path_region_str_split)[[1]][2], "\\\\")[[1]][1],
                                        times=len(sessions_to_use)))
    
    all_df <- rbind(all_df,
                    measurements_df)
  }
   
  colnames(all_df) <- c(colnames(measured_percent_df[[1]]),
                        "Session", "Mice")
  
  
  percent_active_df <- all_df[,c(7:9)]
  percent_place_cells_df <- all_df[,c(1:3)]
  percent_all_df <- all_df[,c(4:6)]
  
  colnames(percent_active_df) <- c("Both", "Right", "Left")
  colnames(percent_place_cells_df) <- c("Both", "Right", "Left")
  colnames(percent_all_df) <- c("Both", "Right", "Left")
  
  melted_percent_active <- melt(percent_active_df)
  melted_percent_place_cells <- melt(percent_place_cells_df)
  melted_percent_all <- melt(percent_all_df)
  
  ylabs = c("Percentage of active cells (%)",
            "Percentage of place cells (% of active cells)",
            "Percentage of place cells (% of all cells)")
  
  dfs <- list(melted_percent_active,
           melted_percent_place_cells,
           melted_percent_all)
  
  bar_plots_list <- list()
  for (df_idx in 1:3) {
    df_to_use <- dfs[[df_idx]]
    colnames(df_to_use) <- c("Direction", "Percent")
    sd_df <- ddply(df_to_use, .(Direction), function(sub_df) {return(c(mean(sub_df[,"Percent"]), sd(sub_df[,"Percent"])))})
    colnames(sd_df) <- c("Direction", "Percent", "Sd")
    gr <- 
    ggplot(sd_df, aes(x=Direction, y=Percent)) + 
      geom_bar(stat="summary", aes(group=Direction), width=0.5,
               color="NA", fill="royalblue4", position = "dodge") + 
      # geom_errorbar(aes(ymin=Percent, ymax=Percent + Sd, group=Direction),
      #               size=2, width=0.3) +
      geom_jitter(data=df_to_use,
                 aes(x=Direction, y=Percent, group=Direction), position=position_jitter(0.2), color="gray70", size=3, alpha=0.7) + 
      theme_light() +     
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_blank(),
            panel.background = element_blank()) + 
      xlab("Running direction") +
      ylab(ylabs[df_idx]) + 
      ylim(0,1)
    
    bar_plots_list <- append(bar_plots_list,
                            list(gr))
  }
  
  bar_plots_list$nrow <- 1
  first_row <- do.call(grid.arrange, bar_plots_list)
  
  estimation_list <- list()
  likelihood_list <- list()
  measured_percent_list <- list()
  estimated_per_percent_list <- list()
  session_metavar_list <- list()
  
  for (working_path in true_data_paths) {
    
    # extract all paths
    session_paths <- 
      list.dirs(working_path, recursive=F)[grep("session", list.dirs(working_path, recursive=F))]
    
    # Run through all paths
    for (tpath in session_paths) {
      
      idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
      if (sum(idx == sessions_to_use) == 0){
        next
      }
      
      if (len((grep("pct_estimate_df.R", list.files(sprintf("%s\\simulations_3", tpath))))) == 0) {
        print(tpath)
        print("Missing")
        next
      }
      
      
      
      load(sprintf("%s\\simulations_3\\%s", tpath,"pct_estimate_df.R"))
      
      
      likelihood <- final_df[,ncol(final_df)] # Put last column on different vector
      final_df <- final_df[,-ncol(final_df)] # remove last column
      
      percents_measured <- c()
      
      for (i in 1:nrow(final_df)) {
        
        kdf_func <- approxfun(density(final_df[i,], width=0.10, from=0, to=1))
        measured_p1 <- sapply(seq(0,0.999, by=0.001), function(v) {euclidean(integrate(kdf_func, v, 1, stop.on.error = F)$value, likelihood[i])})
        measured_p2 <- sapply(seq(0.001,1,  by=0.001), function(v) {euclidean(integrate(kdf_func, 0, v, stop.on.error = F)$value, likelihood[i])})
        
        measured_p1 <- seq(0,0.999, by=0.001)[which.min(measured_p1)]
        measured_p2 <- seq(0.001,1,  by=0.001)[which.min(measured_p2)]
        
        percents_measured <- rbind(percents_measured, c(measured_p1, measured_p2))
      }
      
      percents_tab <- table(percents_measured)
      percents_tab <- percents_tab[as.numeric(names(percents_tab)) > 0.01]
      
      measured_percent_list <- append(measured_percent_list, list(as.numeric(names(which.max(percents_tab)))))
      estimation_list <- append(estimation_list, list(as.numeric(names(which.max(likelihood)))))
      likelihood_list <- append(likelihood_list, list(likelihood))
      estimated_per_percent_list <- append(estimated_per_percent_list, list(rowMeans(final_df)))
      
      
      
      mice_str <- str_split(str_split(tpath, path_region_str_split)[[1]][2], "\\\\Matrices")[[1]]
      
      session_metavar_list <- append(session_metavar_list, list(c(idx, mice_str[1], ifelse(grepl("Left", mice_str[2]), "Left", "Right"))))
      
    }
  }
  
  likelihood_df <- do.call(rbind, likelihood_list)
  estimated_per_percent_df <- do.call(rbind, estimated_per_percent_list)
  session_df <- do.call(rbind, session_metavar_list)
  
  estimation_vec <- unlist(estimation_list)
  measured_vec <- unlist(measured_percent_list)
  
  estimated_df <- data.frame(percent_estimated=estimation_vec)
  estimated_df <- cbind(estimated_df, session_df)
  colnames(estimated_df) <- c("Estimated", "Session", "Mice", "Direction")
  
  # Create weighted estimations of % of place cells for both running directions
  both_dir_df <- ddply(estimated_df, .(Mice), 
                 function(sub_df) {
                    left <- sub_df[which(sub_df[,"Direction"] == "Left"),]
                    right <- sub_df[which(sub_df[,"Direction"] == "Right"),]
                    
                    intersecting_sessions <- as.numeric(intersect(left[,"Session"], right[,"Session"]))
                    
                    mice = sub_df[1,"Mice"]
                    
                    for (ses in intersecting_sessions) {
                      left_ses <- left[which(as.numeric(left[,"Session"]) == ses),]
                      right_ses <- right[which(as.numeric(right[,"Session"]) == ses),]
                      
                      # Get active percentage for weighted averages
                      tmp_all_df_row <- all_df[which(all_df[,"Session"] == ses & all_df[,"Mice"] == mice),]
                      
                      right_p <- tmp_all_df_row[,"Right active"]
                      left_p <- tmp_all_df_row[,"Left active"]

                      
                      weighted_estimate_both_dir <- 
                        (right_ses[,"Estimated"] * right_p + left_ses[,"Estimated"] * left_p) / (left_p + right_p)
                      
                      sub_df <- rbind(sub_df,
                                      c(weighted_estimate_both_dir, ses, mice, "Both"))
                    }
                    
                    return(sub_df)
                   })
  
  estimated_per_percent_df <- cbind(estimated_per_percent_df, measured_vec)
  estimated_per_percent_df <- as.data.frame(estimated_per_percent_df)
  estimated_per_percent_df <- cbind(estimated_per_percent_df, session_df)
  
  
  
  melted_likelihood <- melt(t(apply(likelihood_df[,-1], 1, function(r) {return(r/sum(r))})))
  colnames(melted_likelihood) <- c("#", "Simulated percent", "Normalized likelihood")
  
  both_dir_df <- as.data.frame(both_dir_df)
  both_dir_df[,"Estimated"] <- as.numeric(both_dir_df[,"Estimated"])
  mean_sd_df <- ddply(both_dir_df, .(Direction), 
                       function(sub_df) {
                         return(c(mean(sub_df[,"Estimated"]),
                                  sd(sub_df[,"Estimated"])))
                       })
  
  colnames(mean_sd_df) <- c("Direction", "Estimated", "Sd")
  gest <- 
    ggplot(mean_sd_df, aes(x=Direction, y=Estimated)) + 
    geom_bar(stat="summary", aes(group=Direction), color="NA", 
             width=0.5, fill="gray70", position = "dodge") + 
    # geom_errorbar(aes(ymin=Estimated, ymax=Estimated + Sd, group=Direction),
    #                   size=2, width=0.3) +
    geom_jitter(data=both_dir_df,
      aes(x=Direction, y=Estimated, group=Direction), position=position_jitter(0.2), color="gray40", size=3, alpha=0.5) + 
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank()) + 
    xlab("Running direction") +
    ylab(ylabs[2]) + 
    ylim(0,1.05)
  
  
  likelihood_mean_sd_df <- ddply(melted_likelihood, .(`Simulated percent`), 
                      function(sub_df) {
                        return(c(mean(sub_df[,"Normalized likelihood"]),
                                 sd(sub_df[,"Normalized likelihood"])))
                      })
  
  colnames(likelihood_mean_sd_df) <- c("Simulated percent", "Mean", "Sd")
  
  melted_likelihood$`Simulated percent` <- factor(melted_likelihood$`Simulated percent` * 100,
                                                  levels=as.character(likelihood_mean_sd_df$`Simulated percent`  * 100))
  likelihood_mean_sd_df$`Simulated percent` <- factor(likelihood_mean_sd_df$`Simulated percent` * 100, 
                                                      levels=as.character(likelihood_mean_sd_df$`Simulated percent` * 100))
  
  glikelihood <- 
    ggplot(likelihood_mean_sd_df, aes(x=`Simulated percent`, y=Mean)) + 
    geom_bar(stat="summary", aes(group=`Simulated percent`), color="NA", 
             width=0.5, fill="gray70", position = "dodge") + 
    # geom_errorbar(aes(ymin=Mean, ymax=Mean + Sd, group=`Simulated percent`),
    #               size=2, width=0.3) +
    geom_jitter(data=melted_likelihood,
                aes(x=`Simulated percent`, y=`Normalized likelihood`, group=`Simulated percent`), 
                position=position_jitter(0.05), color="gray40", size=3, alpha=0.5) + 
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank()) + 
    xlab("Simulated percent (%)") +
    ylab("Normalized likelihood") + 
    ylim(0,1.05)
  
  
  scatter_list <- list()
  for (est_p in as.character(seq(0.5, 1, by=0.1))) {
   
    tmp_df <- data.frame(Measured=estimated_per_percent_df[,"measured_vec"],
                         Estimated=estimated_per_percent_df[,est_p],
                         Session=as.numeric(estimated_per_percent_df[,ncol(estimated_per_percent_df) - 2]))
    
    tmp_df$Session[tmp_df$Session > 8] <- tmp_df$Session[tmp_df$Session > 8] - 8
    tmp_df$Session <- factor(sprintf("%dth", tmp_df$Session), c("5th", "6th", "7th", "8th"))
    
    linear_line_df <- data.frame(x=c(0.3,1),
                                 y=c(0.3,1))
    gs <- 
    ggplot(linear_line_df) +
      geom_line(aes(x=x,y=y), linetype="longdash", color="#c0c1c3",
                size=1.5, alpha=0.5) +    
      geom_point(data=tmp_df, aes(x=Measured, y=Estimated, color=Session), size=3) +
      theme_light() +     
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position="NA") +
    scale_color_manual(breaks=c("5th","6th","7th","8th"),
                       values=c("#3f9ab4",
                                "#6cbfa6",
                                "#b0c7e0",
                                "#8a83ad")) + 
      ylim(c(0.3,1)) + 
      xlim(c(0.3,1)) + 
      ylab("") +
      xlab("") + 
      ggtitle(sprintf("r = %f",
              cor(tmp_df$Measured,
                  tmp_df$Estimated)))
    
     scatter_list <- append(scatter_list, list(gs))
  }
  
  scatter_list$nrow = 2
  scatter_list$left="Estimated place cell percentage (%)"
  scatter_list$bottom="Measured place cell percentage (%)"
  scattersp <- do.call(grid.arrange, scatter_list)
  plot_h <- 8
  
  second_row <- grid.arrange(gest, glikelihood, align_plots(scattersp)[[1]], 
                             nrow=1,widths=c(1,1,1.5))
  aligned <- align_plots(first_row, second_row, align="v")
  
  

  
  pdf(file=sprintf("%s\\figure_%s_first_row%s.pdf", figs_path, fig_ext, ext),
      height=(plot_h * 3.5 / 3),
      width=(plot_h * 3.5))
  plot(first_row)
  dev.off() 
    
  pdf(file=sprintf("%s\\figure_%s_second_row%s.pdf", figs_path, fig_ext, ext),
      height=plot_h,
      width=(plot_h * 3.5))
  
  plot(aligned[[2]])
  dev.off() 

  pdf(file=sprintf("%s\\figure_%s_scatter%s.pdf", figs_path, fig_ext, ext),
      height=plot_h,
      width=plot_h * 1.5)
  plot(scattersp)
  dev.off()   
  
  pdf(file=sprintf("%s\\figure_%s_estimate%s.pdf", figs_path, fig_ext, ext),
      height=plot_h,
      width=plot_h)
  plot(gest)
  dev.off()   
  
  pdf(file=sprintf("%s\\figure_%s_likelihood%s.pdf", figs_path, fig_ext, ext),
      height=plot_h,
      width=plot_h)
  plot(glikelihood)
  dev.off()   
  
  
  pdf(file=sprintf("%s\\figure_%s_pct_active%s.pdf", figs_path, fig_ext, ext),
      height=plot_h,
      width=plot_h)
  plot(bar_plots_list[[1]])
  dev.off()   
  
  pdf(file=sprintf("%s\\figure_%s_pct_active_place_cells%s.pdf", figs_path, fig_ext, ext),
      height=plot_h,
      width=plot_h)
  plot(bar_plots_list[[2]])
  dev.off()   
  
  pdf(file=sprintf("%s\\figure_%s_pct_place_cells_of_all%s.pdf", figs_path, fig_ext, ext),
      height=plot_h,
      width=plot_h)
  plot(bar_plots_list[[3]])
  dev.off()   
  
  pdf(file=sprintf("%s\\figure_%s%s.pdf", figs_path, fig_ext, ext),
      height=(plot_h + plot_h * 3.5 / 3),
      width=(plot_h * 3.5))
  plot(grid.arrange(aligned[[1]], aligned[[2]], nrow=2))
  dev.off()  
  
  
}


figure_5_7 <- function(paths, 
                       out_path, 
                       fig_ext = "5", 
                       path_region_str_split = "CA1\\\\", 
                       ext="",
                       non_sim_cp=c(Both="red4", Right="red3", Left="red"),
                       ca3=F) { 
  
  
  if (ca3) {
    fig_ext = "7"
    path_region_str_split = "CA3\\\\"
    non_sim_cp = ca3_non_sim_cp=c(Both="goldenrod4", Right="goldenrod3", Left="goldenrod")
  }
  
  figs_path <- sprintf("%s\\%s_%s", out_path, "figure", fig_ext)
  dir.create(figs_path)
  print(figs_path)
  true_data_paths <- sprintf("%s\\%s", paths, "equalized")
  
  sessions_to_use=c(1:16)
  
  unique_paths_both <- 
    unique(unlist(lapply(str_split(true_data_paths, "Matrices"), function(sr) {return(sr[[1]])})))
  
  percent_df_list <- lapply(unique_paths_both, 
                            function(p) {
                              load(sprintf("%s\\percent_df_both_directions.R", p))
                              return(final_df_list)
                            })
  
  all_df <- c()
  
  for(i in 1:len(percent_df_list)) {
    measured_percent_df <- percent_df_list[[i]] 
    measurements_df <- c()
    
    for (col_idx in 1:ncol(measured_percent_df[[1]])) {
      
      
      col_sessions <- 
        do.call(rbind,
                lapply(measured_percent_df,
                       function(df_rep) {
                         return(df_rep[sessions_to_use, col_idx])
                       }))
      
      measurements_df <- cbind(measurements_df,
                               colMeans(col_sessions))
    }
    
    measurements_df <- cbind(measurements_df,
                             sessions_to_use)
    
    measurements_df <- as.data.frame(measurements_df)
    
    measurements_df <- cbind(measurements_df,
                             rep(str_split(str_split(unique_paths_both[[i]], path_region_str_split)[[1]][2], "\\\\")[[1]][1],
                                 times=len(sessions_to_use)))
    
    all_df <- rbind(all_df,
                    measurements_df)
  }
  
  colnames(all_df) <- c(colnames(measured_percent_df[[1]]),
                        "Session", "Mice")
  
  
  percent_active_df <- all_df[,c(7:9)]
  percent_place_cells_df <- all_df[,c(1:3)]
  percent_all_df <- all_df[,c(4:6)]
  
  sessions_df <- data.frame("Both"=all_df[,10],
                            "Right"=all_df[,10],
                            "Left"=all_df[,10])
  
  colnames(percent_active_df) <- c("Both", "Right", "Left")
  colnames(percent_place_cells_df) <- c("Both", "Right", "Left")
  colnames(percent_all_df) <- c("Both", "Right", "Left")
  
  melted_sessions_df <- melt(sessions_df)
  melted_percent_active <- cbind(melt(percent_active_df), melted_sessions_df[,"value"])
  melted_percent_place_cells <- cbind(melt(percent_place_cells_df), melted_sessions_df[,"value"])
  melted_percent_all <- cbind(melt(percent_all_df), melted_sessions_df[,"value"])
  
  
  ylabs = c("Fraction of active cells (%)",
            "Fraction of place cells (% of active cells)",
            "Fraction of place cells (% of all cells)")
  
  dfs <- list(melted_percent_active,
              melted_percent_place_cells,
              melted_percent_all)
  
  bar_plots_list <- list()
  for (df_idx in 1:3) {
    df_to_use <- dfs[[df_idx]]
    colnames(df_to_use) <- c("Direction", "Percent", "Session")
    
    df_to_use$Session[df_to_use$Session > 8] <- df_to_use$Session[df_to_use$Session > 8] - 8 
    df_to_use$Session <- factor(df_to_use$Session, levels=as.character(1:8))
    
    
    sd_df <- ddply(df_to_use, .(Session), 
                   function(sub_df) {
                      sd_df <- ddply(sub_df, .(Direction), 
                                    function(sub_df_inner) 
                                      {return(c(mean(sub_df_inner[,"Percent"]), 
                                                sd(sub_df_inner[,"Percent"])))})
                      return(sd_df)
                     })
    colnames(sd_df) <- c("Session","Direction", "Percent", "Sd")
    gr <- 
      ggplot(sd_df, aes(x=Session, y=Percent, group=Direction)) + 
      geom_point(data=df_to_use, aes(x=Session, y=Percent, color=Direction, group=Direction), 
                  size=3.5, alpha=0.2,position=position_dodge(0.4)) + 
      geom_line(aes(group=Direction, color=Direction), size=1.5, position=position_dodge(0.4)) + 
      geom_errorbar(aes(ymin=Percent - Sd, ymax=Percent + Sd, group=Direction, color=Direction), 
                    size=2, width=0.3, position=position_dodge(0.4)) +
      geom_point(aes(group=Direction, color=Direction), size=5, position=position_dodge(0.4)) + 
      theme_light() +     
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_blank(),
            panel.background = element_blank()) + 
      xlab("Session") +
      ylab(ylabs[df_idx]) + 
      ylim(0.05,0.85) + 
      scale_color_manual(breaks=c("Both", "Right", "Left"),
                        values=c(non_sim_cp["Both"], non_sim_cp["Right"], non_sim_cp["Left"]))
    
    bar_plots_list <- append(bar_plots_list,
                             list(gr))
  }
  
  bar_plots_list$nrow <- 1
  first_row <- do.call(grid.arrange, bar_plots_list)
  
  estimation_list <- list()
  likelihood_list <- list()
  measured_percent_list <- list()
  estimated_per_percent_list <- list()
  session_metavar_list <- list()
  
  for (working_path in true_data_paths) {
    
    # extract all paths
    session_paths <- 
      list.dirs(working_path, recursive=F)[grep("session", list.dirs(working_path, recursive=F))]
    
    # Run through all paths
    for (tpath in session_paths) {
      
      idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
      if (sum(idx == sessions_to_use) == 0){
        next
      }
      
      if (len((grep("pct_estimate_df.R", list.files(sprintf("%s\\simulations_3", tpath))))) == 0) {
        print(tpath)
        print("Missing")
        next
      }
      
      
      
      load(sprintf("%s\\simulations_3\\%s", tpath,"pct_estimate_df.R"))
      
      
      likelihood <- final_df[,ncol(final_df)] # Put last column on different vector
      final_df <- final_df[,-ncol(final_df)] # remove last column
      
      percents_measured <- c()
      
      for (i in 1:nrow(final_df)) {
        
        kdf_func <- approxfun(density(final_df[i,], width=0.10, from=0, to=1))
        measured_p1 <- sapply(seq(0,0.999, by=0.001), function(v) {euclidean(integrate(kdf_func, v, 1, stop.on.error = F)$value, likelihood[i])})
        measured_p2 <- sapply(seq(0.001,1,  by=0.001), function(v) {euclidean(integrate(kdf_func, 0, v, stop.on.error = F)$value, likelihood[i])})
        
        measured_p1 <- seq(0,0.999, by=0.001)[which.min(measured_p1)]
        measured_p2 <- seq(0.001,1,  by=0.001)[which.min(measured_p2)]
        
        percents_measured <- rbind(percents_measured, c(measured_p1, measured_p2))
      }
      
      
      percents_tab <- table(percents_measured)
      percents_tab <- percents_tab[as.numeric(names(percents_tab)) > 0.01]
      
      measured_percent_list <- append(measured_percent_list, list(as.numeric(names(which.max(percents_tab)))))
      estimation_list <- append(estimation_list, list(as.numeric(names(which.max(likelihood)))))
      likelihood_list <- append(likelihood_list, list(likelihood))
      estimated_per_percent_list <- append(estimated_per_percent_list, list(rowMeans(final_df)))
      
      
      
      mice_str <- str_split(str_split(tpath, path_region_str_split)[[1]][2], "\\\\Matrices")[[1]]
      
      session_metavar_list <- append(session_metavar_list, list(c(idx, mice_str[1], ifelse(grepl("Left", mice_str[2]), "Left", "Right"))))
      
    }
  }
  
  likelihood_df <- do.call(rbind, likelihood_list)
  estimated_per_percent_df <- do.call(rbind, estimated_per_percent_list)
  session_df <- do.call(rbind, session_metavar_list)
  
  estimation_vec <- unlist(estimation_list)
  measured_vec <- unlist(measured_percent_list)
  
  estimated_df <- data.frame(percent_estimated=estimation_vec)
  estimated_df <- cbind(estimated_df, session_df)
  colnames(estimated_df) <- c("Estimated", "Session", "Mice", "Direction")
  estimated_df$Session <- as.numeric(estimated_df$Session)
  
  
  estimated_per_percent_df <- cbind(estimated_per_percent_df, measured_vec)
  estimated_per_percent_df <- as.data.frame(estimated_per_percent_df)
  estimated_per_percent_df <- cbind(estimated_per_percent_df, session_df)
  
  
  
  # Create weighted estimations of % of place cells for both running directions
  both_dir_df <- ddply(estimated_df, .(Mice), 
                       function(sub_df) {
                         left <- sub_df[which(sub_df[,"Direction"] == "Left"),]
                         right <- sub_df[which(sub_df[,"Direction"] == "Right"),]
                         
                         intersecting_sessions <- as.numeric(intersect(left[,"Session"], right[,"Session"]))
                         
                         mice = sub_df[1,"Mice"]
                         
                         for (ses in intersecting_sessions) {
                           left_ses <- left[which(as.numeric(left[,"Session"]) == ses),]
                           right_ses <- right[which(as.numeric(right[,"Session"]) == ses),]
                           
                           # Get active percentage for weighted averages
                           tmp_all_df_row <- all_df[which(all_df[,"Session"] == ses & all_df[,"Mice"] == mice),]
                           
                           right_p <- tmp_all_df_row[,"Right active"]
                           left_p <- tmp_all_df_row[,"Left active"]
                           
                           
                           weighted_estimate_both_dir <- 
                             (right_ses[,"Estimated"] * right_p + left_ses[,"Estimated"] * left_p) / (left_p + right_p)
                           
                           sub_df <- rbind(sub_df,
                                           c(weighted_estimate_both_dir, ses, mice, "Both"))
                         }
                         
                         return(sub_df)
                       })

  
  
  both_dir_df <- as.data.frame(both_dir_df)
  both_dir_df[,"Estimated"] <- as.numeric(both_dir_df[,"Estimated"])
  both_dir_df[,"Session"] <- as.numeric(both_dir_df[,"Session"])
  
  both_dir_df[which(both_dir_df[,"Session"] > 8),"Session"] <- both_dir_df[which(both_dir_df[,"Session"] > 8),"Session"] - 8
  
  both_dir_mean_sd_df <- 
      ddply(both_dir_df, .(Session), 
                 function(sub_df) {
                   sd_df <- ddply(sub_df, .(Direction), 
                                  function(sub_df_inner) 
                                  { est <- as.numeric(sub_df_inner[,"Estimated"])
                                    return(c(mean(est), sd(est)))})
                   return(sd_df)
                 })
  
  colnames(both_dir_mean_sd_df) <- c("Session", "Direction", "Estimated", "Sd")
  both_dir_mean_sd_df$Session <- factor(both_dir_mean_sd_df$Session, levels=as.character(1:8))
  both_dir_df$Session <- factor(both_dir_df$Session, levels=as.character(1:8))
  both_dir_mean_sd_df$Direction <- factor(both_dir_mean_sd_df$Direction, levels=c("Both", "Right", "Left"))
  gest <- 
    ggplot(both_dir_mean_sd_df, aes(x=Session, y=Estimated)) + 
    geom_point(data=both_dir_df, aes(x=Session, y=Estimated, group=Direction, color=Direction), 
               size=3.5, alpha=0.2, position=position_dodge(0.4)) +
    geom_line(aes(group=Direction, color=Direction), size=1.5, position=position_dodge(0.4)) +
    geom_errorbar(aes(ymin=Estimated - Sd, ymax=Estimated + Sd, group=Direction, color=Direction), 
                  size=2, width=0.3, position=position_dodge(0.4)) +
    geom_point(aes(group=Direction, color=Direction), size=5, position=position_dodge(0.4)) + 
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank()) + 
    xlab("Session") +
    ylab(ylabs[df_idx]) + 
    ylim(0,1.1) + 
    scale_color_manual(breaks=c("Both", "Right", "Left"),
                       values=c("royalblue4", "royalblue2", "royalblue3"))
  
  
  likelihood_df <- (t(apply(likelihood_df[,-1], 1, function(r) {return(r/sum(r))})))
  likelihood_df <- as.data.frame(likelihood_df)
  likelihood_sessions <- as.numeric(session_df[,1])
  likelihood_sessions[likelihood_sessions > 8] <-  likelihood_sessions[likelihood_sessions > 8] - 8
  likelihood_df <- cbind(likelihood_df, likelihood_sessions)
  colnames(likelihood_df) <- c(as.character(seq(0.5,1,by=0.1)), "sessions")
  melted_likelihood <- melt(likelihood_df, id.vars="sessions")
  colnames(melted_likelihood) <- c("Session", "Simulated percent", "Normalized likelihood")
  
  likelihood_mean_sd_df <- ddply(melted_likelihood, .(`Simulated percent`), 
                                 function(sub_df) {
                                  ddply(sub_df, .(`Session`), 
                                         function(sub_df_inner) {
                                                mli <- sub_df_inner[,"Normalized likelihood"]
                                                     return(c(mean(mli), sd(mli)))})
                                 })
  
  colnames(likelihood_mean_sd_df) <- c("Simulated percent", "Session", "Mean", "Sd")
  
  melted_likelihood$`Simulated percent` <- as.numeric(as.character(melted_likelihood$`Simulated percent`))
  melted_likelihood$`Simulated percent` <- factor(melted_likelihood$`Simulated percent` * 100,
                                                  as.character(seq(50,100, by=10)))
  
  likelihood_mean_sd_df$`Simulated percent` <- as.numeric(as.character(likelihood_mean_sd_df$`Simulated percent`))
  likelihood_mean_sd_df$`Simulated percent` <- factor(likelihood_mean_sd_df$`Simulated percent` * 100, 
                                                      as.character(seq(50,100, by=10)))
  
  glikelihood <- 
    ggplot(likelihood_mean_sd_df, aes(x=`Simulated percent`, y=Mean, group=Session)) + 
    geom_bar(stat="summary", aes(group=`Session`, fill=Session), color="NA", 
             position = "dodge") + 
    # geom_errorbar(aes(ymin=Mean, ymax=Mean + Sd, group=`Session`),
    #               size=2, width=0.3, position="dodge") +
    # geom_jitter(data=melted_likelihood,
    #             aes(x=`Simulated percent`, y=`Normalized likelihood`, group=`Session`, fill=Session), position=position_jitterdodge(0), color="gray40", size=3, alpha=0.5) + 
    scale_fill_distiller(palette="RdYlBu") +
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank()) + 
    xlab("Simulated percent (%)") +
    ylab("Normalized likelihood") + 
    ylim(0,max(likelihood_mean_sd_df[,"Mean"]) + 0.02) 


  
  
  plot_h <- 8
  
  first_row <- grid.arrange(bar_plots_list[[1]], bar_plots_list[[2]], nrow=1)
  second_row <- grid.arrange(bar_plots_list[[3]], gest, nrow=1)
  aligned <- align_plots(first_row, second_row, glikelihood, align="v")
  
  
  plot_h = 8
  
  pdf(file=sprintf("%s\\figure_%s_first_row%s.pdf", figs_path, fig_ext, ext),
      height=plot_h,
      width=plot_h * 3)
  plot(aligned[[1]])
  dev.off() 
  
  pdf(file=sprintf("%s\\figure_%s_second_row%s.pdf", figs_path, fig_ext, ext),
      height=plot_h,
      width=plot_h * 3)
  
  plot(aligned[[2]])
  dev.off() 
  
  
  pdf(file=sprintf("%s\\figure_%s_estimate%s.pdf", figs_path,fig_ext, ext),
      height=plot_h,
      width=plot_h * 1.5)
  plot(gest)
  dev.off()   
  
  pdf(file=sprintf("%s\\figure_%s_likelihood%s.pdf", figs_path, fig_ext, ext),
      height=plot_h,
      width=plot_h * 3)
  plot(aligned[[3]])
  dev.off()   
  
  
  pdf(file=sprintf("%s\\figure_%s_pct_active%s.pdf", figs_path, fig_ext, ext),
      height=plot_h,
      width=plot_h * 1.5)
  plot(bar_plots_list[[1]])
  dev.off()   
  
  pdf(file=sprintf("%s\\figure_%s_pct_active_place_cells%s.pdf", figs_path, fig_ext, ext),
      height=plot_h,
      width=plot_h * 1.5)
  plot(bar_plots_list[[2]])
  dev.off()   
  
  pdf(file=sprintf("%s\\figure_%s_pct_place_cells_of_all%s.pdf", figs_path, fig_ext, ext),
      height=plot_h,
      width=plot_h * 1.5)
  plot(bar_plots_list[[3]])
  dev.off()   
  
  pdf(file=sprintf("%s\\figure_%s%s.pdf", figs_path, fig_ext, ext),
      height=plot_h * 3,
      width=plot_h * 3)
  plot(grid.arrange(aligned[[1]], aligned[[2]], aligned[[3]], nrow=3))
  dev.off()  
  
  
}


get_fit_params_figures <- function(path, likelihood_path="simulations_3", estimation_path="simulations_2_new", best_fit=T, new_simulation=F) {
  print(sprintf("Loading from %s - %s (%s)", likelihood_path, estimation_path, path))
  load(sprintf("%s\\%s\\pct_estimate_df.R", path, likelihood_path))
  if (best_fit) { 
    if (!new_simulation) {
      estimated_pct <- as.numeric(names(which.max(final_df[,"likelihood"]))) 
    } else {
      estimated_pct <- rev(seq(0.5,1,by=0.1))[which.max(final_df[,"likelihood"])]
    }
    
    #fit_indx <- order(final_df[,"likelihood"], decreasing = T)[3]
    #estimated_pct <- as.numeric(names(final_df[,"likelihood"])[fit_indx])
    print(estimated_pct)
  } else {
    if (!new_simulation) {
      estimated_pct <- as.numeric(names(which.min(final_df[,"likelihood"]))) 
    } else {
      estimated_pct <- rev(seq(0.5,1,by=0.1))[which.min(final_df[,"likelihood"])]
    }
  }
  
  #estimated_pct <- 1
  
  extr_df <- extract_simulation_params(path, est_path = estimation_path,pct_range = as.numeric(names(final_df[,"likelihood"])))
  print(sprintf("Got pct! %.2f", estimated_pct))
  return(list(params=extr_df[as.character(estimated_pct),],
             pct=estimated_pct))
  
}



simulation_data_frame <- do.call(rbind,
                                 list(c("simulations_2", "simulations_3", F, 5,  "C5M1", "R"),
                                      c("simulations_2", "simulations_3", F, 6,  "C5M1", "R"),
                                      c("simulations_2", "simulations_3", F, 7,  "C5M1", "R"),
                                      c("simulations_2", "simulations_3", F, 8,  "C5M1", "R"),
                                      c("simulations_2", "simulations_3", F, 13, "C5M1", "R"),
                                      c("simulations_2", "simulations_3", F, 14, "C5M1", "R"),
                                      c("simulations_2", "simulations_3", F, 15, "C5M1", "R"),
                                      c("simulations_2", "simulations_3", F, 16, "C5M1", "R"),
                                      c("simulations_2", "simulations_3", F, 5,  "C5M1", "L"),
                                      c("simulations_2", "simulations_3", F, 6,  "C5M1", "L"),
                                      c("simulations_2", "simulations_3", F, 7,  "C5M1", "L"),
                                      c("simulations_2", "simulations_3", F, 8,  "C5M1", "L"),
                                      c("simulations_2", "simulations_3", F, 13, "C5M1", "L"),
                                      c("simulations_2", "simulations_3", F, 14, "C5M1", "L"),
                                      c("simulations_2", "simulations_3", F, 15, "C5M1", "L"),
                                      c("simulations_2", "simulations_3", F, 16, "C5M1", "L"),
                                      c("simulations_new_func", "simulations_new_estimate", T, 5,  "C6M3", "R"),
                                      c("simulations_2", "simulations_3", F, 6,  "C6M3", "R"),
                                      c("simulations_2", "simulations_3", F, 7,  "C6M3", "R"),
                                      c("simulations_new_func", "simulations_new_estimate", T, 8,  "C6M3", "R"),
                                      c("simulations_new_func", "simulations_new_estimate", T, 13, "C6M3", "R"),
                                      c("simulations_2", "simulations_3", F, 14, "C6M3", "R"),
                                      c("simulations_2", "simulations_3", F, 15, "C6M3", "R"),
                                      c("simulations_new_func", "simulations_new_estimate", T, 16, "C6M3", "R"),
                                      c("simulations_2", "simulations_3", F, 5,  "C6M3", "L"),
                                      c("simulations_new_func_SI_1.5", "simulations_new_estimate_SI_1.5", T, 6,  "C6M3", "L"),
                                      c("simulations_2", "simulations_3", F, 7,  "C6M3", "L"),
                                      c("simulations_new_func_SI_1.5", "simulations_new_estimate_SI_1.5", T, 8,  "C6M3", "L"),
                                      c("simulations_new_func_SI_1.5", "simulations_new_estimate_SI_1.5", T, 13, "C6M3", "L"),
                                      c("simulations_2", "simulations_3", F, 14, "C6M3", "L"),
                                      c("simulations_2", "simulations_3", F, 15, "C6M3", "L"),
                                      c("simulations_2", "simulations_3", F, 16, "C6M3", "L"),
                                      c("simulations_2_new", "simulations_3", F, 5,  "C6M4", "R"),
                                      c("simulations_2_new", "simulations_3", F, 6,  "C6M4", "R"),
                                      c("simulations_2_new", "simulations_3", F, 7,  "C6M4", "R"),
                                      c("simulations_2_new", "simulations_3", F, 8,  "C6M4", "R"),
                                      c("simulations_2_new", "simulations_3", F, 13, "C6M4", "R"),
                                      c("simulations_2_new", "simulations_3", F, 14, "C6M4", "R"),
                                      c("simulations_2_new", "simulations_3", F, 15, "C6M4", "R"),
                                      c("simulations_2_new", "simulations_3", F, 16, "C6M4", "R"),
                                      c("simulations_2_new", "simulations_3", F, 5,  "C6M4", "L"),
                                      c("simulations_2_new", "simulations_3", F, 6,  "C6M4", "L"),
                                      c("simulations_2_new", "simulations_3", F, 7,  "C6M4", "L"),
                                      c("simulations_2_new", "simulations_3", F, 8,  "C6M4", "L"),
                                      c("simulations_2_new", "simulations_3", F, 13, "C6M4", "L"),
                                      c("simulations_2_new", "simulations_3", F, 14, "C6M4", "L"),
                                      c("simulations_2_new", "simulations_3", F, 15, "C6M4", "L"),
                                      c("simulations_2_new", "simulations_3", F, 16, "C6M4", "L"),
                                      c("simulations_2_new_new", "simulations_3_new", T, 5,  "C8M2", "R"),
                                      c("simulations_2_new_new", "simulations_3_new", T, 6,  "C8M2", "R"),
                                      c("simulations_2_new_new", "simulations_3_new", T, 7,  "C8M2", "R"),
                                      c("simulations_2_new_new", "simulations_3_new", T, 8,  "C8M2", "R"),
                                      c("simulations_2_new_new", "simulations_3_new", T, 13, "C8M2", "R"),
                                      c("simulations_2_new_new", "simulations_3_new", T, 14, "C8M2", "R"),
                                      c("simulations_2_new_new", "simulations_3_new", T, 15, "C8M2", "R"),
                                      c("simulations_2_new_new", "simulations_3_new", T, 16, "C8M2", "R"),
                                      c("simulations_2_new_new", "simulations_3_new", T, 5,  "C8M2", "L"),
                                      c("simulations_2_new_new", "simulations_3_new", T, 6,  "C8M2", "L"),
                                      c("simulations_2_new_new", "simulations_3_new", T, 7,  "C8M2", "L"),
                                      c("simulations_2_new_new", "simulations_3_new", T, 8,  "C8M2", "L"),
                                      c("simulations_2_new_new", "simulations_3_new", T, 13, "C8M2", "L"),
                                      c("simulations_2_new_new", "simulations_3_new", T, 14, "C8M2", "L"),
                                      c("simulations_2_new_new", "simulations_3_new", T, 15, "C8M2", "L"),
                                      c("simulations_2_new_new", "simulations_3_new", T, 16, "C8M2", "L")))


simulation_data_frame <- as.data.frame(simulation_data_frame)
simulation_data_frame[,3] <- as.logical(simulation_data_frame[,3])
simulation_data_frame[,4] <- as.numeric(simulation_data_frame[,4])
colnames(simulation_data_frame) <- c("Estimated", "Likelihood", "New", "Ses", "Mice", "Dir")

get_fit_params_figures_from_path <- function(path, best=T) {
  
  base <- str_split(str_split(path, "CA1\\\\")[[1]][2], "\\\\Matrices\\\\")
  base_2 <- str_split(base[[1]][[2]],  "\\\\equalized")
  mice <- base[[1]][1]
  direc <- base_2[[1]][1]
  session <- as.numeric(str_split(base_2[[1]][2], "session_")[[1]][2])
  direction <- ifelse(direc == "Left", "L", "R")
  
  row_ind <- which(simulation_data_frame[,"Mice"] == mice  &
                     simulation_data_frame[,"Dir"] == direction & 
                     simulation_data_frame[,"Ses"] == session)
  
  params_path <- simulation_data_frame[row_ind,]
  
  print(sprintf("Loading from %s - %s (%s)", params_path["Likelihood"], params_path["Estimated"], path))
  
  print(sprintf("%s\\%s\\pct_estimate_df.R", path, params_path["Likelihood"]))
  load(sprintf("%s\\%s\\pct_estimate_df.R", path, params_path["Likelihood"]))
  
  if (best) { 
    if (!params_path["New"]) {
      estimated_pct <- as.numeric(names(which.max(final_df[,"likelihood"]))) 
    } else {
      estimated_pct <- rev(seq(0.5,1,by=0.1))[which.max(final_df[,"likelihood"])]
    }
    
    #fit_indx <- order(final_df[,"likelihood"], decreasing = T)[3]
    #estimated_pct <- as.numeric(names(final_df[,"likelihood"])[fit_indx])
    print(estimated_pct)
  } else {
    if (!params_path["New"]) {
      estimated_pct <- as.numeric(names(which.min(final_df[,"likelihood"]))) 
    } else {
      estimated_pct <- rev(seq(0.5,1,by=0.1))[which.min(final_df[,"likelihood"])]
    }
  }
  
  #estimated_pct <- 1
  
  extr_df <- extract_simulation_params(path, est_path = params_path["Estimated"],pct_range = as.numeric(names(final_df[,"likelihood"])))
  print(sprintf("Got pct! %.2f", estimated_pct))
  return(list(params=extr_df[as.character(estimated_pct),],
              pct=estimated_pct))
  
}


get_cost_compare_vectors_by_path <- function(path, true_spike_train, stim_trace, pval_func=wilcox.test, ret=F, smooth_lambda=-1, best=T, new_simulation=F, np="simulations_2") {
  
  firing_rate <-rowMeans(true_spike_train) / dt
  processed_real <- preprocess_spike_train(true_spike_train, stim_trace=stim_trace, verbose=F)
  true_cells_spike_train <- processed_real$working_cells_spike_train
  true_firing_rate <- processed_real$working_firing_rate
  true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
  
  tmp <- compute_tuning(true_cells_spike_train, stim_trace)
  stim_prob <- tmp[[1]]
  true_tuning_curve <- tmp[[2]]
  rm(tmp)
  
  
  if (smooth_lambda == - 1){
    tuning_for_peaks <- true_tuning_curve
  } else {
    tuning_for_peaks <- 
      t(apply(true_tuning_curve, 1, function(trc) {
        barp <- barplot(trc, plot=F)
        smt <- smooth.spline(barp, trc, all.knots = T, lambda=smooth_lambda)
        smt$y[smt$y < 0] <- 0
        return(smt$y);
      }))
  }
  
  true_peaks <- unlist(apply(tuning_for_peaks, 1, 
                             function(n){return(len(get_peaks(n, threshold_from_max = peaks_significance_threshold)))}))
  true_active_spatial_bins_per_cell <- apply(true_tuning_curve, 1, function(n){sum(n>0)})
  #print(true_active_spatial_bins_per_cell)
  true_SI <- compute_SI(stim_prob, true_tuning_curve, true_firing_rate)
  
  params <- get_fit_params_figures(path,
                                   estimation_path=ifelse(new_simulation, "simulations_new_func", "simulations_2_new"),
                                   likelihood_path=ifelse(new_simulation, "simulations_new_estimate", "simulations_3"),
                                   best_fit = best,
                                   new_simulation=new_simulation)


#  params <- get_fit_params_figures_from_path(path, best = best)
  SI_f <- c()
  tb_f <- c()
  peaks_f <- c()
  sb_f <- c()
  simulated_tuning_curve <-
    generate_tuning_curves_cost(n = nrow(spike_train),
                                percentage = params$pct,
                                average_width = params$params["average_width"], 
                                sd_width = params$params["sd_width"],
                                fixed_fr=firing_rate,
                                noise=params$params["noise"],
                                double_peak_pct = params$params["double_peak_pct"],
                                plot=F)
  pois_factor <- currate_spike_train_cost(simulated_tuning_curve, 
                                          true_time_bins_per_cells,
                                          stim_trace,
                                          verbose = F)
  
  
  for (i in 1:10) { 
    simulated_tuning_curve <-
      generate_tuning_curves_cost(n = nrow(spike_train),
                                  percentage = params$pct,
                                  average_width = params$params["average_width"], 
                                  sd_width = params$params["sd_width"],
                                  fixed_fr=firing_rate,
                                  noise=params$params["noise"],
                                  double_peak_pct = params$params["double_peak_pct"],
                                  plot=F)
  

  generated_spike_train <- 
    generate_spike_trains_cost(tuning_curves = simulated_tuning_curve,
                               stim_trace = stim_trace,
                               factor=pois_factor,
                               fs=1)

  
  processed_generated <- preprocess_spike_train(generated_spike_train, stim_trace, verbose=F)
  generated_cells_active_st <- processed_generated$working_cells_spike_train
  generated_firing_rate <- processed_generated$working_firing_rate
  generated_time_bins <- processed_generated$working_time_bins_per_cells
  
  tmp <- compute_tuning(generated_cells_active_st, stim_trace)
  gen_stim_prob <- tmp[[1]]
  gen_tuning_curve <- tmp[[2]]
  rm(tmp)
  
  if (smooth_lambda == - 1){
    tuning_for_peaks <- gen_tuning_curve
  } else {
    tuning_for_peaks <- 
      t(apply(gen_tuning_curve, 1, function(trc) {
        barp <- barplot(trc, plot=F)
        smt <- smooth.spline(barp, trc, all.knots = T, lambda=smooth_lambda)
        smt$y[smt$y < 0] <- 0
        return(smt$y);
      }))
  }
  
  generated_peaks <- unlist(apply(tuning_for_peaks, 1, 
                                  function(n){return(len(get_peaks(n, threshold_from_max = peaks_significance_threshold)))}))
  generated_active_spatial_bins_per_cell <- apply(gen_tuning_curve, 1, function(n){sum(n>0)})
  generated_SI <- compute_SI(gen_stim_prob, gen_tuning_curve, generated_firing_rate)
  
  
  if (ret) {
  return(list(true_SI=true_SI[[1]],
              true_tb=true_time_bins_per_cells,
              true_peaks=true_peaks,
              true_sb=true_active_spatial_bins_per_cell,
              gen_SI=generated_SI[[1]],
              gen_tb=generated_time_bins,
              gen_peaks=generated_peaks,
              gen_sb=generated_active_spatial_bins_per_cell))
  }

  SI=pval_func(true_SI[[1]], generated_SI[[1]])$p.value
  tb=pval_func(true_time_bins_per_cells, generated_time_bins)$p.value
  peaks=pval_func(true_peaks, generated_peaks)$p.value
  sb=pval_func(true_active_spatial_bins_per_cell, generated_active_spatial_bins_per_cell)$p.value
  
  print(sprintf("SI %.3f, tb %.3f, peaks %.3f, sb %.3f",
        SI,
        tb,
        peaks,
        sb))
  SI_f <- c(SI_f, SI)
  tb_f <- c(tb_f, tb)
  peaks_f <- c(peaks_f, peaks)
  sb_f <- c(sb_f, sb)
  
  }
  
  return(c(mean(SI_f), mean(tb_f), mean(peaks_f), mean(sb_f)))
  
}




figure_2_all_histograms <- function(paths, out_path, sessions_to_use=c(5:8,13:16), verbose=F) {
  stimulus_trace_filename = "\\stim_trace.mat"
  spike_train_filename = "\\spike_train.mat"
  
  equalized_paths <- sprintf("%s\\%s", paths, "equalized")
  
  most_likely <- list()
  least_likely <- list()
  
  for (path_idx in 1:len(paths)) {
    
    #if (path_idx %in%  5:8) {
      np = "simulations_2_new"
    #} else {
      #np = "simulations_2"
    #}
    
    working_path <- paths[path_idx]
    working_path <- paths[path_idx]
    equalized_path <- equalized_paths[path_idx]
    
    m <- readMat(paste(working_path, stimulus_trace_filename, sep=""))
    
    if(len(grep("Left", working_path)) > 0) {
      stim_trace <- m$position.left.trials
    } else {
      stim_trace <- m$position.right.trials
    }
    
    
    m <- readMat(paste(working_path, spike_train_filename, sep=""))
    
    if(len(grep("Left", working_path)) > 0) {
      spike_train <- m$spikes.left.trials
    } else {
      spike_train <- m$spikes.right.trials
    }
    
    
    session_paths <- 
      list.dirs(equalized_path, recursive=F)[grep("session", 
                                                list.dirs(equalized_path, 
                                                          recursive=F))]
    
    
    tmp_spike_train <- spike_train
    tmp_stimulus_trace <- stim_trace
    
    # Run through all paths
    for (tpath in session_paths) {
      
      
      ses_idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
      
      # Ignore session
      if (sum(ses_idx == sessions_to_use) == 0){
        next
      }
      
      spike_train <- tmp_spike_train
      stim_trace <- tmp_stimulus_trace
      spike_train <- spike_train[[ses_idx]][[1]]
      spike_train <- t(as.matrix(spike_train))
      stim_trace <- stim_trace[[ses_idx]][[1]]
      stim_trace <-  as.vector(stim_trace)
      
      ###### equalize spike train
      equalize_bins <- c(1:3,22:24)
      indices <- 1:len(stim_trace)
      
      
      final_ind <- rep(T, times=len(stim_trace))
      
      for (ebin in equalize_bins) {
        print(ebin)
        
        num_to_draw <- gen_num_frames(stim_trace)
        if (num_to_draw > len(which(stim_trace == ebin))) {
          print(sprintf("No need to equalize bin %d", ebin))
          next
        }
        
        sampled_ind <- sample(which(stim_trace==ebin), num_to_draw)
        tmp <- indices %in% sampled_ind | stim_trace[indices] != ebin       
        final_ind <- final_ind & tmp
      }
      
      spike_train <- spike_train[,final_ind]
      stim_trace <- stim_trace[final_ind]

      most <- get_cost_compare_vectors_by_path(tpath, spike_train, stim_trace, ret=F, smooth_lambda = 1e-5, pval_func = wilcox.test, np=np)
      #least <- get_cost_compare_vectors_by_path(tpath, spike_train, stim_trace, ret=F, smooth_lambda = 1e-5, pval_func = wilcox.test, np=np, best = F)
      print(tpath)
      #least <- get_cost_compare_vectors_by_path(tpath, spike_train, stim_trace, ret=F, smooth_lambda = 1e-5, pval_func = wilcox.test, new_simulation = T)

      most_likely <- append(most_likely, list(c(most, ses_idx)))
      least_likely <- append(least_likely, list(c(most, ses_idx)))
    }
  }
  
  df_most <- do.call(rbind, most_likely)
  df_least <- do.call(rbind, least_likely)
  
  melted_dfs <- list()
  for (fdf in list(df_most[,1:4], df_least[,1:4])) {
    colnames(fdf) <- c("SI", "Time bins", "Peaks", "Spatial bins")
    fdf_2 <- fdf
    #fdf_2 <- fdf[,c("SI", "Time bins", "Spatial bins")]
    melted_df <- melt(fdf_2)
    colnames(melted_df) <- c("#", "Comparsion", "Pval")
    
    melted_df$Ses <- melted_df$`#` %% 8 + 8
    melted_df$Ses[melted_df$Ses == 8] <- melted_df$Ses[melted_df$Ses == 8] + 8
    melted_df$Ses <- as.character(melted_df$Ses)
    
    melted_dfs <- append(melted_dfs, list(melted_df))
    
  }
    
    fdf <- rbind(melted_dfs[[1]], melted_dfs[[2]])
    fdf$Fit = c(rep("Most likely", times=nrow(melted_dfs[[1]])),rep("Least likely", times=nrow(melted_dfs[[2]])))
    fdf$Comparsion <- factor(fdf$Comparsion, levels=c("Peaks", "SI", "Spatial bins", "Time bins"))
    g <-
    ggplot(melted_dfs[[1]], aes(x=Comparsion, y=Pval)) +
      geom_violin(alpha=0.1, size=0.75, position=position_dodge(), aes(x=Comparsion, y=Pval),
                  color="darkmagenta", fill="darkmagenta") +
      geom_jitter(position=position_jitterdodge(0.2), aes(col=`Ses`), size=0.5) +
      geom_hline(yintercept=0.05, linetype="dashed", size=0.7) +
      stat_summary(fun=mean, geom="point", size=2.5, color="red", aes(x=Comparsion)) + 
      theme_light() +
      #scale_color_manual(values=c("gray80", "gray60", "gray40", "gray20"), "Session") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            #legend.position="NA",
            axis.text.x = element_text(angle = 45),
            panel.border = element_blank(),
            panel.background = element_blank())
    
    ggplot(fdf, aes(x=Comparsion, y=Pval, color=Fit, fill=Fit)) +
      geom_jitter(data=fdf, position=position_jitterdodge(0.2), aes(x=Comparsion, y=Pval,col=`Ses`, group=Fit)) +
      geom_hline(yintercept=0.05, linetype="dashed", size=1.25) + 
      geom_boxplot(size=1, alpha=0.1) + 
      #stat_summary(fun=mean, geom="point", size=5, color="red", aes(x=Comparsion, group=Fit)) + 
      scale_color_manual(breaks=c("5","6","7","8","Most likely", "Least likely"), 
                         values=c("gray80", "gray60", "gray40", "gray20", "darkmagenta", "goldenrod")) +
      scale_fill_manual(breaks=c("Most likely", "Least likely"), 
                         values=c("darkmagenta", "goldenrod")) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border = element_blank(),
              panel.background = element_blank())

     gf <- grid.arrange(gb, g, g2, nrow=1)
}




get_spike_train_and_stim_trace_from_path <- function(path, session, old=F, equalize_frames=T, simulate=F) {

  
  if (old) {
    m <- readMat(paste(path, "\\stimulus_trace.mat", sep=""))
    stim_trace <- m$stimulus.trace
    
    
    m <- readMat(paste(path, spike_train_filename, sep=""))
    spike_train <- m$spike.train
  } else {  
    m <- readMat(paste(path, stimulus_trace_filename, sep=""))
    
    if(len(grep("Left", path)) > 0) {
      stim_trace <- m$position.left.trials
    } else {
      stim_trace <- m$position.right.trials
    }
    
    
    m <- readMat(paste(path, spike_train_filename, sep=""))
    
    if(len(grep("Left", path)) > 0) {
      spike_train <- m$spikes.left.trials
    } else {
      spike_train <- m$spikes.right.trials
    }
  }
  
  tmp_path <- path
  tmp_spike_train <- spike_train
  tmp_stimulus_trace <- stim_trace
  
  
  spike_train <- tmp_spike_train
  stim_trace <- tmp_stimulus_trace
  spike_train <- spike_train[[session]][[1]]
  spike_train <- t(as.matrix(spike_train))
  stim_trace <- stim_trace[[session]][[1]]
  stim_trace <-  as.vector(stim_trace)
    
    

    if (!old && equalize_frames ) {
      equalize_bins <- c(1:3,22:24)
      indices <- 1:len(stim_trace)
      
      
      final_ind <- rep(T, times=len(stim_trace))
      
      for (ebin in equalize_bins) {
        print(ebin)
        
        num_to_draw <- gen_num_frames(stim_trace)
        if (num_to_draw > len(which(stim_trace == ebin))) {

          next
        }
        
        sampled_ind <- sample(which(stim_trace==ebin), num_to_draw)
        tmp <- indices %in% sampled_ind | stim_trace[indices] != ebin       
        final_ind <- final_ind & tmp
      }
      
      spike_train <- spike_train[,final_ind]
      stim_trace <- stim_trace[final_ind]
    }
    
    if (simulate) {
      fr <- rowMeans(spike_train) / dt
      processed_real <- preprocess_spike_train(spike_train, stim_trace)
      
      true_cells_spike_train <- processed_real$working_cells_spike_train
      true_firing_rate <- processed_real$working_firing_rate
      true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
      
      params <- get_fit_params(path, estimation_path="simulations_2_new")
      simulated_tuning_curve <-
        generate_tuning_curves_cost(n = nrow(spike_train),
                                    percentage = params$pct,
                                    average_width = params$params["average_width"], 
                                    sd_width = params$params["sd_width"],
                                    fixed_fr=fr,
                                    noise=params$params["noise"],
                                    double_peak_pct = params$params["double_peak_pct"],
                                    plot=F)
      
      pois_factor <- currate_spike_train_cost(simulated_tuning_curve, 
                                              true_time_bins_per_cells,
                                              stim_trace,
                                              verbose = T)
      
      generated_spike_train <- 
        generate_spike_trains_cost(tuning_curves = simulated_tuning_curve,
                                   stim_trace = stim_trace,
                                   factor=pois_factor,
                                   fs=1)
      
      orig_st <- spike_train
      orig_path <- path
      spike_train <- generated_spike_train

    }

  
  return (list(spike_train, stim_trace))
}

