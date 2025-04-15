library(RColorBrewer)
library(dunn.test)

supp_figure_npcs_tuning_corr_heatmaps <- function() {
  write_path <- sprintf("%s\\supp_figure_npcs\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\supp_figure_npcs\\tuning_coor_heatmaps\\",figures_path)
  dir.create(write_path)
  
  load(sprintf("%s\\samples\\sliding_tuning\\sliding_df_all_by_time.Rda", 
               base_path), 
       verbose=T)

  dir.create(write_path)
  dir.create(sprintf("%s\\big\\",write_path))
  dir.create(sprintf("%s\\big_1\\",write_path))
  dir.create(sprintf("%s\\medium_1\\",write_path))
  dir.create(sprintf("%s\\medium_2\\",write_path))
  dir.create(sprintf("%s\\medium\\",write_path))
  dir.create(sprintf("%s\\small\\",write_path))
  
  sizes <- c(big=3,
             big_1=2.75,
             medium=2.5,
             medium_1=2,
             medium_2=1.75,
             small=1.5)
  
  
  mt <- 
  do.call(rbind,
          lapply(str_split(windows, "_"), 
                 function(wind_jump_str) {as.numeric(unlist(str_split(unlist(str_split(wind_jump_str, "W")), "J"))[c(2,4)])}))  
  
  windows <- unique(mt[,1])
  jumps <- unique(mt[,2])
  
  cor_range <- seq(-.05,.25, length.out=50)
  cor_colors <- colorRampPalette(brewer.pal(11, "RdYlBu"))(len(cor_range))
  
  sliding_tuning_df <- res$df
  median_sliding_vals <- res$median
  mean_sliding_vals <- res$mean
  
  
  session_groups <- list(all=1:16)
  
  
  df_res <- data.frame()
  
  metrics <- unique(sliding_tuning_df$Metric)
  subfields <- unique(sliding_tuning_df$Subfield)
  #windows <- unique(sliding_tuning_df$WindowP)
  
  ylabs_to_use <- c(peaks="Peak movement (Bins)",
                    cor="Tuning correlation (sliding)",
                    cosine="Tuning cosine distance (sliding)",
                    eucl="Eucledian distance (sliding)",
                    wasserstein="Wasserstein distance (sliding)")
  
  
  for (session_group_name in names(session_groups)) {
    sessions_to_use <- session_groups[[session_group_name]]
    
    
    for (metric in c("cor")) {
      for (sbf in unique(sliding_tuning_df$Subfield)) {
        
        mean_plots <- list()
        median_plots <- list()
        sd_plots_mean <- list()
        sd_plots_median <- list()
        mean_sd_plots_mean <- list()
        mean_sd_plots_median <- list()
        
        for (grp in c("Pcs", "Npcs")) {
          for (sim_group in unique(sliding_tuning_df$isSim)) {
            
            mean_metric_mat_mean <-
              matrix(rep(NA, times=len(windows) * len(jumps)),
                     nrow = len(windows))

            sd_metric_mat_mean <-
              matrix(rep(NA, times=len(windows) * len(jumps)),
                     nrow = len(windows))
            
            mean_metric_mat_median <-
              matrix(rep(NA, times=len(windows) * len(jumps)),
                     nrow = len(windows))
            
            sd_metric_mat_median <-
              matrix(rep(NA, times=len(windows) * len(jumps)),
                     nrow = len(windows)) 
            
            
            mean_sd_metric_mat_mean <- 
              matrix(rep(NA, times=len(windows) * len(jumps)),
                     nrow = len(windows)) 
            mean_sd_metric_mat_median <- 
              matrix(rep(NA, times=len(windows) * len(jumps)),
                     nrow = len(windows)) 
            
            rownames(mean_metric_mat_mean) <- windows
            rownames(sd_metric_mat_mean) <- windows
            rownames(mean_sd_metric_mat_mean) <- windows
            rownames(mean_metric_mat_median) <- windows
            rownames(sd_metric_mat_median) <- windows
            rownames(mean_sd_metric_mat_median) <- windows
            
            colnames(mean_metric_mat_mean) <- jumps
            colnames(sd_metric_mat_mean) <- jumps
            colnames(mean_sd_metric_mat_mean) <- jumps
            colnames(mean_metric_mat_median) <- jumps
            colnames(sd_metric_mat_median) <- jumps
            colnames(mean_sd_metric_mat_median) <- jumps
            
            for (window in windows) {
              for (jmp in jumps) {
                wind_to_use <- sprintf("W%d_J%d", window, jmp)
                
                ind <- sliding_tuning_df$Session %in% sessions_to_use &
                  sliding_tuning_df$Metric == metric &
                  sliding_tuning_df$WindowP == wind_to_use &
                  sliding_tuning_df$Group == grp &
                  sliding_tuning_df$Subfield == sbf & 
                  sliding_tuning_df$isSim == sim_group
                
                df_to_use <- sliding_tuning_df[ind, ]
                median_sliding_vals_to_use <- median_sliding_vals[ind]
                mean_sliding_vals_to_use <- mean_sliding_vals[ind]
                
                mean_res_df <- data.frame()
                median_res_df <- data.frame()
                  
                for (trace_idx in 1:sum(ind)) {
                  median_values <- rep(NA, times=40)
                  mean_values <- rep(NA, times=40)
                    
                  
                  # As sessions are in different length, 
                  # Cut them into 40 pieces according to their relative size
                  median_values_ind <- as.numeric(cut(seq(1.025,
                                                          20,
                                                          length.out=len(median_sliding_vals_to_use[[trace_idx]])),
                                                       breaks=seq(1,20,length.out=40)))
                    
                  mean_values_ind <- as.numeric(cut(seq(1.025,
                                                        20,
                                                        length.out=len(mean_sliding_vals_to_use[[trace_idx]])),
                                                    breaks=seq(1,20,length.out=40)))
                    
                  if (len(mean_values_ind) > 1) {
                      mean_values[mean_values_ind] <- mean_sliding_vals_to_use[[trace_idx]]
                  }
                    
                  if (len(median_values_ind) > 1) {
                      median_values[median_values_ind] <- median_sliding_vals_to_use[[trace_idx]]
                  }
                    
                  mean_res_df <- rbind(mean_res_df,mean_values)
                  median_res_df <- rbind(median_res_df,median_values)              
                }
                
                # As each session is in different lenght,
                # Get the mean trace of that session 
                #
                
                
                sd_res_mean <- sd(apply(mean_res_df, 1, mean, na.rm=T), na.rm=T)
                mean_res_mean <- mean(apply(mean_res_df, 1, mean , na.rm=T), na.rm=T)
                sd_res_median <- sd(apply(median_res_df, 1, mean, na.rm=T), na.rm=T)
                mean_res_median <- mean(apply(median_res_df, 1, mean, na.rm=T), na.rm=T)
                
                mean_metric_mat_mean[as.character(window),
                                     as.character(jmp)] <- mean_res_mean
                
                mean_sd_metric_mat_mean[as.character(window),
                                     as.character(jmp)] <- mean_res_mean + sd_res_mean
                
                sd_metric_mat_mean[as.character(window),
                                   as.character(jmp)] <- sd_res_mean                
                
                mean_metric_mat_median[as.character(window),
                                     as.character(jmp)] <- mean_res_median 
                
                mean_sd_metric_mat_median[as.character(window),
                                        as.character(jmp)] <- mean_res_median + sd_res_mean
                
                sd_metric_mat_median[as.character(window),
                                   as.character(jmp)] <- sd_res_median                
                
                
              }
            }
            
            mean_mean_ph <- 
            pheatmap(mean_metric_mat_mean, 
                     cluster_rows=F, 
                     cluster_cols=F, 
                     breaks=cor_range, 
                     col=rev(cor_colors), 
                     border_col=NA,
                     main=sprintf("%s-%s", grp, sim_group))
            
            mean_median_ph <- 
              pheatmap(mean_metric_mat_median, 
                       cluster_rows=F, 
                       cluster_cols=F, 
                       breaks=cor_range, 
                       col=rev(cor_colors), 
                       border_col=NA,
                       main=sprintf("%s-%s", grp, sim_group))
            
            
            mean_sd_mean_ph <- 
              pheatmap(mean_sd_metric_mat_mean, 
                       cluster_rows=F, 
                       cluster_cols=F, 
                       breaks=cor_range, 
                       col=rev(cor_colors), 
                       border_col=NA,
                       main=sprintf("%s-%s", grp, sim_group))
            
            mean_sd_median_ph <- 
              pheatmap(mean_sd_metric_mat_median, 
                       cluster_rows=F, 
                       cluster_cols=F, 
                       breaks=cor_range, 
                       col=rev(cor_colors), 
                       border_col=NA,
                       main=sprintf("%s-%s", grp, sim_group))            
            
            sd_mean_ph <- 
              pheatmap(sd_metric_mat_mean, 
                       cluster_rows=F, 
                       cluster_cols=F, 
                       breaks=cor_range, 
                       col=rev(cor_colors), 
                       border_col=NA,
                       main=sprintf("%s-%s", grp, sim_group))
            
            sd_median_ph <- 
              pheatmap(sd_metric_mat_median, 
                       cluster_rows=F, 
                       cluster_cols=F, 
                       breaks=cor_range, 
                       col=rev(cor_colors), 
                       border_col=NA,
                       main=sprintf("%s-%s", grp, sim_group))
            
            mean_plots <- append(mean_plots,
                                 list(mean_mean_ph[[4]]))
            median_plots <- append(median_plots,
                                 list(mean_median_ph[[4]]))       
            
            sd_plots_mean <- append(sd_plots_mean,
                                    list(sd_mean_ph[[4]]))       
            sd_plots_median <- append(sd_plots_median,
                                      list(sd_median_ph[[4]]))      
            
            mean_sd_plots_mean <- append(mean_sd_plots_mean,
                                         list(mean_sd_mean_ph[[4]]))      
            mean_sd_plots_median <- append(mean_sd_plots_median,
                                           list(mean_sd_median_ph[[4]]))       
          }
        }
        
        mean_plots$nrow <- 1
        median_plots$nrow <- 1
        sd_plots_mean$nrow <- 1
        sd_plots_median$nrow <- 1
        mean_sd_plots_mean$nrow <- 1
        mean_sd_plots_median$nrow <- 1
        
        median_f <- do.call(arrangeGrob, median_plots)
        mean_f <- do.call(arrangeGrob, mean_plots)
        sd_mean_f <- do.call(arrangeGrob, sd_plots_mean)
        sd_median_f <- do.call(arrangeGrob, sd_plots_median)
        mean_sd_mean_f <- do.call(arrangeGrob, mean_sd_plots_mean)
        mean_sd_median_f <- do.call(arrangeGrob, mean_sd_plots_median)

        for (size_name in names(sizes)) {
          size = sizes[[size_name]]
          write_path_session <- sprintf("%s\\%s\\",write_path, size_name)
          dir.create(write_path_session)
          
          pdf(sprintf("%s//mean_corr_heatmap_%s.pdf",
                      write_path_session,
                      sbf),
              height=size * 1,
              width=size * 4)
          plot(mean_f) 
          dev.off()
          
          pdf(sprintf("%s//median_corr_heatmap_%s.pdf",
                      write_path_session,
                      sbf),
              height=size * 1,
              width=size * 4)
          plot(median_f) 
          dev.off()
          
          pdf(sprintf("%s//mean_sd_mean_corr_heatmap_%s.pdf",
                      write_path_session,
                      sbf),
              height=size * 1,
              width=size * 4)
          plot(mean_sd_mean_f) 
          dev.off()
          
          pdf(sprintf("%s//mean_sd_median_corr_heatmap_%s.pdf",
                      write_path_session,
                      sbf),
              height=size * 1,
              width=size * 4)
          plot(mean_sd_median_f) 
          dev.off()

          pdf(sprintf("%s//sd_mean_corr_heatmap_%s.pdf",
                      write_path_session,
                      sbf),
              height=size * 1,
              width=size * 4)
          plot(sd_mean_f) 
          dev.off()
          
          pdf(sprintf("%s//sd_median_corr_heatmap_%s.pdf",
                      write_path_session,
                      sbf),
              height=size * 1,
              width=size * 4)
          plot(sd_median_f) 
          dev.off()          
        }
      }
    }
  }
}

supp_figure_npcs_properties_comparisions <- function() {
  
  write_path <- sprintf("%s\\supp_figure_npcs\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\supp_figure_npcs\\properties_comp_plots\\",figures_path)
  dir.create(write_path)
  
  groups_to_use = c(1:3)
  true_data_paths <- sprintf("%s\\%s", all_data_paths, "equalized")
  file_name="properties.Rda"
  sessions_to_use=c(1:16)
  
  dir.create(sprintf("%s\\big\\", write_path))
  dir.create(sprintf("%s\\medium_1\\", write_path))
  dir.create(sprintf("%s\\medium_2\\", write_path))
  dir.create(sprintf("%s\\medium\\", write_path))
  dir.create(sprintf("%s\\small\\", write_path))
  
  sizes = c(big=2,
            medium=1.75,
            medium_1=1.5,
            medium_2=1.25,
            small=1.15)
  
  comp_groups <-
    c("PC_REAL_PC_SIM",
      "NPC_REAL_NPC_SIM",
      "PC_REAL_NPC_REAL",
      "PC_SIM_NPC_SIM")
  
  
  all_properties_list <- list()
  indices_df <- data.frame()
  

  
  for (working_path in true_data_paths) {
    
    # extract all paths
    session_paths <- 
      list.dirs(working_path, recursive=F)[grep("session", list.dirs(working_path, recursive=F))]
    
    for (ext in c("", "\\fit_simulation")) { 
      
      # Run through all paths
      for (tpath in session_paths) {
        
        
        idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
        if (sum(idx == sessions_to_use) == 0){
          next
        }
        
        if (len((grep(file_name, list.files(sprintf("%s%s", tpath, ext))))) == 0) {
          print(tpath)
          missing_paths <- c(missing_paths, tpath)
          print("Missing")
          next
        }
        
        if (verbose) {
          print(sprintf("Loading %s", sprintf("%s%s\\%s", tpath, ext, file_name)))
        }
        load(sprintf("%s%s\\%s", tpath, ext, file_name), verbose=T)
        
        
        if (grepl("CA3", working_path)) {
          mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', tpath))
          mice_str <- substr(tpath, mice_str_index, mice_str_index+4) 
        } else {
          mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', tpath))
          mice_str <- substr(tpath, mice_str_index, mice_str_index+3) 
        }
        
        subfield <- ifelse(grepl("CA3", working_path), "CA3", "CA1")
        sim <- ifelse(grepl("fit_simulation", ext), "Simulation", "Real")
        
        for (group in names(result$tuning_properties_by_group)) {
          properties <- result$tuning_properties_by_group[[group]]
          properties <- properties[,-ncol(properties)]
          properties_list <- list(properties)
          
          indices_df <- rbind(indices_df,
                              data.frame(
                                nrow(indices_df) + 1, 
                                idx, 
                                ifelse(grepl("Left", tpath), "Left", "Right"), 
                                mice_str, 
                                subfield,
                                group,
                                sim))
          
          all_properties_list <- append(all_properties_list,
                                        list(properties_list))
          
        }
      }
    }
  }
  
  metric_prop_comp_dfs_list <- list()
  
  properties_to_use <- c("SI",
                         "peaks",
                         "spatial_bins",
                         "time_bins")
  metric_functions <- list(KL=KL_dist,
                           JSD=JSD_dist,
                           dprime=dprime_dist,
                           wasserstein=wasserstein_dist)
  
  for (metric_name in names(metric_functions)) {
    metric_prop_comp_dfs_list[[metric_name]] <- data.frame()
  }
  
  comp_list <- list(`PC_REAL_PC_SIM`=
                     list(group_1=list(isSim="Real", Group="cyclic_pcs"),
                          group_2=list(isSim="Simulation", Group="cyclic_pcs")),
                    `NPC_REAL_NPC_SIM`=
                      list(group_1=list(isSim="Real", Group="cyclic_non_pcs"),
                           group_2=list(isSim="Simulation", Group="cyclic_non_pcs")),
                    `PC_REAL_NPC_REAL`=
                      list(group_1=list(isSim="Real", Group="cyclic_pcs"),
                           group_2=list(isSim="Real", Group="cyclic_non_pcs")),
                    `PC_SIM_NPC_SIM`=
                      list(group_1=list(isSim="Simulation", Group="cyclic_pcs"),
                           group_2=list(isSim="Simulation", Group="cyclic_non_pcs")))
  colnames(indices_df) <- c("Index",
                            "Session",
                            "Direction",
                            "Mice",
                            "Subfield",
                            "Group",
                            "isSim")
  for (Subfield in unique(indices_df$Subfield)) {
    for (Mish in unique(indices_df$Mice)) {
     for (Session in unique(indices_df$Session)) {
       for (Direc in unique(indices_df$Direction)) {
          work_df <- indices_df[indices_df$Subfield == Subfield &
                                indices_df$Mice == Mish &
                                indices_df$Session == Session &
                                indices_df$Direction == Direc,]
          
          if (nrow(work_df) == 0) {
            next
          }
          
          for (comp_group_name in names(comp_list)){
            comp_group <- comp_list[[comp_group_name]]
            
            grp_1_df <- work_df[work_df$Group == comp_group$group_1$Group &
                                work_df$isSim == comp_group$group_1$isSim,]

            grp_2_df <- work_df[work_df$Group == comp_group$group_2$Group &
                                  work_df$isSim == comp_group$group_2$isSim,]    
            
            prop_1 <- all_properties_list[[grp_1_df$Index]][[1]]
            prop_2 <- all_properties_list[[grp_2_df$Index]][[1]]
            
            for (prop in properties_to_use) {
              for (metric_name in names(metric_functions)) {
                metric <- metric_functions[[metric_name]]
                

                
                metric_prop_comp_dfs_list[[metric_name]] <- 
                  rbind(metric_prop_comp_dfs_list[[metric_name]],
                        c(Dist=metric(prop_1[,prop], prop_2[,prop]),
                          Subfiled=Subfield,
                          Property=prop,
                          Session=Session,
                          group=comp_group_name))
              }
            }
          }
        }
      } 
    }
  }
  
  metric_prop_plt_lists <- list()
  
  for (metric_name in names(metric_functions)) {
    metric_prop_plt_lists[[metric_name]] <- list()
    colnames(metric_prop_comp_dfs_list[[metric_name]]) <- 
      c("Dist", "Subfield", "Prop", "Session", "Group")
    metric_prop_comp_dfs_list[[metric_name]] <- as.data.frame(metric_prop_comp_dfs_list[[metric_name]])
    metric_prop_comp_dfs_list[[metric_name]]$Dist <- as.numeric(metric_prop_comp_dfs_list[[metric_name]]$Dist)
  }
  pooled_sessions_statistics_df <- data.frame()
  
  for (metric_name in names(metric_functions)) {

    
    subfields_plots <- list()
    compare_boxplots <- list()
    
    for (subfield in c("CA1", "CA3")) {
      
      property_plots <- list()      
      
      for (prop in properties_to_use) {
        
        work_df <- metric_prop_comp_dfs_list[[metric_name]]
        work_df <- work_df[work_df$Prop == prop & 
                             work_df$Group %in% comp_groups[1:2],]
        
        work_df$Group <- factor(work_df$Group, levels=comp_groups[groups_to_use])
        work_df$Dist <- (work_df$Dist - min(work_df$Dist)) / 
          (max(work_df$Dist) - min(work_df$Dist))
        
        
        work_df_subfield <- work_df[work_df$Subfield == subfield,]
        g <- 
          ggplot(work_df_subfield, aes(x=Group, y=Dist)) +
          #geom_jitter(aes(color=Group), alpha=0.3, size=1,
          #            position=position_jitterdodge(.5)) + 
          geom_violin(aes(fill=Group), size=.5, alpha=0.8) + 
          stat_summary(fun=mean, geom="point", size=1.5, color="black", aes(x=Group, y=Dist)) + 
          base_plot_theme +
          theme(axis.text.x  = element_blank(),#element_text(angle = 45, vjust = 1, hjust=1),
                axis.title.y=element_text(size=9.5),
                axis.text.y=element_text(size=12.5, color="black"),
                plot.margin=unit(c(5.5,0,5.5, ifelse(prop != "SI", 0, 5.5)), "pt")) + 
          ylab(ifelse(prop != "SI", "", sprintf("Distance (%s)", metric_name))) + 
          xlab("") +
          ggtitle(sprintf("%s - %s", subfield, prop)) +

          geom_line(data=data.frame(x=c(1.15,1.15,
                                        1.85,1.85),
                                    y=c(.65,.675,.675,.65)),
                    aes(x=x,y=y)) +
          geom_text(label=signif.num(wilcox.test(work_df[work_df$Group == "PC_REAL_PC_SIM",]$Dist,
                                                 work_df[work_df$Group == "NPC_REAL_NPC_SIM",]$Dist,
                                                 alternative="less")$p.value),
                    x=(1.85 + 1.15)/2,
                    y=.7)
          scale_fill_manual(values=c("#652D90", "#F05A28", "gray20"))
        
        property_plots <- append(property_plots,
                                 list(g))
        
        gpaired = ggplot() +
          base_plot_theme + 
          ggtitle(sprintf("%s - %s", subfield, prop)) +
          theme(text=element_text(size=13.5),
                plot.title=element_text(size=11),
                plot.margin=unit(c(0,0,0,0), "cm")) + 
          ylab(metric_name) +
          xlab("")
        
        
        wilc_prop_df <- 
                    ddply(work_df, 
                         .(Session),  
                         function(session_df) {
                           ddply(session_df,
                                 .(Group),
                                 function(group_df) {
                                   return(c(mean(group_df$Dist),
                                            sd(group_df$Dist),
                                            sem(group_df$Dist)))
                                 })
                         })
        
        for (session in unique(wilc_prop_df$Session)) {

          #for (direc in unique(property_df$Direction)) {
          
          
          wilc_df <- wilc_prop_df[wilc_prop_df$Session == session,]
          
          pcs_df <- wilc_df[wilc_df$Group == "PC_REAL_PC_SIM",]
          npcs_df <- wilc_df[wilc_df$Group == "NPC_REAL_NPC_SIM",]
          
        
          
          gpaired <- gpaired +
            geom_line(data=data.frame(x=c("PCs", "NPCs"),
                                      y=c(pcs_df$V1, npcs_df$V1)),
                      aes(x=x,y=y, group=1), alpha=.3) +
            geom_point(data=data.frame(x=c("PCs", "NPCs"),
                                       y=c(pcs_df$V1, npcs_df$V1)),
                       aes(x=x,y=y, group=1))
          
        }
        
        wilk_test <- wilcox.test(V1~Group, data=wilc_prop_df, alternative="less", correct=F)
        
        max_val <- max(wilc_prop_df$V1)
        
        
        signif_df = data.frame(y=max_val, x=1.5, label=signif.num(wilk_test$p.value))
        
        gpaired <- 
          gpaired + 
          geom_text(data=signif_df, aes(x=x,y=y,label=label))
        
        
        
        compare_boxplots <- append(compare_boxplots, list(gpaired))
        
        pooled_sessions_statistics_df <- 
          
        rbind(pooled_sessions_statistics_df, 
              data.frame(property=prop,
                         metric=metric_name,
                         subfield=subfield,
                         pvalue=wilk_test$p.value,
                         statistic=wilk_test$statistic,
                         alternative=wilk_test$alternative,
                         method=wilk_test$method,
                         signif_code=signif.num(wilk_test$p.value),
                         mean_pc = mean(wilc_prop_df[wilc_prop_df$Group == "PC_REAL_PC_SIM",]$V1),
                         sd_pc = sd(wilc_prop_df[wilc_prop_df$Group == "PC_REAL_PC_SIM",]$V1),
                         sem_pc = sem(wilc_prop_df[wilc_prop_df$Group == "PC_REAL_PC_SIM",]$V1),
                         mean_npc = mean(wilc_prop_df[wilc_prop_df$Group != "PC_REAL_PC_SIM",]$V1),
                         sd_npc = sd(wilc_prop_df[wilc_prop_df$Group != "PC_REAL_PC_SIM",]$V1),
                         sem_npc = sem(wilc_prop_df[wilc_prop_df$Group != "PC_REAL_PC_SIM",]$V1)))
        
      }
      
      property_plots$nrow <- 1
      gf <- do.call(arrangeGrob, property_plots)
      

      
      for (size_name in names(sizes)) {
        size = sizes[[size_name]]
        dir.create(sprintf("%s\\%s\\%s", write_path, size_name, metric_name))
        pdf(sprintf("%s\\%s\\%s\\%s_prop_comp.pdf", write_path, size_name, metric_name, subfield),
            height=size,
            width=size * 4)
        plot(gf)
        dev.off()
      }
      
      subfields_plots <- append(subfields_plots,
                                list(gf))
    }
    
    compare_boxplots$nrow <- 2
    gf_boxplots <- do.call(plot_grid, compare_boxplots)
    
    ap <- align_plots(subfields_plots[[1]],
                      subfields_plots[[2]],
                      align="h",
                      axis="y")
    
    ap$nrow <- 2
    gf_both_rows <- do.call(plot_grid, ap)
    
    for (size_name in names(sizes)) {
      size = sizes[[size_name]]
      dir.create(sprintf("%s\\%s\\%s", write_path, size_name, metric_name))
      pdf(sprintf("%s\\%s\\%s\\all_prop_comp.pdf", write_path, size_name, metric_name),
          height=size * 2,
          width=size * 4)
      plot(gf_both_rows)
      dev.off()
      
      pdf(sprintf("%s\\%s\\%s\\all_prop_comp_boxplots.pdf", write_path, size_name, metric_name),
          height=size * 2,
          width=size * 4)
      plot(gf_boxplots)
      dev.off()
    }
  }
  
  write.csv(file=sprintf("%s//statistics.csv",write_path),
            pooled_sessions_statistics_df)
}

supp_figure_npcs_cyclic_vs_random_diff_sim_comp <- function() {
  
  write_path <- sprintf("%s\\supp_figure_npcs\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\supp_figure_npcs\\cyclic_vs_random_sim_comp\\",figures_path)
  dir.create(write_path)
  
  groups_to_use = c(1:3)
  true_data_paths <- sprintf("%s\\%s", all_data_paths, "equalized")
  file_name="properties.Rda"
  sessions_to_use=c(1:16)
  
  # sizes=c(big=3.5,
  #         big_1=3,
  #         small=2.5,
  #         medium_1=2,
  #         medium_2=1.75,
  #         small=1.5)
  

  
  
  sizes=c(big=3.5,
          medium=3,
          small=2.5)
  
  diff_df <- data.frame()
  
  
  
  
  for (working_path in true_data_paths) {
    
    # extract all paths
    session_paths <- 
      list.dirs(working_path, recursive=F)[grep("session", list.dirs(working_path, recursive=F))]
    
    for (ext in c("", "\\fit_simulation")) { 
      
      # Run through all paths
      for (tpath in session_paths) {
        
        
        idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
        if (sum(idx == sessions_to_use) == 0){
          next
        }
        
        if (len((grep(file_name, list.files(sprintf("%s%s", tpath, ext))))) == 0) {
          print(tpath)
          missing_paths <- c(missing_paths, tpath)
          print("Missing")
          next
        }
        
        if (verbose) {
          print(sprintf("Loading %s", sprintf("%s%s\\%s", tpath, ext, file_name)))
        }
        load(sprintf("%s%s\\%s", tpath, ext, file_name), verbose=T)
        
        
        if (grepl("CA3", working_path)) {
          mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', tpath))
          mice_str <- substr(tpath, mice_str_index, mice_str_index+4) 
        } else {
          mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', tpath))
          mice_str <- substr(tpath, mice_str_index, mice_str_index+3) 
        }
        
        
        subfield <- ifelse(grepl("CA3", working_path), "CA3", "CA1")
        if (grepl("fit_simulation", ext)) {
          params_df <- get_fit_params(tpath)
          sim <- as.character(params_df$pct)
        } else {
          sim <- "Real"
        }
        
        diff_frac <- 
          result$place_cells_metadata$random_fraction - 
          result$place_cells_metadata$cyclic_fraction
          diff_df <- rbind(diff_df,
                           data.frame(diff_frac,
                                ifelse(grepl("Left", tpath), "Left", "Right"), 
                                mice_str, 
                                subfield,
                                sim))

          
      }
    }
  }
  
  gf <- 
  ggplot(diff_df, aes(x=sim, y=diff_frac)) + 
    geom_boxplot(color="gray65", fill=adjustcolor("gray80", alpha=0.8), size=1, width=0.5) +
    xlab("") +
    ylab("Fraction difference (Random - Cyclic)") +
    base_plot_theme +
    theme(text=element_text(size=14, color="black"),
          axis.ticks=element_line(color="black")) 
    
  
  
  for (size_name in names(sizes)) {
    size = sizes[[size_name]]
    dir.create(sprintf("%s\\%s", write_path, size_name))
    pdf(sprintf("%s\\%s\\comp_cyc_rand_diff.pdf", write_path, size_name),
        height=size,
        width=size)
    plot(gf)
    dev.off()
  }
  
  dtr <- dunn.test(diff_df$diff_frac, diff_df$sim, method="bonferroni")
  statistics_df = as.data.frame(do.call(cbind, dtr))
  statistics_df$signif_corr = signif.num(as.numeric(statistics_df$P.adjusted))
  statistics_df$signif = signif.num(as.numeric(statistics_df$P))
  
  write.csv(file=sprintf("%s//dunn_statistics.csv", write_path), statistics_df)
}

supp_figure_npcs_cyclic_vs_random_real_diff_across_learning <- function() {
  write_path <- sprintf("%s\\supp_figure_npcs\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\supp_figure_npcs\\cyc_vs_random_real_across_learning\\",figures_path)
  dir.create(write_path)
  
  dir.create(sprintf("%s\\big\\", write_path))
  dir.create(sprintf("%s\\medium\\", write_path))
  dir.create(sprintf("%s\\small\\", write_path))
  
  sizes=c(big=3.5,
          medium=3,
          small=2.5)
  
  true_data_paths <- sprintf("%s\\%s", all_data_paths, "equalized")
  sessions_to_use=c(1:16)
  
  all_df <- get_all_df(true_data_paths, sessions_to_use)
  
  
  subfield_idx = which(colnames(all_df) == "Subfield")
  session_idx = which(colnames(all_df) == "Session")
  
  percent_active_df <- all_df[,c(7:9, session_idx, subfield_idx)]
  percent_place_cells_df <- all_df[,c(1:3, session_idx, subfield_idx)]
  percent_all_df <- all_df[,c(4:6,session_idx, subfield_idx)]
  
  colnames(percent_active_df) <- c("Both", "Right", "Left", "Session", "Subfield")
  colnames(percent_place_cells_df) <- c("Both", "Right", "Left", "Session", "Subfield")
  colnames(percent_all_df) <- c("Both", "Right", "Left", "Session", "Subfield")
  
  melted_percent_active <- melt(percent_active_df, measure.vars = c("Both", "Left", "Right"))
  melted_percent_place_cells <- melt(percent_place_cells_df, measure.vars = c("Both", "Left", "Right"))
  melted_percent_all <- melt(percent_all_df, measure.vars = c("Both", "Left", "Right"))
  
  ylabs = c("Percentage of active cells (%)",
            "Percentage of place cells (% of active cells)",
            "Percentage of place cells (% of all cells)")
  
  dfs <- list(melted_percent_active,
              melted_percent_place_cells,
              melted_percent_all)
  
  
  estimated_res_cyclic <- get_estimated_df(all_df, true_data_paths, sessions_to_use, estimated_df_folder="KS_simulations_likelihood_c")
  estimated_res_rand <- get_estimated_df(all_df, true_data_paths, sessions_to_use, estimated_df_folder="KS_simulations_likelihood_r")
  
  
  
  
  both_dir_cyc <- get_both_dir_df(estimated_df=estimated_res_cyclic$estimated_df, 
                                  all_df=all_df,
                                  dev_func = sem)
  both_dir_rand <- get_both_dir_df(estimated_df=estimated_res_rand$estimated_df, 
                                   all_df=all_df,
                                   dev_func=sem)
  
  
  
  sessions <- unique(estimated_res_rand$estimated_df$Session)
  directions <- unique(estimated_res_rand$estimated_df$Direction)
  mice <- unique(estimated_res_rand$estimated_df$Mice)
  
  delta_df <- c()
  
  for (mouse in mice) {
    for (session in sessions) {
      for (direction in directions) {
        
        
        cyclic_ind <- 
          estimated_res_cyclic$estimated_df$Session == session &
          estimated_res_cyclic$estimated_df$Direction == direction &
          estimated_res_cyclic$estimated_df$Mice == mouse
        
        rand_ind <- 
          estimated_res_rand$estimated_df$Session == session &
          estimated_res_rand$estimated_df$Direction == direction &
          estimated_res_rand$estimated_df$Mice == mouse
        
        
        if(sum(cyclic_ind) > 0 && sum(rand_ind) > 0) {
          rand_entry <- estimated_res_rand$estimated_df[rand_ind,]
          cyclic_entry <- estimated_res_cyclic$estimated_df[cyclic_ind,]
          delta_entry <- rand_entry
          delta_entry["Estimated"] <- delta_entry["Estimated"] - cyclic_entry["Estimated"]
          names(delta_entry["Estimated"]) <- "Delta"
          delta_df <- rbind(delta_df, delta_entry)
        }
      }
    }
  }
  
  both_dir_df <- rbind(both_dir_cyc$both_dir_df,
                       both_dir_rand$both_dir_df)
  
  mean_sd_df <- rbind(both_dir_cyc$mean_sd_df,
                      both_dir_rand$mean_sd_df)
  
  session_mean_sd_df <- rbind(both_dir_cyc$session_mean_sd_df,
                              both_dir_rand$session_mean_sd_df)
  
  both_dir_df$Shuffle <- c(rep("Cyclic", times=(nrow(both_dir_df) / 2)),
                           rep("Random", times=(nrow(both_dir_df) / 2)))
  
  mean_sd_df$Shuffle <- c(rep("Cyclic", times=(nrow(mean_sd_df) / 2)),
                          rep("Random", times=(nrow(mean_sd_df) / 2)))
  
  session_mean_sd_df$Shuffle <- c(rep("Cyclic", times=(nrow(session_mean_sd_df) / 2)),
                                  rep("Random", times=(nrow(session_mean_sd_df) / 2)))
  
  
  # session_mean_sd_df$Session <- factor(session_mean_sd_df$Session, levels=as.character(1:16))
  # both_dir_df$Session <- factor(both_dir_df$Session, levels=as.character(1:16))
  mean_sd_df <- mean_sd_df[!mean_sd_df$Direction %in% c("Left", "Right"),]
  both_dir_df <- both_dir_df[!both_dir_df$Direction %in% c("Left", "Right"),]
  
  
  session_mean_sd_df <- session_mean_sd_df[session_mean_sd_df$Direction == "Both",]
  both_dir_df <- both_dir_df[both_dir_df$Direction == "Both",]
  
  # session_mean_sd_df$Session <- as.character(as.numeric(session_mean_sd_df$Session) %% 8)
  # session_mean_sd_df$Session[session_mean_sd_df$Session == "0"] <- "8"
  
  delta_both_dir <- get_both_dir_df(estimated_df=delta_df, 
                                    all_df=all_df,
                                    dev_func = sem)
  
  
  session_mean_sd_df <- delta_both_dir$session_mean_sd_df
  both_dir_df <- delta_both_dir$both_dir_df
  session_mean_sd_df <- session_mean_sd_df[session_mean_sd_df$Direction == "Both",]
  both_dir_df <- both_dir_df[both_dir_df$Direction == "Both",]
  
  # session_mean_sd_df$Session <- as.character(as.numeric(session_mean_sd_df$Session) %% 8)
  # session_mean_sd_df$Session[session_mean_sd_df$Session == "0"] <- "8"
  
  
  gdelta <- 
    ggplot(session_mean_sd_df, aes(x=Session, y=Estimated)) + 
    geom_line(aes(group=interaction(Subfield), 
                  color=interaction(Subfield)), 
              size=1, 
              position=position_dodge(0.4)) +
    
    geom_errorbar(aes(ymin=Estimated - Sd, ymax=Estimated + Sd, 
                      group=interaction(Subfield), 
                      color=interaction(Subfield)), 
                  size=.75, 
                  width=1, 
                  position=position_dodge(0.4)) +
    geom_point(aes(group=interaction(Subfield),
                   color=interaction(Subfield)),
               size=2,
               position=position_dodge(0.4)) +
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position="NA",
          text=element_text(size=15),
          plot.title=element_text(size=11))   + 
    scale_color_manual(values=c(CA1="#CE3736", 
                                CA3="#3F5DAB")) +
    xlab("Session") +
    ylab("Fraction") + 
    ggtitle("Unstable cells (estimated)")
  ylim(0,0.4)
  
  for (size_name in names(sizes)) {
    size = sizes[[size_name]]
    
    pdf(file=sprintf("%s\\%s\\delta_plot.pdf",
                     write_path,
                     size_name),
        height=size,
        width=size)
    
    plot(gdelta)
    dev.off()
    
  }
  
  
    delta_statistcs_df <- data.frame()
    dfs_list = list(Rand=both_dir_rand$both_dir_df, Delta=delta_both_dir$both_dir_df)
    for (df_name in names(dfs_list)) {
      
      anova_df = dfs_list[[df_name]]
      anova_df <- anova_df[anova_df$Direction == "Both",]
      
      #anova_df <- work_df[work_df$variable == direc,]
      anova_df$SubID <- 1:nrow(anova_df)
      anova_df$Session <- as.character(as.numeric(anova_df$Session) %% 8)
      anova_df$Session[anova_df$Session == "0"] <- "8"
      
      
      
      ca3_anova_df <- anova_df[anova_df$Subfield == "CA3",]
      ca1_anova_df <- anova_df[anova_df$Subfield == "CA1",]
      
      model.aov <- aov(data=anova_df,
                       formula=Estimated ~ 
                         Subfield * Session + 
                         Error(Mice/(Subfield*Session)))
      
      pvals <- summary(model.aov)
      pvals_df <- as.data.frame(pvals$`Error: Mice:Session`[[1]])
      pvals_df$Group = "Two way"
      pvals_df$Measurement <- df_name
      pvals_df$Conf = naive_conf_name
      pvals_df$Direction = direc
      pvals_df$Factor <- rownames(pvals_df)
      rownames(pvals_df) <- c()
      delta_statistcs_df <- rbind(delta_statistcs_df, pvals_df)
      
      print("2way")
      print(signif.num(pvals$`Error: Mice`[[1]]$`Pr(>F)`[1:2]))
      print(signif.num(pvals$`Error: Mice:Session`[[1]]$`Pr(>F)`[1:2]))
      
      print("CA1")
      model.aov <- aov(data=ca1_anova_df, formula=Estimated~Session + Error(factor(Mice)))
      pvals <- summary(model.aov)
      print(signif.num(pvals$`Error: Within`[[1]]$`Pr(>F)`[1]))
      
      pvals_df <- as.data.frame(pvals$`Error: Within`[[1]])
      pvals_df$Group <- "CA1"
      pvals_df$Measurement <-  df_name
      pvals_df$Conf = naive_conf_name
      pvals_df$Direction = direc
      pvals_df$Factor <- rownames(pvals_df)
      rownames(pvals_df) <- c()
      delta_statistcs_df <- rbind(delta_statistcs_df, pvals_df)
      
      print("CA3")
      model.aov <- aov(data=ca3_anova_df, formula=Estimated~Session + Error(factor(Mice)))
      pvals <- summary(model.aov)
      print(signif.num(pvals$`Error: Within`[[1]]$`Pr(>F)`[1]))
      
      pvals_df <- as.data.frame(pvals$`Error: Within`[[1]])
      pvals_df$Group <- "CA3"
      pvals_df$Measurement <- df_name
      pvals_df$Conf = naive_conf_name
      pvals_df$Direction = direc
      pvals_df$Factor <- rownames(pvals_df)
      rownames(pvals_df) <- c()
      delta_statistcs_df <- rbind(delta_statistcs_df, pvals_df)
      
    }
    delta_statistcs_df$signif = signif.num(delta_statistcs_df$`Pr(>F)`)
    
    tmp_df <- rbind(both_dir_rand$session_mean_sd_df,
                    delta_both_dir$session_mean_sd_df)
    tmp_df$Group <- rep(c("Random", "Delta"), each=nrow(both_dir_rand$session_mean_sd_df))

    write.csv(file=sprintf("%s//delta_statistics.csv", write_path),
              delta_statistcs_df)
    write.csv(file=sprintf("%s//values.csv", write_path),
              tmp_df)
  
}
