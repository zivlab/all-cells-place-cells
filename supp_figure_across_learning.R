library(pheatmap)
library(networkD3)

supp_figure_naive_learning_place_cells_all_sessions <- function()
{
  write_path <- sprintf("%s\\supp_figure_naive_learning\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\supp_figure_naive_learning\\pie_plots\\",figures_path)
  dir.create(write_path)
  
  directions_to_use = c("Both", "Left", "Right")
  sizes=c(big=2.5,
          medium=2,
          medium_1=1.75,
          medium_2=1.5,
          small=1.25)
  
  df_res <- c()
  df_non_pcs_sessions <- c()
  path_list <- all_data_paths
  for (path_idx in 1:len(mice_ids)) {
    
    first_path <- all_data_paths[(path_idx * 2 - 1)]
    second_path <- all_data_paths[(path_idx * 2)]
    
    mat_path <- sprintf("%s\\%s", first_path, "registered_cells_categories.Rda")
    load(mat_path)
    
    cells_mat_both_directions <- cells_mat
    
    mat_path <- sprintf("%s\\%s", second_path, "registered_cells_categories.Rda")
    load(mat_path)
    cells_mat_both_directions <- cbind(cells_mat_both_directions,
                                       cells_mat)
    work <- cells_mat_both_directions
    for (env_ind in list(c(1:8,17:24), c(9:16,25:32))) {

    cells_mat_both_directions <- work[,env_ind]
    most_advent_status <- unlist(lapply(apply(cells_mat_both_directions, 1, table), function(tb) {max(as.numeric(names(tb)))}))
    non_pcs <- which(most_advent_status == 1)
    
    non_pcs_sessions_distribution <- 
      colMeans(t(apply(cells_mat_both_directions[non_pcs,], 1, function(cell) {c(sum(cell == -1), sum(cell == 0), sum(cell == 1))})))
    non_pcs_sessions_distribution <- non_pcs_sessions_distribution / sum(non_pcs_sessions_distribution)
    
    # measurement of which are the sessions that cells that were at best non-pcs
    sessions_of_non_pcs <- table(unlist(apply(cells_mat_both_directions[non_pcs,], 1, 
                                              function(r) {
                                                tmp <- which(r==1) %% 8; 
                                                tmp[tmp ==0] <- 8; 
                                                rep(0, times=8)
                                                return(tmp)})))
    
    final_percentages <- 
      c(sum(most_advent_status == 1 | most_advent_status == 2),
        sum(most_advent_status == 3 | most_advent_status == 4))
    
    final_percentages  <- final_percentages / sum(final_percentages)
    
    df_res <- rbind(df_res,
                    c(final_percentages,
                      non_pcs_sessions_distribution,
                      path_idx,
                      nrow(cells_mat_both_directions)))
    
    
    ntmp <- rep(0, times=8)
    names(ntmp) <- 1:8
    ntmp[names(sessions_of_non_pcs)] <- sessions_of_non_pcs
    
    df_non_pcs_sessions <- rbind(df_non_pcs_sessions,
                                 ntmp)
    }
  }
  
  blank_theme <-   theme_minimal()+
    theme(axis.text.x=element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank(),
          axis.ticks = element_blank(),
          plot.title=element_text(size=14, face="bold"))
  
  blank_theme_no_legend <-   theme_minimal()+
    theme(axis.text.x=element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank(),
          axis.ticks = element_blank(),
          plot.title=element_text(size=14, face="bold"),
          legend.position="NA")  
  ca1_df <- df_res[1:8,]
  ca3_df <- df_res[9:18,]
  
  ca1_df_pcs <- t(apply(ca1_df, 1, function(mice) {mice[1:2] * mice[7]}))
  ca1_pcs_all <- colSums(ca1_df_pcs) / sum(ca1_df_pcs)
  
  ca3_df_pcs <- t(apply(ca3_df, 1, function(mice) {mice[1:2] * mice[7]}))
  ca3_pcs_all <- colSums(ca3_df_pcs) / sum(ca3_df_pcs)
  
  ca1_df_sessions_non_pcs <- t(apply(ca1_df, 1, function(mice) {mice[3:5]}))
  ca3_df_sessions_non_pcs <- t(apply(ca3_df, 1, function(mice) {mice[3:5]}))
  
  
  ca1_non_pcs_sessions <- colSums(df_non_pcs_sessions[1:8,])
  ca3_non_pcs_sessions <- colSums(df_non_pcs_sessions[9:18,])
  
  ca1_plot_df <- data.frame(Fraction=ca1_pcs_all,
                            Labels=c(#"Insuffieciently active",
                                     "Non place cells",
                                     "Place cells"))
  
  ca1_plot_df_sessions <- data.frame(Fraction=colMeans(ca1_df_sessions_non_pcs),
                                     Labels=c("Silent",
                                              "Inactive",
                                              "Non place cells"))
  
  ca1_plot_df_npcs_which_sessions <- 
    data.frame(Fraction=c(sum(ca1_non_pcs_sessions[1:4]), 
                          sum(ca1_non_pcs_sessions[5:8])) / sum(ca1_non_pcs_sessions),
               Labels=c("First four sessions","Last four sessions"))
  
  ca3_plot_df <- data.frame(Fraction=ca3_pcs_all,
                            Labels=c(#"Insuffieciently active",
                                     "Non place cells",
                                     "Place cells"))
  
  ca3_plot_df_sessions <- data.frame(Fraction=colMeans(ca3_df_sessions_non_pcs),
                                     Labels=c("Silent",
                                              "Inactive",
                                              "Non place cells"))
  
  ca3_plot_df_npcs_which_sessions <- 
    data.frame(Fraction=c(sum(ca3_non_pcs_sessions[1:4]), 
                          sum(ca3_non_pcs_sessions[5:8])) / sum(ca3_non_pcs_sessions),
               Labels=c("First four sessions","Last four sessions"))
  
  gfraction_ca1 <- 
  ggplot(ca1_plot_df, aes(x="", y=Fraction, fill=Labels)) +
    geom_bar(width=1, stat="identity") +
    coord_polar("y", start=+180) + 
    geom_text(aes(y = Fraction/4 + c(0, cumsum(Fraction)[-length(Fraction)]), 
                  label = round(Fraction * 100)), size=5) +
    scale_fill_grey() + 
    blank_theme
  
  gnpc_sessions_ca1 <- 
    ggplot(ca1_plot_df_sessions, aes(x="", y=Fraction, fill=Labels)) +
    geom_bar(width=1, stat="identity") +
    coord_polar("y", start=+180) + 
    geom_text(aes(y = Fraction/4 + c(0, cumsum(Fraction)[-length(Fraction)]), 
                  label = round(Fraction * 100)), size=5) +
    scale_fill_manual(values=c("#652D90", "#F05A28", "gray20"),
                      breaks=c("Inactive", "Silent", "Non place cells"))
  
  gplot_df_npcs_which_sessions_ca1 <- 
    ggplot(ca1_plot_df_npcs_which_sessions, aes(x="", y=Fraction, fill=Labels)) +
    geom_bar(width=1, stat="identity") +
    coord_polar("y", start=+180) + 
    geom_text(aes(y = Fraction/4 + c(0, cumsum(Fraction)[-length(Fraction)]), 
                  label = round(Fraction * 100)), size=5) +
    scale_fill_manual(values=c("#3F5DAB", "#CE3736"))
  
  gfraction_ca3 <- 
    ggplot(ca3_plot_df, aes(x="", y=Fraction, fill=Labels)) +
    geom_bar(width=1, stat="identity") +
    coord_polar("y", start=+180) + 
    geom_text(aes(y = Fraction/4 + c(0, cumsum(Fraction)[-length(Fraction)]), 
                  label = round(Fraction * 100)), size=5) +
    scale_fill_grey()
  
  gnpc_sessions_ca3 <- 
    ggplot(ca3_plot_df_sessions, aes(x="", y=Fraction, fill=Labels)) +
    geom_bar(width=1, stat="identity") +
    coord_polar("y", start=+180) + 
    geom_text(aes(y = Fraction/4 + c(0, cumsum(Fraction)[-length(Fraction)]), 
                  label = round(Fraction * 100)), size=5) +
    scale_fill_manual(values=c("#652D90", "#F05A28", "gray20"),
                      breaks=c("Inactive", "Silent", "Non place cells"))
  
  gplot_df_npcs_which_sessions_ca3 <- 
    ggplot(ca3_plot_df_npcs_which_sessions, aes(x="", y=Fraction, fill=Labels)) +
    geom_bar(width=1, stat="identity") +
    coord_polar("y", start=+180) + 
    geom_text(aes(y = Fraction/4 + c(0, cumsum(Fraction)[-length(Fraction)]), 
                  label = round(Fraction * 100)), size=5) +
    scale_fill_manual(values=c("#3F5DAB", "#CE3736"))
  
  
  
  
  plots <- list(list(plot=gfraction_ca1,
                     name="PCS_fraction_CA1"),
                
                list(plot=gnpc_sessions_ca1,
                     name="NPCS_sessions_dist_CA1"),
                
                list(plot=gplot_df_npcs_which_sessions_ca1,
                     name="NPCS_which_sessions_CA1"),
                
                list(plot=gfraction_ca3,
                     name="PCS_fraction_CA3"),
                
                list(plot=gnpc_sessions_ca3,
                     name="NPCS_sessions_dist_CA3"),
                
                list(plot=gplot_df_npcs_which_sessions_ca3,
                     name="NPCS_which_sessions_CA3"))
  
  
  ca1_tree_plot <- 
    ggplot() +
    geom_line(aes(x=x,y=y), data=data.frame(x=c(1,1,2,2), y=c(5,6,6,7))) +
    geom_line(aes(x=x,y=y), data=data.frame(x=c(2,3,3,4,4), y=c(6,6,4.5,4.5,3.5))) +
    geom_line(aes(x=x,y=y), data=data.frame(x=c(2,2,3,3), y=c(3,4.5,4.5,3.5))) +
    #geom_line(aes(x=x,y=y), data=data.frame(x=c(1,1,2), y=c(2,3,3))) + 
    geom_line(aes(x=x,y=y), data=data.frame(x=c(1.5,1.5,2.5,2.5), y=c(2,3,3,2))) +
    geom_line(aes(x=x,y=y), data=data.frame(x=c(1.5,1.5,2,2,3,3), 
                                            y=c(2,3,3,4.5,4.5,6)),
              col="red") +
    geom_line(aes(x=x,y=y), data=data.frame(x=c(2,2,3), 
                                            y=c(7,6,6)),
              col="red") +      
    geom_text(data=data.frame(x=1.35, 
                              y=5, 
                              label=sprintf("%.3f",ca1_plot_df[ca1_plot_df$Labels == "Place cells",]$Fraction)), 
              aes(x=x,y=y,label=label)) +
    
    geom_text(data=data.frame(x=3.35, 
                              y=5, 
                              label=sprintf("%.3f",ca1_plot_df[ca1_plot_df$Labels != "Place cells",]$Fraction)), 
              aes(x=x,y=y,label=label)) + 
    geom_text(data=data.frame(x=3.65, 
                              y=3.5, 
                              label=sprintf("%.3f",ca1_plot_df_sessions[ca1_plot_df_sessions$Labels == "Silent",]$Fraction)), 
              aes(x=x,y=y,label=label)) +
    geom_text(data=data.frame(x=2.65, 
                              y=3.5, 
                              label=sprintf("%.3f",ca1_plot_df_sessions[ca1_plot_df_sessions$Labels == "Inactive",]$Fraction)), 
              aes(x=x,y=y,label=label)) +
    geom_text(data=data.frame(x=1.65, 
                              y=3.5, 
                              label=sprintf("%.3f",ca1_plot_df_sessions[ca1_plot_df_sessions$Labels == "Non place cells",]$Fraction)), 
              aes(x=x,y=y,label=label)) +
    geom_text(data=data.frame(x=1.15, 
                              y=2, 
                              label=sprintf("%.3f",ca1_plot_df_npcs_which_sessions[ca1_plot_df_npcs_which_sessions$Labels != "First four sessions",]$Fraction)), 
              aes(x=x,y=y,label=label)) +
    geom_text(data=data.frame(x=2.85, 
                              y=2, 
                              label=sprintf("%.3f",ca1_plot_df_npcs_which_sessions[ca1_plot_df_npcs_which_sessions$Labels == "First four sessions",]$Fraction)), 
              aes(x=x,y=y,label=label)) + 
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title=element_blank(),
          axis.ticks=element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position="NA",
          text=element_text(size=15),
          plot.title=element_text(size=11))
  
  
  ca3_tree_plot <- 
    ggplot() +
    geom_line(aes(x=x,y=y), data=data.frame(x=c(1,1,2,2), y=c(5,6,6,7))) +
    geom_line(aes(x=x,y=y), data=data.frame(x=c(2,3,3,4,4), y=c(6,6,4.5,4.5,3.5))) +
    geom_line(aes(x=x,y=y), data=data.frame(x=c(2,2,3,3), y=c(3,4.5,4.5,3.5))) +
    #geom_line(aes(x=x,y=y), data=data.frame(x=c(1,1,2), y=c(2,3,3))) + 
    geom_line(aes(x=x,y=y), data=data.frame(x=c(1.5,1.5,2.5,2.5), y=c(2,3,3,2))) +
    geom_line(aes(x=x,y=y), data=data.frame(x=c(1.5,1.5,2,2,3,3), 
                                            y=c(2,3,3,4.5,4.5,6)),
              col="red") +
    geom_line(aes(x=x,y=y), data=data.frame(x=c(2,2,3), 
                                            y=c(7,6,6)),
              col="red") +      
    geom_text(data=data.frame(x=1.35, 
                              y=5, 
                              label=sprintf("%.3f",ca3_plot_df[ca3_plot_df$Labels == "Place cells",]$Fraction)), 
              aes(x=x,y=y,label=label)) +
    
    geom_text(data=data.frame(x=3.35, 
                              y=5, 
                              label=sprintf("%.3f",ca3_plot_df[ca3_plot_df$Labels != "Place cells",]$Fraction)), 
              aes(x=x,y=y,label=label)) + 
    geom_text(data=data.frame(x=3.65, 
                              y=3.5, 
                              label=sprintf("%.3f",ca3_plot_df_sessions[ca3_plot_df_sessions$Labels == "Silent",]$Fraction)), 
              aes(x=x,y=y,label=label)) +
    geom_text(data=data.frame(x=2.65, 
                              y=3.5, 
                              label=sprintf("%.3f",ca3_plot_df_sessions[ca3_plot_df_sessions$Labels == "Inactive",]$Fraction)), 
              aes(x=x,y=y,label=label)) +
    geom_text(data=data.frame(x=1.65, 
                              y=3.5, 
                              label=sprintf("%.3f",ca3_plot_df_sessions[ca3_plot_df_sessions$Labels == "Non place cells",]$Fraction)), 
              aes(x=x,y=y,label=label)) +
    geom_text(data=data.frame(x=1.15, 
                              y=2, 
                              label=sprintf("%.3f",ca3_plot_df_npcs_which_sessions[ca3_plot_df_npcs_which_sessions$Labels != "First four sessions",]$Fraction)), 
              aes(x=x,y=y,label=label)) +
    geom_text(data=data.frame(x=2.85, 
                              y=2, 
                              label=sprintf("%.3f",ca3_plot_df_npcs_which_sessions[ca3_plot_df_npcs_which_sessions$Labels == "First four sessions",]$Fraction)), 
              aes(x=x,y=y,label=label)) + 
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title=element_blank(),
          axis.ticks=element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position="NA",
          text=element_text(size=15),
          plot.title=element_text(size=11))
  
  
  
  
  for (size_name in names(sizes)) {
    size = sizes[[size_name]]
    
    dir.create(sprintf("%s\\%s",
                       write_path,
                       size_name))
    
    for (p in plots) {
      pdf(sprintf("%s\\%s\\%s.pdf",
                  write_path,
                  size_name,
                  p$name),
          height=size,
          width=size)
      plot(p$plot + blank_theme_no_legend)
      dev.off()
      
      pdf(sprintf("%s\\%s\\legend_%s.pdf",
                  write_path,
                  size_name,
                  p$name),
          height=size,
          width=size)
      plot(p$plot + blank_theme)
      dev.off()
      
      pdf(sprintf("%s\\%s\\tree_plot_ca1.pdf",
                  write_path,
                  size_name),
          height=size,
          width=size)
      plot(ca1_tree_plot)
      dev.off()
      
      pdf(sprintf("%s\\%s\\tree_plot_ca3.pdf",
                  write_path,
                  size_name),
          height=size,
          width=size)
      plot(ca3_tree_plot)
      dev.off()
      
    }
  }
}


supp_figure_naive_learning_lineplots_equal_traversals  <- function() 
{
  
  write_path <- sprintf("%s\\supp_figure_naive_learning\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\supp_figure_naive_learning\\lineplots_equal_traversals\\",figures_path)
  dir.create(write_path)
  sizes=c(big=3.5,
          medium=3,
          small=2.5)
  
  dir.create(write_path)
  
  true_data_paths <- sprintf("%s\\%s", all_data_paths[1:18], "equalized")
  
  sessions_to_use=c(1:16)
  
  
  directions_to_use = c("Both")
  
  naive_statistics_df <- data.frame()
  values_df <- data.frame()
  
  naive_configurations <- list(cyclic_eq_prior_either_dir=list(path="multiple_pvals_percent_df_cyclic_SI_eq_trav.R", 
                                                               weighted=F),
                               cyclic_eq_prior_weighted=list(path="multiple_pvals_percent_df_cyclic_SI_eq_trav.R", 
                                                             weighted=T))
  
  
  for (naive_conf_name in names(naive_configurations)) {
    naive_conf <- naive_configurations[[naive_conf_name]]
    
    all_df <- get_all_df(true_data_paths, 
                         sessions_to_use,
                         file_name = naive_conf$path,
                         multiple_pvals = T)
    
    all_df <- all_df$`0.05`
    
    estimated_res <- get_estimated_df(all_df, true_data_paths, sessions_to_use)
    
    both_dir_res <- get_both_dir_df(estimated_df = estimated_res$estimated_df, all_df=all_df, dev_func = sem)
    
    both_dir_df <- both_dir_res$both_dir_df
    mean_sd_df <- both_dir_res$mean_sd_df
    session_mean_sd_df <- both_dir_res$session_mean_sd_df
    
    mean_sd_df <- mean_sd_df[mean_sd_df$Direction %in% directions_to_use, ]
    both_dir_df <- both_dir_df[both_dir_df$Direction %in% directions_to_use, ]
    
    both_dir_df_all <- both_dir_res$both_dir_df_all
    session_mean_sd_df_all <- both_dir_res$session_mean_sd_df_al
    
    percent_active_df <- all_df[,c(7:9)]
    percent_place_cells_df <- all_df[,c(1:3)]
    percent_all_df <- all_df[,c(4:6)]
    
    metadata_df <- all_df[,c("Session", "Mice", "Subfield")]
    percent_active_df <- cbind(percent_active_df, metadata_df)
    percent_all_df <- cbind(percent_all_df, metadata_df)
    
    percent_place_cells_df <- cbind(percent_place_cells_df, metadata_df)  
    colnames(percent_active_df) <- c("Both", "Right", "Left", "Session", "Mice", "Subfield")
    colnames(percent_place_cells_df) <- c("Both", "Right", "Left", "Session", "Mice", "Subfield")
    colnames(percent_all_df) <- c("Both", "Right", "Left", "Session", "Mice", "Subfield")
    
    melted_percent_active <- melt(percent_active_df, id.vars = c("Session", "Mice", "Subfield"))
    melted_percent_place_cells <- melt(percent_place_cells_df, id.vars = c("Session", "Mice", "Subfield"))
    melted_percent_all <- melt(percent_all_df, id.vars = c("Session", "Mice", "Subfield"))
    
    ylabs = c("Percentage of active cells (%)",
              "Percentage of place cells (% of active cells)",
              "Percentage of place cells (% of all cells)")
    
    dfs <- list(melted_percent_active,
                melted_percent_place_cells,
                melted_percent_all)
    
    weighted_avg_both_df_colnames <- c("Active",
                                       "Of_Active",
                                       "Of_All")
    
    plot_names <- c("ActiveCells",
                    "PlaceCellsOfActive",
                    "PlaceCellsOfAll")
    
    
    plot_titles <- c("Active cells",
                     "Place cells (of active)",
                     "Place cells (of all)")
    
    ylabs <- c("Fraction",
               "Naive fraction",
               "Naive fraction")
    
    
    for (i  in 1:len(dfs)) {
      df_idx <- i
      pn <- weighted_avg_both_df_colnames[i]
      
      work_df <- dfs[[i]]
      
      
      for(direc in unique(work_df$variable)) {
        
        
        anova_df <- work_df[work_df$variable == direc,]
        anova_df$SubID <- 1:nrow(anova_df)
        anova_df$Session <- as.character(as.numeric(anova_df$Session) %% 8)
        anova_df$Session[anova_df$Session == "0"] <- "8"
        colnames(anova_df) <- c("Session","Mice", "Subfield", "#", "Frac", "SubID")
        
        
        ca3_anova_df <- anova_df[anova_df$Subfield == "CA3",]
        ca1_anova_df <- anova_df[anova_df$Subfield == "CA1",]
        
        model.aov <- aov(data=anova_df,
                         formula=Frac ~ 
                           Subfield * Session + 
                           Error(Mice/(Subfield*Session)))
        
        pvals <- summary(model.aov)
        pvals_df <- as.data.frame(pvals$`Error: Mice:Session`[[1]])
        pvals_df$Group = "Two way"
        pvals_df$Measurement <- plot_titles[i]
        pvals_df$Conf = naive_conf_name
        pvals_df$Direction = direc
        pvals_df$Factor <- rownames(pvals_df)
        rownames(pvals_df) <- c()
        naive_statistics_df <- rbind(naive_statistics_df, pvals_df)
        
        print("2way")
        print(signif.num(pvals$`Error: Mice`[[1]]$`Pr(>F)`[1:2]))
        print(signif.num(pvals$`Error: Mice:Session`[[1]]$`Pr(>F)`[1:2]))
        
        print("CA1")
        model.aov <- aov(data=ca1_anova_df, formula=Frac~Session + Error(factor(Mice)))
        pvals <- summary(model.aov)
        print(signif.num(pvals$`Error: Within`[[1]]$`Pr(>F)`[1]))
        
        pvals_df <- as.data.frame(pvals$`Error: Within`[[1]])
        pvals_df$Group <- "CA1"
        pvals_df$Measurement <- plot_titles[i]
        pvals_df$Conf = naive_conf_name
        pvals_df$Direction = direc
        pvals_df$Factor <- rownames(pvals_df)
        rownames(pvals_df) <- c()
        naive_statistics_df <- rbind(naive_statistics_df, pvals_df)
        
        print("CA3")
        model.aov <- aov(data=ca3_anova_df, formula=Frac~Session + Error(factor(Mice)))
        pvals <- summary(model.aov)
        print(signif.num(pvals$`Error: Within`[[1]]$`Pr(>F)`[1]))
        
        pvals_df <- as.data.frame(pvals$`Error: Within`[[1]])
        pvals_df$Group <- "CA3"
        pvals_df$Measurement <- plot_titles[i]
        pvals_df$Conf = naive_conf_name
        pvals_df$Direction = direc
        pvals_df$Factor <- rownames(pvals_df)
        rownames(pvals_df) <- c()
        naive_statistics_df <- rbind(naive_statistics_df, pvals_df)
      }
      
      if (naive_conf$weighted) {
        session_mean_work_df <- session_mean_sd_df_all[,c("Subfield", "Session", "Direction", sprintf("%s_%s", pn, c("mean", "sd")))]
        colnames(session_mean_work_df) <- c("Subfield", "Session", "Direction", "Mean", "Sd")
        value_df <- session_mean_work_df[,c("Subfield", "Session", "Mean", "Sd")]
        
      } else {
        work_df <- dfs[[i]]
        work_df <- work_df[work_df$variable %in% directions_to_use,]
        work_df$Session <- as.character(as.numeric(work_df$Session) %% 8)
        work_df$Session[work_df$Session == "0"] <- "8"
        
        session_mean_work_df <- 
          ddply(work_df, .(Subfield), function(subfield_df) {
            both_df  <- ddply(subfield_df, .(Session), function(session_df){
              
              session_df
              return(c(mean(session_df[,"value"]), 
                       sd(session_df[,"value"])))
            })
          })
        colnames(session_mean_work_df) <- c("Subfield", "Session", "Mean", "Sd")
        value_df <- session_mean_work_df
      }
      
      
      value_df$title <-  plot_titles[i]
      value_df$conf <- naive_conf_name
      
      values_df <- rbind(values_df,
                         value_df)
      
      g <- 
        ggplot(session_mean_work_df, aes(x=Session, y=Mean)) + 
        geom_line(aes(group=Subfield, color=Subfield), size=1, position=position_dodge(0.4)) +
        geom_errorbar(aes(ymin=Mean - Sd, ymax=Mean + Sd, group=Subfield, color=Subfield), 
                      size=.75, width=1, position=position_dodge(0.4)) +
        geom_point(aes(group=Subfield, color=Subfield), size=2, position=position_dodge(0.4)) + 
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
        ylab(ylabs[i]) + 
        ylim(0,1.1) + 
        ggtitle(plot_titles[i])
      
      
      
      for (size_name in names(sizes)) {
        size = sizes[[size_name]]
        
        dir.create(sprintf("%s\\%s",
                           write_path,
                           size_name))
        
        dir.create(sprintf("%s\\%s\\%s\\",
                           write_path,
                           size_name,
                           plot_names[[df_idx]]))
        
        pdf(file=sprintf("%s\\%s\\%s\\%s.pdf",
                         write_path,
                         size_name,
                         plot_names[[df_idx]],
                         naive_conf_name),
            height=size,
            width=size)
        
        plot(g)
        dev.off()
      }
    }
  }
  
  naive_statistics_df$signif_code <- signif.num(naive_statistics_df$`Pr(>F)`)
  
  write.csv(naive_statistics_df,
            file=sprintf("%s\\naive_measurements_statistics.csv",write_path))
  write.csv(values_df,
            file=sprintf("%s\\naive_measurements_values.csv",write_path))
}

supp_figure_naive_learning_npcs_across_days <- function()
{
  write_path <- sprintf("%s\\supp_figure_naive_learning\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\supp_figure_naive_learning\\fraction_plots\\",figures_path)
  dir.create(write_path)
  
  sizes=c(big=3,
          medium=2.5,
          medium_1=2,
          medium_2=1.75,
          small=1.5)
  
  plt_list <- list()
  path_list <- list(CA1=all_data_paths[1:8],
                    CA3=all_data_paths[9:17])
  
  switchers = F
  
  DISAPPEARED = "-1";
  DORMENT = "0";
  NON_RANDOM_NON_CYCLIC = "1";
  RANDOM_NON_CYLIC = "2"; # Swithcers
  NON_RANDOM_CYCLIC = "3";
  RANDOM_CYCLIC = "4";
  
  categoreis_all <- c(DISAPPEARED,
                      DORMENT,
                      NON_RANDOM_NON_CYCLIC,
                      RANDOM_NON_CYLIC,
                      NON_RANDOM_CYCLIC,
                      RANDOM_CYCLIC)
  
  
  all_tabulation <- c()
  for (i1 in 1:len(categoreis_all)) {
    for (i2 in 1:len(categoreis_all)) {
      all_tabulation <- c(all_tabulation, paste(categoreis_all[i1], categoreis_all[i2], sep="-"))
    } 
  }
  
  
  
  for (subfield in names(path_list)) {
    work_paths = path_list[[subfield]]
    tabulated_vec <- rep(0, times=len(categoreis_all) ** 2)
    names(tabulated_vec) <- all_tabulation
    
    final_mat <- c();
    
    for (p in work_paths) {
      
      mat_path <- sprintf("%s\\%s", p, "registered_cells_categories.Rda")
      load(mat_path)
      
      tabulated_mice_dir_transitions <- 
        table(unlist(lapply(1:(ncol(cells_mat) - 1), 
                            function(idx) {c(paste(cells_mat[,idx], 
                                                   cells_mat[,idx +1], sep="-"))})))
      tabulated_vec[names(tabulated_mice_dir_transitions)] <- 
        tabulated_vec[names(tabulated_mice_dir_transitions)] + 
        tabulated_mice_dir_transitions
      
      full_cells_mat  <- cells_mat
      for (session in list(c(1:8), c(9:16))) {
        
        cells_indices_list <- list()
        
        
        for (idx in rev(session)) {
          i <- idx
          print(sprintf("Running from %d to %d", i, min(session)))
          
          if (switchers) {
            npc_group = c(2)
          } else {
            npc_group = c(1,2)
          }
          ind <- cells_mat[,i] %in% npc_group
          print(sprintf("Number of non place cells in session %d : %d", i, sum(ind)))
          i <- i - 1
          
          while (i >= min(session)) {
            
            ind <- ind & cells_mat[,i] == -1
            print(sprintf("- That did not appear on session %d %d", i ,
                          sum(ind)))
            
            i <- i - 1
            
            
          }
          
          print("######")
          print("-")
          print(which(ind))
          print("-")
          print("######")
          
          cells_indices_list <- append(list(list(ind=which(ind), session=idx)),
                                       cells_indices_list)
        }
        
        
        #### MAKE SURE THEY ARE ALL ORTHOGONAL!!! 
        for(i in 1:len(cells_indices_list)) {
          for(j in 1:len(cells_indices_list)) {
            if(i != j) {
              assert(len(intersect(cells_indices_list[[i]]$ind, 
                                   cells_indices_list[[j]]$ind)) == 0)
            }
          }
        }
        
        for (npc_indices in cells_indices_list) { 
          
          if (npc_indices$session == max(session)) {
            print(sprintf("Session %d continuing!", npc_indices$session))
            next
          }
          
          if(len(npc_indices$ind) <= 1) {
            print("No new non-place cells")
            next
          }
          
          npc_mat <- cells_mat[npc_indices$ind, npc_indices$session:max(session)]
          #tabled <- apply(npc_mat, 2, table)
          tabled <- lapply(1:ncol(npc_mat), function(col_idx){table(npc_mat[,col_idx])})
          
          
          res <- 
            lapply(1:len(tabled),
                   function(idx) {
                     fractions <- tabled[[idx]]
                     disappeared <- 0
                     inactive <- 0
                     non_place_cells <- 0
                     switchers  <- 0
                     place_cells <- 0
                     
                     if (DISAPPEARED %in% names(fractions)) {
                       disappeared <- disappeared + fractions[DISAPPEARED]
                     } 
                     
                     if (DORMENT %in% names(fractions)) {
                       inactive <- inactive + fractions[DORMENT]
                     } 
                     
                     if (NON_RANDOM_NON_CYCLIC %in% names(fractions)) {
                       non_place_cells <- non_place_cells + fractions[NON_RANDOM_NON_CYCLIC]
                     } 
                     
                     if (RANDOM_NON_CYLIC %in% names(fractions)) {
                       if (switchers) {
                         switchers <- switchers + fractions[RANDOM_NON_CYLIC]
                       } else {
                         non_place_cells <- non_place_cells + fractions[RANDOM_NON_CYLIC]
                       }
                     }
                     
                     if (NON_RANDOM_CYCLIC %in% names(fractions)) {
                       place_cells <- place_cells + fractions[NON_RANDOM_CYCLIC]
                     }
                     
                     if (RANDOM_CYCLIC %in% names(fractions)) {
                       place_cells <- place_cells + fractions[RANDOM_CYCLIC]
                     } 
                     
                     
                     if (switchers) {
                       fracs <- c(disappeared,
                                  inactive,
                                  place_cells,
                                  switchers,
                                  non_place_cells)
                     } else {
                       fracs <- c(disappeared,
                                  inactive,
                                  place_cells,
                                  non_place_cells)          
                     }
                     
                     fracs <- fracs / sum(fractions)
                     
                     res <- c(fracs, idx, sum(fractions))
                     print("###############")
                     #print(fracs)
                     print(">>>>>>>>>>>>")
                     print(res)
                     print("<<<<<<<<<<<<<<<<<")
                     
                     
                     return(res)
                     
                   })
          
          res <- do.call(rbind, res)
          if (switchers) {
            colnames(res) <- c("Disappeared", "Inactive", "Place cells", "Switchers", "Non place cells", "Session", "Num of cells")
          } else {
            colnames(res) <- c("Disappeared", "Inactive", "Place cells", "Non place cells", "Session", "Num of cells")
          }
          
          final_mat <- rbind(final_mat,
                             res)
        }
      }
      
    }
    
    final_df <- as.data.frame(final_mat)
    
    if (switchers) {
      n_groups = 5
    } else {
      n_groups = 4
    }
    df <- ddply(final_df, .(Session), 
                function(sub_df) {
                  #a <- colMeans(sub_df)
                  weighted <- t(apply(sub_df, 1, function(r) { return(r[1:n_groups] * r["Num of cells"])}))
                  return(colSums(weighted) / sum(sub_df[,"Num of cells"]))
                })
    
    melted <- melt(df[,1:(n_groups + 1)], id.vars = "Session")
    melted$Session <- melted$Session - 1
    colnames(melted) <- c("Session", "Group", "Fraction")
    g <- 
      ggplot(melted, aes(x=Session, y=Fraction, fill=Group)) + 
      geom_bar(position="stack", stat="identity", color="black", size=0.25, width=1) + 
      scale_fill_brewer(palette="Blues") +
      xlab("Time elapsed (sessions)") + theme_classic() +
      ggtitle(subfield)
    
    plt_list <- append(plt_list,
                       list(g))
    
    
    ### Transition mats
    
    prob_transitions <- tabulated_vec / sum(tabulated_vec)
    transition_mat <- matrix(rep(0, times=len(categoreis_all) ** 2),
                             nrow=len(categoreis_all))
    
    simplify_category <- function(cat) 
    {ifelse(cat == 0, "Inactive", ifelse(cat == 1 | cat == 2, "NPC", "PC"))}
    
    simple_categoreis <- c("Inactive", "NPC", "PC")
    simple_transition_mat <- matrix(rep(0, times=len(simple_categoreis) ** 2),
                                    nrow=len(simple_categoreis))
    
    colnames(simple_transition_mat) <- simple_categoreis
    rownames(simple_transition_mat) <- simple_categoreis
    
    
    colnames(transition_mat) <- categoreis_all
    rownames(transition_mat) <- categoreis_all
    
    
    
    for (i1 in 1:len(categoreis_all)) {
      for (i2 in 1:len(categoreis_all)) {
        cat_i1 <- categoreis_all[i1]
        cat_i2 <- categoreis_all[i2]
        transition_name <- paste(cat_i1,
                                 cat_i2, 
                                 sep="-")
        transition_mat[cat_i1, cat_i2] <- tabulated_vec[transition_name]
        
        simple_transition_mat[simplify_category(cat_i1),
                              simplify_category(cat_i2)] <- 
          simple_transition_mat[simplify_category(cat_i1),
                                simplify_category(cat_i2)] + 
          tabulated_vec[transition_name]
      } 
    }
    
    complex_prob <- transition_mat / sum(transition_mat)
    simple_prob <- simple_transition_mat / sum(simple_transition_mat)
    
    pheatmap(simple_prob , cluster_rows=F, cluster_cols=F, breaks=seq(0,0.1, length.out=100), col=rev(blc(100)))
  }
  
  
  plt_list_no_legend <- lapply(plt_list, function(p) {p + base_plot_theme})
  plt_list$nrow <- 1
  plt_list_no_legend$nrow <- 1
  gf <- do.call(grid.arrange, plt_list)
  gf_no_legend <- do.call(grid.arrange, plt_list_no_legend)
  
  
  for (size_name in names(sizes)) {
    
    dir.create(sprintf("%s\\%s",
                       write_path, 
                       
                       size_name))
    size = sizes[size_name]
    pdf(file=sprintf("%s\\%s\\frac_plots_no_legend%s.pdf",
                     write_path,
                     size_name,
                     ifelse(switchers, "_switchers", "")),
        height=size,
        width=size *2)
    plot(gf_no_legend)
    dev.off()  
    
    pdf(file=sprintf("%s\\%s\\frac_plots%s.pdf",
                     write_path,
                     size_name,
                     ifelse(switchers, "_switchers", "")),
        height=size,
        width=size *2)
    plot(gf)
    dev.off()  
  }
}


supp_figure_naive_learning_npcs_across_sessions <- function()
{
  write_path <- sprintf("%s\\supp_figure_naive_learning\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\supp_figure_naive_learning\\npcs_across_sessions\\",figures_path)
  dir.create(write_path)
  
  directions_to_use = c("Both", "Left", "Right")
  sizes=c(big=7,
          medium=6,
          medium_1=5,
          medium_2=4,
          small=3)
  
  surv_sizes=c(big=3,
               medium=2.75,
               medium_1=2.5,
               medium_2=2.25,
               medium_3=2,
               small=1.75,
               small_1=1.5)
  
  df_res <- c()
  df_non_pcs_sessions <- c()
  path_list <- all_data_paths
  old=F
  
  for (path_idx in 1:len(mice_ids)) {
    if (T) {
      
      first_path <- all_data_paths[(path_idx * 2 - 1)]
      second_path <- all_data_paths[(path_idx * 2)]
      
      path <- first_path
      if (grepl("CA3", path)) {
        mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', path))
        mice_str <- substr(path, mice_str_index, mice_str_index+4) 
      } else {
        mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', path))
        mice_str <- substr(path, mice_str_index, mice_str_index+3) 
      }
      
      mat_path <- sprintf("%s\\%s", first_path, "registered_cells_categories.Rda")
      load(mat_path)
      cells_mat_dir_1 <- cells_mat
      mt_dir_1 <- t(apply(cells_mat_dir_1, 1, function (cell) {as.numeric(cell != 3 & cell != 4)}))
      
      mat_path <- sprintf("%s\\%s", second_path, "registered_cells_categories.Rda")
      load(mat_path)
      cells_mat_dir_2 <- cells_mat
      mt_dir_2 <- t(apply(cells_mat_dir_2, 1, function (cell) {as.numeric(cell != 3 & cell != 4)}))
      
      
      
      
      mt_both_dir <- t(apply((mt_dir_1 & mt_dir_2), 1, as.numeric))
      
      suff_active_env_a <- which(apply(cbind(cells_mat_dir_1[,1:8],
                                             cells_mat_dir_2[,1:8]), 1, max) >= 1)
      
      suff_active_env_b <- which(apply(cbind(cells_mat_dir_1[,9:16],
                                             cells_mat_dir_2[,9:16]), 1, max) >= 1)
      
      
      mt_both_env_a <- mt_both_dir[suff_active_env_a, 1:8]
      mt_both_env_b <- mt_both_dir[suff_active_env_b, 9:16]
      
      both_alphabet <- c()
      both_sumord <- c()
      
      envs <- list(A=mt_both_env_a, B=mt_both_env_b)
      for (env_name in names(envs)) {
        env = envs[[env_name]]
        
        rownames(env) <- 1:nrow(env)
        ordered_mt <- env[order(env[,1]),]
        work_omt_1 <- ordered_mt[which(ordered_mt[,1] == 0),]
        
        alphabet_ord_1 <- order(apply(work_omt_1[,-1], 1, function(c) {paste(as.character(c), collapse = "")}))
        sumord_1 <- order(apply(work_omt_1[,-1], 1, sum))
        colnames(work_omt_1) <- c(1:8)
        
        work_omt_2 <- ordered_mt[which(ordered_mt[,1] == 1),]
        
        alphabet_ord_2 <- order(apply(work_omt_2[,-1], 1, function(c) {paste(as.character(c), collapse = "")}))
        sumord_2 <- order(apply(work_omt_2[,-1], 1, sum))
        colnames(work_omt_2) <- c(1:8)
        
        work_omt <- rbind(work_omt_1, work_omt_2)
        sumord <- c(sumord_1, sumord_2 + len(sumord_1))
        alphabet_ord <- c(alphabet_ord_1, alphabet_ord_2 + len(alphabet_ord_1))
        
        
        
        ph_sumord <- pheatmap(work_omt[sumord,], cluster_rows=F, cluster_cols=F,
                              col=c("white", "black"),
                              legend = F,
                              main=mice_str)
        
        ph_alphord <- pheatmap(work_omt[alphabet_ord,], cluster_rows=F, cluster_cols=F,
                               col=c("white", "black"),
                               legend = F,
                               main=mice_str)
        
        
        ph_sumord_pc <- pheatmap(work_omt_1[sumord_1,], cluster_rows=F, cluster_cols=F,
                                 col=c("white", "black"),
                                 legend = F,
                                 main=mice_str)
        
        ph_alphord_pc <- pheatmap(work_omt_1[alphabet_ord_1,], cluster_rows=F, cluster_cols=F,
                                  col=c("white", "black"),
                                  legend = F,
                                  main=mice_str)
        
        
        ph_sumord_npc <- pheatmap(work_omt_2[sumord_2,], cluster_rows=F, cluster_cols=F,
                                  col=c("white", "black"),
                                  legend = F,
                                  main=mice_str)
        
        ph_alphord_npc <- pheatmap(work_omt_2[alphabet_ord_2,], cluster_rows=F, cluster_cols=F,
                                   col=c("white", "black"),
                                   legend = F,
                                   main=mice_str)
        
        for (size_name in names(surv_sizes)) {
          size = surv_sizes[[size_name]]
          
          dir.create(sprintf("%s\\%s",
                             write_path,
                             size_name))
          
          dir.create(sprintf("%s\\%s\\survival\\",
                             write_path,
                             size_name))
          
          dir.create(sprintf("%s\\%s\\survival\\%s",
                             write_path,
                             size_name,
                             mice_str))
          
          
          pdf(sprintf("%s\\%s\\survival\\%s\\all_sum_env_%s.pdf",
                      write_path,
                      size_name,
                      mice_str,
                      env_name),
              height=size,
              width=size)
          plot(ph_sumord[[4]])
          dev.off()
          
          pdf(sprintf("%s\\%s\\survival\\%s\\all_alph_env_%s.pdf",
                      write_path,
                      size_name,
                      mice_str,
                      env_name),
              height=size,
              width=size)
          plot(ph_alphord[[4]])
          dev.off()
          
          pdf(sprintf("%s\\%s\\survival\\%s\\npc_sum_env_%s.pdf",
                      write_path,
                      size_name,
                      mice_str,
                      env_name),
              height=size,
              width=size)
          plot(ph_sumord_npc[[4]])
          dev.off()
          
          pdf(sprintf("%s\\%s\\survival\\%s\\npc_alph_env_%s.pdf",
                      write_path,
                      size_name,
                      mice_str,
                      env_name),
              height=size,
              width=size)
          plot(ph_alphord_npc[[4]])
          dev.off()        
          
          pdf(sprintf("%s\\%s\\survival\\%s\\pc_sum_env_%s.pdf",
                      write_path,
                      size_name,
                      mice_str,
                      env_name),
              height=size,
              width=size)
          plot(ph_sumord_pc[[4]])
          dev.off()
          
          pdf(sprintf("%s\\%s\\survival\\%s\\pc_alph_%s.pdf",
                      write_path,
                      size_name,
                      mice_str,
                      env_name),
              height=size,
              width=size)
          plot(ph_alphord_pc[[4]])
          dev.off()        
        }
        print(1-(sum(rowSums(work_omt) == 32) / nrow(work_omt)))
      }
    }
  }
  
  for (path_idx in 2:len(all_data_paths)) {
    
    path <- all_data_paths[[path_idx]]
    
    mat_path <- sprintf("%s\\%s", path, "registered_cells_categories.Rda")
    load(mat_path)
    
    
    
    most_advent_status <- unlist(lapply(apply(cells_mat, 1, table), function(tb) {max(as.numeric(names(tb)))}))
    non_pcs <- which(most_advent_status == 1)
    
    
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
    registration_path <- sprintf("%s%s",
                                 str_split(path, "Matrices\\\\")[[1]][1], registration_mat_filename)
    
    registration_mat <- readMat(registration_path)
    registration_mat <- registration_mat$cell.to.index.map
    
    if (grepl("CA3", path)) {
      mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', path))
      mice_str <- substr(path, mice_str_index, mice_str_index+4) 
    } else {
      mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', path))
      mice_str <- substr(path, mice_str_index, mice_str_index+3) 
    }
    
    
    direc <- ifelse(grepl("Left", path), "L", "R")
    
    for (npc_idx in non_pcs) {
      cells_indices_vec <- registration_mat[npc_idx,]
      
      cell_across_sessions <- list()
      for (session in 1:len(cells_indices_vec)) {
        
        ind <- cells_indices_vec[session]
        
        st <- t(spike_train[[session]][[1]])
        stim <- stim <- as.vector(stim_trace[[session]][[1]])
        
        
        
        tmp <- equalize_prior(st, stim, verbose=T)
        
        
        st <- tmp$spike_train
        stim <- tmp$stim_trace
        
        run_df <- data.frame(Time=((1:len(stim) / len(stim)) * 20), 
                             Position=stim * 4) 
        run_df$Position[which(run_df$Position == 4)] <- 0
        firing_ind <- which(st[ind,] > 0)
        
        
        grundf <-
          ggplot(run_df, aes(x=Time, y=Position)) + 
          #geom_line(col="gray30") + 
          theme_light() +  
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                legend.position="NA",
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks=element_line(color="black"),
                axis.text=element_text(color="black"),
                panel.border = element_rect(color="black"),
                panel.background = element_blank(),
                plot.title=element_text(size=7)) +
          xlim(20,0) +#plot.margin = margin(t=0,r=0,b=0,l=0,unit="cm")) +
          ylim(0,96) +
          ylab("Position (cm)") +
          xlab("Time (minutes)") +  
          coord_flip() +
          ggtitle(sprintf("Session %d", session))
        # geom_vline(xintercept = polygon_df$x[1],
        #            linetype="dashed",
        #            color="black") 
        
        if (ind != 0 ) {
          
          cell_df <- data.frame(Positions=stim[firing_ind] * 4,
                                Times=run_df$Time[firing_ind])
          
          grundf <- grundf +
            geom_point(data=cell_df, aes(x=Times,y=Positions),
                       fill="#DA1C5C",
                       color="#DA1C5C",
                       stroke=0,
                       size=2,
                       alpha=.8)
        }
        
        
        
        
        cell_across_sessions <- append(cell_across_sessions,
                                       list(grundf))
      }
      
      cell_across_sessions$nrow <- 4
      across_plot <- do.call(plot_grid, cell_across_sessions)
      
      for (size_name in names(sizes)) {
        size = sizes[[size_name]]
        
        dir.create(sprintf("%s\\%s",
                           write_path,
                           size_name))
        
        dir.create(sprintf("%s\\%s\\%s_%s",
                           write_path,
                           size_name,
                           mice_str,
                           direc))
        
        
        pdf(sprintf("%s\\%s\\%s_%s\\npc_across_%d.pdf",
                    write_path,
                    size_name,
                    mice_str,
                    direc,
                    npc_idx),
            height=size,
            width=size * .6)
        plot(across_plot)
        dev.off()
        
      }
      
    }
  }
}


supp_figure_naive_learning_npcs_across_sessions <- function()
{
  write_path <- sprintf("%s\\supp_figure_naive_learning\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\supp_figure_naive_learning\\first_day_pcs_across_sessions\\",figures_path)
  dir.create(write_path)
  
  directions_to_use = c("Both", "Left", "Right")
  sizes=c(big=3.5,
          medium=3,
          medium_1=2.5,
          medium_2=2,
          small=1.5)
  
  surv_sizes=c(big=3,
               medium=2.75,
               medium_1=2.5,
               medium_2=2.25,
               medium_3=2,
               small=1.75,
               small_1=1.5)
  
  df_res <- c()
  df_non_pcs_sessions <- c()
  path_list <- all_data_paths
  old=F
  
  
  for (path_idx in 1:len(all_data_paths)) {
    
    path <- all_data_paths[[path_idx]]
    
    mat_path <- sprintf("%s\\%s", path, "registered_cells_categories.Rda")
    load(mat_path)
    
    
    
    
    
    
    
    non_pcs <- which(cells_mat[,1] == 4)# | cells_mat[,1] == 1 | cells_mat[,1] == 2 | cells_mat[,1] == 3)
    #non_pcs <- which(cells_mat[,1] == 2 | cells_mat[,1] == 1 | cells_mat[,1] == 2 | cells_mat[,1] == 3)
    
    
    
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
    registration_path <- sprintf("%s%s",
                                 str_split(path, "Matrices\\\\")[[1]][1], registration_mat_filename)
    
    registration_mat <- readMat(registration_path)
    registration_mat <- registration_mat$cell.to.index.map
    
    if (grepl("CA3", path)) {
      mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', path))
      mice_str <- substr(path, mice_str_index, mice_str_index+4) 
    } else {
      mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', path))
      mice_str <- substr(path, mice_str_index, mice_str_index+3) 
    }
    
    
    direc <- ifelse(grepl("Left", path), "L", "R")
    
    for (npc_idx in non_pcs) {
      cells_indices_vec <- registration_mat[npc_idx,]
      
      cell_across_sessions <- list()
      for (session in (1:len(cells_indices_vec))[1:8]) {
        
        ind <- cells_indices_vec[session]
        
        st <- t(spike_train[[session]][[1]])
        stim <- stim <- as.vector(stim_trace[[session]][[1]])
        
        
        
        tmp <- equalize_prior(st, stim, verbose=T)
        
        
        st <- tmp$spike_train
        stim <- tmp$stim_trace
        
        run_df <- data.frame(Time=((1:len(stim) / len(stim)) * 20), 
                             Position=stim * 4) 
        run_df$Position[which(run_df$Position == 4)] <- 0
        firing_ind <- which(st[ind,] > 0)
        
        
        grundf <-
          ggplot(run_df, aes(x=Time, y=Position)) + 
          #geom_line(col="gray30") + 
          theme_light() +  
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                legend.position="NA",
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks=element_line(color="black"),
                axis.text=element_text(color="black"),
                panel.border = element_rect(color="black"),
                panel.background = element_blank(),
                plot.title=element_text(size=7)) +
          xlim(20,0) +#plot.margin = margin(t=0,r=0,b=0,l=0,unit="cm")) +
          ylim(0,96) +
          ylab("Position (cm)") +
          xlab("Time (minutes)") +  
          coord_flip() +
          ggtitle(sprintf("Session %d", session))
        # geom_vline(xintercept = polygon_df$x[1],
        #            linetype="dashed",
        #            color="black") 
        
        if (ind != 0 ) {
          
          cell_df <- data.frame(Positions=stim[firing_ind] * 4,
                                Times=run_df$Time[firing_ind])
          
          grundf <- grundf +
            geom_point(data=cell_df, aes(x=Times,y=Positions),
                       fill="#DA1C5C",
                       color="#DA1C5C",
                       stroke=0,
                       size=2,
                       alpha=.8)
        }
        
        
        
        
        cell_across_sessions <- append(cell_across_sessions,
                                       list(grundf))
      }
      
      cell_across_sessions$nrow <- 2
      across_plot <- do.call(plot_grid, cell_across_sessions)
      
      for (size_name in names(sizes)) {
        size = sizes[[size_name]]
        
        dir.create(sprintf("%s\\%s",
                           write_path,
                           size_name))
        
        dir.create(sprintf("%s\\%s\\%s_%s",
                           write_path,
                           size_name,
                           mice_str,
                           direc))
        
        
        prfx = ifelse(npc_idx %in% which(cells_mat[,1] == 3), "SWITCHER_", "")
        
        png(sprintf("%s\\%s\\%s_%s\\%snpc_across_%d.png",
                    write_path,
                    size_name,
                    mice_str,
                    direc,
                    prfx,
                    npc_idx),
            height=size,
            width=size * 1.2,
            units="in",
            res=300)
        plot(across_plot)
        dev.off()
        
      }
      
    }
  }
}
