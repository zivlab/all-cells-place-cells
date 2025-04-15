
library(cowplot)

wilcox_sample_test <- function(a,b) {res <- list(); res$p.value <- mean(unlist(sapply(1:5, function(idx) {wilcox.test(sample(a, 30), sample(b,30))}$p.value))); res$statistic <- mean(unlist(sapply(1:5, function(idx) {wilcox.test(sample(a, 30), sample(b,30))}$statistic))); return(res)}

permutation_test <- function(a,b, central_tend_f=mean, plot=F) {
  sample_diff <- central_tend_f(a) - central_tend_f(b)
  
  pooled <- c(a,b)    
  all_perms  <- c()
  for (i in 1:1000) {
    s1 <-  sample(pooled, len(a))
    s2 <-  sample(pooled, len(b))
    all_perms <- c(all_perms,  central_tend_f(s1) - central_tend_f(s2))
  }
  
  if (plot) {
    hist(all_perms, xlab="Diff", main="");
    abline(v=sample_diff, col="red", lty=2, lwd=2);
  }
  
  cumulative_func <- ecdf(all_perms)
  return(list(statistic=sample_diff, p.value=cumulative_func(sample_diff)))
}


get_counts_vec <- function(a,b,nbreaks) {
  vec_of_breaks = seq(min(a,b), max(a,b), length.out=nbreaks) 
  h1 <- hist(a, breaks=vec_of_breaks, plot=F)$counts
  h2 <- hist(b, breaks=vec_of_breaks, plot=F)$counts
  
  return(list(h1, h2))
  
}

chisq_goodness_of_fit <- function(a, b, nbreaks, plot=F) {
  
  
  counts_list <- get_counts_vec(a,b,nbreaks)
  
  
  if(plot) {
    barplot(counts_list[[1]], col="red")
    barplot(counts_list[[2]], col=adjustcolor("blue", alpha=0.3), add=T)
  }
  
  prob_vec <- counts_list[[2]] / sum(counts_list[[2]])
  return(chisq.test(counts_list[[1]], p=prob_vec + 10^-30))
}

#figures_path <- "~/..\\Desktop\\allcells_figures"
figures_path= "/Users/itayta/Downloads/yaniv_place_cells_project/allcells_figures"

# figure_1 <- function() {
#   
#   write_path <- sprintf("%s\\figure_1\\",figures_path)
#   dir.create(write_path)
#   
#   a <- get_spike_train_and_stim_trace_from_path(data_path_3, 2)
#   cells_to_plot <- c(430, 355, 170, 189, 10, 26, 75, 135, 405, 738, 459, 218, 227, 80, 136, 146, 177, 184, 194, 227, 256, 281)
#   
#   b <- get_spike_train_and_stim_trace_from_path(test_path3, 2, simulate = T)
#   ext = "c6m3_R2_simulation_"
#   
#   
#   a <- get_spike_train_and_stim_trace_from_path(test_path2, 5)
#   figure_1_top_row(stim_trace=a[[2]], spike_train = a[[1]], ext = "c5m1_L5_", all_cells = T)
#   
#   
#   ca1_paths <- sprintf("%s\\%s", all_data_paths[1:8], "equalized")
#   ca3_paths <- sprintf("%s\\%s", all_data_paths[9:18], "equalized")
#   #by_activity <- plot_place_cell_pct_by_activity_plot(ca1_paths, cyclic_only=T, )
#   by_duration_ca1 <- plot_place_cell_pct_by_sample_duration_figure(ca1_paths, bin_size = 2.5, 
#                                                                    ext = "", 
#                                                                    cyclic_only=T, 
#                                                                    fit_sim = "", 
#                                                                    verbose=T, 
#                                                                    ses_ind = 1:16,
#                                                                    bwidth = 1)
#   
#   by_duration_ca3 <- plot_place_cell_pct_by_sample_duration_figure(ca3_paths, bin_size = 2.5, 
#                                                                    ext = "", 
#                                                                    cyclic_only=T, 
#                                                                    fit_sim = "", 
#                                                                    verbose=T, 
#                                                                    ses_ind = 1:16)
#   by_pval_ca1 <- plot_pval_by_sample_duration_figure(ca1_paths, 
#                                                      ext = "", 
#                                                      fit_sim = "", 
#                                                      verbose=F, 
#                                                      type1 = F, 
#                                                      ses_ind=1:16,
#                                                      bwidth = 1)
#   
#   by_pval_ca3 <- plot_pval_by_sample_duration_figure(ca3_paths, 
#                                                      ext = "", 
#                                                      fit_sim = "", 
#                                                      verbose=F, 
#                                                      type1 = F, 
#                                                      ses_ind=1:16)
#   
#   
#   overlaid_pval <- by_pval_ca1[[5]] + 
#                    geom_line(data=by_pval_ca3[[2]], 
#                              aes(x=time, y=pvd), 
#                              size=1, 
#                              color="black") + 
#                     geom_errorbar(size=1, 
#                                   data=by_pval_ca3[[2]], 
#                                   aes(ymin=pvd-pvdsd, ymax=pvd+pvdsd), 
#                                   color=adjustcolor("black", alpha=1), width=1)
#   overlaid_duration <- by_duration_ca1[[6]] + 
#                        geom_line(data=by_duration_ca3[[4]], 
#                                  aes(x=time, y=fraction), 
#                                  size=1, color="black") + 
#                        geom_errorbar(size=1, 
#                                      data=by_duration_ca3[[4]], 
#                                      aes(ymin=fraction-sd, ymax=fraction+sd), 
#                                      color=adjustcolor("black", alpha=1), width=1)
#   
#   
#   
#   overlaid_pval_ribbon <- 
#     by_pval_ca1[[1]] + 
#     geom_line(data=by_pval_ca3[[2]], 
#               aes(x=time, y=pvd), 
#               size=1, 
#               color="black") + 
#     geom_ribbon(size=1, 
#                 data=by_pval_ca3[[2]], 
#                 aes(ymin=pvd-pvdsd, ymax=pvd+pvdsd), 
#                 fill=adjustcolor("black", alpha=0.2), 
#                 width=1,
#                 color=NA)
#   overlaid_duration_ribbon <- 
#     by_duration_ca1[[1]] + 
#     geom_line(data=by_duration_ca3[[4]], 
#               aes(x=time, y=fraction), 
#               size=1, color="black") + 
#     geom_ribbon(size=1, 
#                 data=by_duration_ca3[[4]], 
#                 aes(ymin=fraction-sd, ymax=fraction+sd), 
#                 fill=adjustcolor("black", alpha=0.2), 
#                 width=1,
#                 color=NA)
#   
#   w_error_bar_p <- arrangeGrob(overlaid_duration, overlaid_pval, nrow=1)
#   w_ribbon_p <- arrangeGrob(overlaid_duration_ribbon, overlaid_pval_ribbon, nrow=1)
#   
#   pdf(sprintf("%s\\bottom_row_errors.pdf",write_path), height=3, width=6); 
#   plot(w_error_bar_p); dev.off()   
#   
#   pdf(sprintf("%s\\bottom_row_ribbons.pdf",write_path), height=2.5, width=5); 
#   plot(w_ribbon_p); dev.off()   
# }
# 
# 
# 
# figure_2 <- function() {
#   write_path <- sprintf("%s\\figure_2\\",figures_path)
#   dir.create(write_path)  
#   
#   figure_2_top_row()
#   figure_2_middle_row()
#   figure_2_bottom_row()
#   #pdf(sprintf("%s\\simulation_hists_%s.pdf",write_path, ext), height=1.9 * 4/3, width=1.9 * 4/3); plot(hists_p); dev.off() 
#   
# }
# 
# 
# figure_3 <- function(paths,  estimated_df_folder="KS_simulations_likelihood_c")  {
#   
#   write_path <- sprintf("%s\\figure_3\\",figures_path)
#   dir.create(write_path)
#   
#   true_data_paths <- sprintf("%s\\%s", all_data_paths[1:8], "equalized")
#   
#   sessions_to_use=c(5:8,13:16)
#   
#   all_df <- get_all_df(true_data_paths, sessions_to_use)
#   
#   percent_active_df <- all_df[,c(7:9)]
#   percent_place_cells_df <- all_df[,c(1:3)]
#   percent_all_df <- all_df[,c(4:6)]
#   
#   colnames(percent_active_df) <- c("Both", "Right", "Left")
#   colnames(percent_place_cells_df) <- c("Both", "Right", "Left")
#   colnames(percent_all_df) <- c("Both", "Right", "Left")
#   
#   melted_percent_active <- melt(percent_active_df)
#   melted_percent_place_cells <- melt(percent_place_cells_df)
#   melted_percent_all <- melt(percent_all_df)
#   
#   ylabs = c("Percentage of active cells (%)",
#             "Percentage of place cells (% of active cells)",
#             "Percentage of place cells (% of all cells)")
#   
#   dfs <- list(melted_percent_active,
#               melted_percent_place_cells,
#               melted_percent_all)
#   
#   bar_plots_list <- list()
#   for (df_idx in 1:3) {
#     df_to_use <- dfs[[df_idx]]
#     colnames(df_to_use) <- c("Direction", "Percent")
#     sd_df <- ddply(df_to_use, .(Direction), function(sub_df) {return(c(mean(sub_df[,"Percent"]), sd(sub_df[,"Percent"])))})
#     colnames(sd_df) <- c("Direction", "Percent", "Sd")
#     gr <- 
#       ggplot(sd_df, aes(x=Direction, y=Percent)) + 
#       geom_jitter(data=df_to_use,
#                   aes(x=Direction, y=Percent, group=Direction), position=position_jitter(0.2), color="gray20", size=1, alpha=0.1) + 
#       geom_boxplot(data=df_to_use, aes(group=Direction), width=0.5, fill=adjustcolor("gray80", alpha=0.8), color="gray65", size=0.75) + 
#       theme_light() +     
#       theme(panel.grid.major = element_blank(), 
#             panel.grid.minor = element_blank(),
#             axis.line = element_line(colour = "black"),
#             panel.border = element_blank(),
#             panel.background = element_blank()) + 
#       xlab("Running direction") +
#       ylab(ylabs[df_idx]) + 
#       ylim(0,1)
#     
#     bar_plots_list <- append(bar_plots_list,
#                              list(gr))
#   }
#   
#   bar_plots_list$nrow <- 1
#   first_row <- do.call(grid.arrange, bar_plots_list)
#   
#   estimated_res <- get_estimated_df(all_df, true_data_paths, sessions_to_use)
#   
#   estimation_list <- estimated_res$estimation_list
#   likelihood_list <- estimated_res$likelihood_list
#   measured_percent_list <- estimated_res$measured_percent_list
#   estimated_per_percent_list <- estimated_res$estimated_per_percent_list
#   session_metavar_list <- estimated_res$session_metavar_list
#   missing_paths <- estimated_res$session_metavar_list 
#   likelihood_df <- estimated_res$likelihood_df
#   estimated_per_percent_df <- estimated_res$estimated_per_percent_df
#   session_df <- estimated_res$session_df
#   estimation_vec <- estimated_res$estimation_vec
#   measured_vec <- estimated_res$measured_vec
#   estimated_df <- estimated_res$estimated_df
#   
#   
#   both_dir_res <- get_both_dir_df(estimated_df = estimated_df, all_df=all_df)
#   
#   both_dir_df <- both_dir_res$both_dir_df
#   mean_sd_df <- both_dir_res$mean_sd_df
#   session_mean_sd_df <- both_dir_res$session_mean_sd_df
#   
#   gest <- 
#     ggplot(mean_sd_df, aes(x=Direction, y=Estimated)) + 
#     geom_boxplot(data=both_dir_df, 
#                  aes(group=Direction), width=0.5, fill=adjustcolor("gray80", alpha=0.8), color="gray65", size=1) + 
#     geom_jitter(data=both_dir_df,
#                 aes(x=Direction, y=Estimated, group=Direction), position=position_jitter(0.2), color="gray20", size=0.75, alpha=0.1) + 
#     theme_light() +     
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           axis.line = element_line(colour = "black"),
#           panel.border = element_blank(),
#           panel.background = element_blank()) + 
#     xlab("Running direction") +
#     ylab(ylabs[2]) + 
#     ylim(0,1.05)
#   
#   
#   likelihood_mean_sd_df <- ddply(melted_likelihood, .(`Simulated percent`), 
#                                  function(sub_df) {
#                                    return(c(mean(sub_df[,"Normalized likelihood"]),
#                                             sd(sub_df[,"Normalized likelihood"])))
#                                  })
#   
#   colnames(likelihood_mean_sd_df) <- c("Simulated percent", "Mean", "Sd")
#   
#   melted_likelihood$`Simulated percent` <- factor(melted_likelihood$`Simulated percent` * 100,
#                                                   levels=as.character(likelihood_mean_sd_df$`Simulated percent`  * 100))
#   likelihood_mean_sd_df$`Simulated percent` <- factor(likelihood_mean_sd_df$`Simulated percent` * 100, 
#                                                       levels=as.character(likelihood_mean_sd_df$`Simulated percent` * 100))
#   
#   glikelihood <- 
#     ggplot(likelihood_mean_sd_df, aes(x=`Simulated percent`, y=Mean)) + 
#     geom_bar(stat="summary", aes(group=`Simulated percent`), 
#              width=0.5, fill="gray50", color="black", position = "dodge") + 
#     geom_errorbar(aes(ymin=Mean, ymax=Mean + Sd, group=`Simulated percent`),
#                   size=2, width=0.3) +
#     # geom_jitter(data=melted_likelihood,
#     #             aes(x=`Simulated percent`, y=`Normalized likelihood`, group=`Simulated percent`), 
#     #             position=position_jitter(0.05), color="gray40", size=3, alpha=0.5) + 
#     theme_light() +     
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           axis.line = element_line(colour = "black"),
#           panel.border = element_blank(),
#           panel.background = element_blank()) + 
#     xlab("Simulated percent (%)") +
#     ylab("Normalized likelihood") + 
#     ylim(0,1.05)
#   
#   
#   scatter_list <- list()
#   for (est_p in as.character(seq(0.5, 1, by=0.1))) {
#     
#     tmp_df <- data.frame(Measured=estimated_per_percent_df[,"measured_vec"],
#                          Estimated=estimated_per_percent_df[,est_p],
#                          Session=as.numeric(estimated_per_percent_df[,ncol(estimated_per_percent_df) - 3]))
#     
#     tmp_df$Session[tmp_df$Session > 8] <- tmp_df$Session[tmp_df$Session > 8] - 8
#     tmp_df$Session <- factor(sprintf("%dth", tmp_df$Session), c("5th", "6th", "7th", "8th"))
#     
#     linear_line_df <- data.frame(x=c(0.3,1),
#                                  y=c(0.3,1))
#     gs <- 
#       ggplot(linear_line_df) +
#       geom_line(aes(x=x,y=y), linetype="longdash", color="#c0c1c3",
#                 size=1.5, alpha=0.5) +    
#       geom_point(data=tmp_df, aes(x=Measured, y=Estimated, color=Session), size=1.5) +
#       theme_light() +     
#       theme(panel.grid.major = element_blank(), 
#             panel.grid.minor = element_blank(),
#             axis.line = element_line(colour = "black"),
#             panel.border = element_blank(),
#             panel.background = element_blank(),
#             legend.position="NA",
#             plot.margin=unit(c(0,0,0,0), "cm")) +
#       scale_color_manual(breaks=c("5th","6th","7th","8th"),
#                          values=c("#3f9ab4",
#                                   "#6cbfa6",
#                                   "#b0c7e0",
#                                   "#8a83ad")) + 
#       ylim(c(0.3,1)) + 
#       xlim(c(0.3,1)) + 
#       ylab("") +
#       xlab("")
#     print(sprintf("MSE = %f",
#                   mean((tmp_df$Estimated - tmp_df$Measured) ** 2)))
#     
#     scatter_list <- append(scatter_list, list(gs))
#   }
#   
#   scatter_list$nrow = 2
#   #scatter_list$left="Estimated place cell percentage (%)"
#   #scatter_list$bottom="Measured place cell percentage (%)"
#   scattersp <- do.call(grid.arrange, scatter_list)
#   plot_h <- 1.9 * 4 / 3
#   
#   second_row <- grid.arrange(gest, glikelihood, align_plots(scattersp)[[1]], 
#                              nrow=1,widths=c(plot_h,plot_h,plot_h * 1.5))
#   aligned <- align_plots(first_row, second_row, align="v")
#   
#   ext =""  
#   
#   
#   
#   pdf(file=sprintf("%s\\first_row%s.pdf", write_path,  ext),
#       height=plot_h,
#       width=plot_h * 3)
#   plot(first_row)
#   dev.off() 
#   
#   pdf(file=sprintf("%s\\second_row%s.pdf", write_path, ext),
#       height=plot_h,
#       width=plot_h * 3.5)
#   
#   plot(aligned[[2]])
#   dev.off() 
# }
# 
# figure_3_ca3 <- function(paths)  {
#   write_path <- sprintf("%s\\figure_3_ca3\\",figures_path)
#   dir.create(write_path)
#   
#   true_data_paths <- sprintf("%s\\%s", paths, "equalized")
#   
#   sessions_to_use=c(5:8,13:16)
#   
#   all_df <- get_all_df(true_data_paths, sessions_to_use)
#   
#   percent_active_df <- all_df[,c(7:9)]
#   percent_place_cells_df <- all_df[,c(1:3)]
#   percent_all_df <- all_df[,c(4:6)]
#   
#   colnames(percent_active_df) <- c("Both", "Right", "Left")
#   colnames(percent_place_cells_df) <- c("Both", "Right", "Left")
#   colnames(percent_all_df) <- c("Both", "Right", "Left")
#   
#   melted_percent_active <- melt(percent_active_df)
#   melted_percent_place_cells <- melt(percent_place_cells_df)
#   melted_percent_all <- melt(percent_all_df)
#   
#   ylabs = c("Percentage of active cells (%)",
#             "Percentage of place cells (% of active cells)",
#             "Percentage of place cells (% of all cells)")
#   
#   dfs <- list(melted_percent_active,
#               melted_percent_place_cells,
#               melted_percent_all)
#   
#   bar_plots_list <- list()
#   for (df_idx in 1:3) {
#     df_to_use <- dfs[[df_idx]]
#     colnames(df_to_use) <- c("Direction", "Percent")
#     sd_df <- ddply(df_to_use, .(Direction), function(sub_df) {return(c(mean(sub_df[,"Percent"]), sd(sub_df[,"Percent"])))})
#     colnames(sd_df) <- c("Direction", "Percent", "Sd")
#     gr <- 
#       ggplot(sd_df, aes(x=Direction, y=Percent)) + 
#       geom_jitter(data=df_to_use,
#                   aes(x=Direction, y=Percent, group=Direction), position=position_jitter(0.2), color="gray20", size=1, alpha=0.1) + 
#       geom_boxplot(data=df_to_use, aes(group=Direction), width=0.5, fill=adjustcolor("royalblue2", alpha=0.8), color="royalblue4", size=0.75) + 
#       theme_light() +     
#       theme(panel.grid.major = element_blank(), 
#             panel.grid.minor = element_blank(),
#             axis.line = element_line(colour = "black"),
#             panel.border = element_blank(),
#             panel.background = element_blank()) + 
#       xlab("Running direction") +
#       ylab(ylabs[df_idx]) + 
#       ylim(0,1)
#     
#     bar_plots_list <- append(bar_plots_list,
#                              list(gr))
#   }
#   
#   bar_plots_list$nrow <- 1
#   first_row <- do.call(grid.arrange, bar_plots_list)
#   
#   estimated_res <- get_estimated_df(all_df, true_data_paths, sessions_to_use)
#   
#   estimation_list <- estimated_res$estimation_list
#   likelihood_list <- estimated_res$likelihood_list
#   measured_percent_list <- estimated_res$measured_percent_list
#   estimated_per_percent_list <- estimated_res$estimated_per_percent_list
#   session_metavar_list <- estimated_res$session_metavar_list
#   missing_paths <- estimated_res$session_metavar_list 
#   likelihood_df <- estimated_res$likelihood_df
#   estimated_per_percent_df <- estimated_res$estimated_per_percent_df
#   session_df <- estimated_res$session_df
#   estimation_vec <- estimated_res$estimation_vec
#   measured_vec <- estimated_res$measured_vec
#   estimated_df <- estimated_res$estimated_df
#   
#   both_dir_res <- get_both_dir_df(estimated_df = estimated_df, all_df=all_df)
#   
#   both_dir_df <- both_dir_res$both_dir_df
#   mean_sd_df <- both_dir_res$mean_sd_df
#   session_mean_sd_df <- both_dir_res$session_mean_sd_df
#   
#   gest <- 
#     ggplot(mean_sd_df, aes(x=Direction, y=Estimated)) + 
#     geom_jitter(data=both_dir_df,
#                 aes(x=Direction, y=Estimated, group=Direction), position=position_jitter(0.2), color="gray20", size=1, alpha=0.1) + 
#     geom_boxplot(data=both_dir_df, aes(group=Direction), width=0.5, fill=adjustcolor("royalblue2", alpha=0.8), color="royalblue4", size=0.75) + 
#     theme_light() +     
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           axis.line = element_line(colour = "black"),
#           panel.border = element_blank(),
#           panel.background = element_blank()) + 
#     xlab("Running direction") +
#     ylab(ylabs[2]) + 
#     ylim(0,1.05)
#   
#   
#   likdf_tmp <- likelihood_df
#   likelihood_df <- (t(apply(likelihood_df, 1, function(r) {return(r/sum(r))})))
#   likelihood_df <- as.data.frame(likelihood_df)
#   likelihood_sessions <- as.numeric(session_df[,1])
#   likelihood_sessions[likelihood_sessions > 8] <-  likelihood_sessions[likelihood_sessions > 8] - 8
#   likelihood_df <- cbind(likelihood_df, likelihood_sessions)
#   colnames(likelihood_df) <- c(as.character(seq(0.5,1,by=0.1)), "sessions")
#   melted_likelihood <- melt(likelihood_df, id.vars="sessions")
#   colnames(melted_likelihood) <- c("Session", "Simulated percent", "Normalized likelihood")
#   
#   likelihood_mean_sd_df <- ddply(melted_likelihood, .(`Simulated percent`), 
#                                  function(sub_df) {
#                                    return(c(mean(sub_df[,"Normalized likelihood"]),
#                                             sd(sub_df[,"Normalized likelihood"])))
#                                  })
#   
#   colnames(likelihood_mean_sd_df) <- c("Simulated percent", "Mean", "Sd")
#   
#   melted_likelihood$`Simulated percent` <- factor(melted_likelihood$`Simulated percent` * 100,
#                                                   levels=as.character(likelihood_mean_sd_df$`Simulated percent`  * 100))
#   likelihood_mean_sd_df$`Simulated percent` <- factor(likelihood_mean_sd_df$`Simulated percent` * 100, 
#                                                       levels=as.character(likelihood_mean_sd_df$`Simulated percent` * 100))
#   
#   glikelihood <- 
#     ggplot(likelihood_mean_sd_df, aes(x=`Simulated percent`, y=Mean)) + 
#     geom_bar(stat="summary", aes(group=`Simulated percent`), 
#              width=0.5, fill="gray50", color="black", position = "dodge") + 
#     geom_errorbar(aes(ymin=Mean, ymax=Mean + Sd, group=`Simulated percent`),
#                   size=2, width=0.3) +
#     # geom_jitter(data=melted_likelihood,
#     #             aes(x=`Simulated percent`, y=`Normalized likelihood`, group=`Simulated percent`), 
#     #             position=position_jitter(0.05), color="gray40", size=3, alpha=0.5) + 
#     theme_light() +     
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           axis.line = element_line(colour = "black"),
#           panel.border = element_blank(),
#           panel.background = element_blank()) + 
#     xlab("Simulated percent (%)") +
#     ylab("Normalized likelihood") + 
#     ylim(0,1.05)
#   
#   
#   scatter_list <- list()
#   for (est_p in as.character(seq(0.5, 1, by=0.1))) {
#     
#     tmp_df <- data.frame(Measured=estimated_per_percent_df[,"measured_vec"],
#                          Estimated=estimated_per_percent_df[,est_p],
#                          Session=as.numeric(estimated_per_percent_df[,ncol(estimated_per_percent_df) - 3]))
#     
#     tmp_df$Session[tmp_df$Session > 8] <- tmp_df$Session[tmp_df$Session > 8] - 8
#     tmp_df$Session <- factor(sprintf("%dth", tmp_df$Session), c("5th", "6th", "7th", "8th"))
#     
#     linear_line_df <- data.frame(x=c(0.3,1),
#                                  y=c(0.3,1))
#     gs <- 
#       ggplot(linear_line_df) +
#       geom_line(aes(x=x,y=y), linetype="longdash", color="#c0c1c3",
#                 size=1.5, alpha=0.5) +    
#       geom_point(data=tmp_df, aes(x=Measured, y=Estimated, color=Session), size=1.5) +
#       theme_light() +     
#       theme(panel.grid.major = element_blank(), 
#             panel.grid.minor = element_blank(),
#             axis.line = element_line(colour = "black"),
#             panel.border = element_blank(),
#             panel.background = element_blank(),
#             legend.position="NA",
#             plot.margin=unit(c(0,0,0,0), "cm")) +
#       scale_color_manual(breaks=c("5th","6th","7th","8th"),
#                          values=c("#3f9ab4",
#                                   "#6cbfa6",
#                                   "#b0c7e0",
#                                   "#8a83ad")) + 
#       ylim(c(0.3,1)) + 
#       xlim(c(0.3,1)) + 
#       ylab("") +
#       xlab("")
#     print(sprintf("MSE = %f",
#                   mean((tmp_df$Estimated - tmp_df$Measured) ** 2)))
#     
#     scatter_list <- append(scatter_list, list(gs))
#   }
#   
#   scatter_list$nrow = 2
#   #scatter_list$left="Estimated place cell percentage (%)"
#   #scatter_list$bottom="Measured place cell percentage (%)"
#   scattersp <- do.call(grid.arrange, scatter_list)
#   plot_h <- 1.9 * 4 / 3
#   
#   second_row <- grid.arrange(gest, glikelihood, align_plots(scattersp)[[1]], 
#                              nrow=1,widths=c(plot_h,plot_h,plot_h * 1.5))
#   aligned <- align_plots(first_row, second_row, align="v")
#   
#   
#   
#   
#   pdf(file=sprintf("%s\\first_row%s.pdf", write_path,  ext),
#       height=plot_h,
#       width=plot_h * 3)
#   plot(first_row)
#   dev.off() 
#   
#   pdf(file=sprintf("%s\\second_row%s.pdf", write_path, ext),
#       height=plot_h,
#       width=plot_h * 3.5)
#   
#   plot(aligned[[2]])
#   dev.off() 
# }
# 
# figure_4_extended <- function(paths, estimated_df_folder="KS_simulations_likelihood_c")  {
#   
#   paths <- all_data_paths
#   write_path <- sprintf("%s\\figure_4_extended\\",figures_path)
#   dir.create(write_path)
#   
#   true_data_paths <- sprintf("%s\\%s", paths, "equalized")
#   
#   sessions_to_use=c(1:16)
#   all_df <- get_all_df(true_data_paths, sessions_to_use)
#   
#   subfield_idx = which(colnames(all_df) == "Subfield")
#   session_idx = which(colnames(all_df) == "Session")
#   
#   percent_active_df <- all_df[,c(7:9, session_idx, subfield_idx)]
#   percent_place_cells_df <- all_df[,c(1:3, session_idx, subfield_idx)]
#   percent_all_df <- all_df[,c(4:6,session_idx, subfield_idx)]
#   
#   colnames(percent_active_df) <- c("Both", "Right", "Left", "Session", "Subfield")
#   colnames(percent_place_cells_df) <- c("Both", "Right", "Left", "Session", "Subfield")
#   colnames(percent_all_df) <- c("Both", "Right", "Left", "Session", "Subfield")
#   
#   melted_percent_active <- melt(percent_active_df, measure.vars = c("Both", "Left", "Right"))
#   melted_percent_place_cells <- melt(percent_place_cells_df, measure.vars = c("Both", "Left", "Right"))
#   melted_percent_all <- melt(percent_all_df, measure.vars = c("Both", "Left", "Right"))
#   
#   ylabs = c("Percentage of active cells (%)",
#             "Percentage of place cells (% of active cells)",
#             "Percentage of place cells (% of all cells)")
#   
#   dfs <- list(melted_percent_active,
#               melted_percent_place_cells,
#               melted_percent_all)
#   
#   bar_plots_list <- list()
#   bar_plots_list_2 <- list()
#   for (df_idx in 1:3) {
#     
#     df_to_use <- dfs[[df_idx]]
#     
#     ca1 <- df_to_use[df_to_use$Subfield == "CA1",]
#     ca3 <- df_to_use[df_to_use$Subfield == "CA3",]
#     # print(t.test(ca3[ca3$variable == "Both",4], ca1[ca1$variable == "Both",4]))
#     # print(t.test(ca3[ca3$variable == "Left",4], ca1[ca1$variable == "Left",4]))
#     # print(t.test(ca3[ca3$variable == "Right",4], ca1[ca1$variable == "Right",4]))
#     
#     colnames(df_to_use) <- c("Session", "Subfield", "Direction", "Percent")
#     sd_df <- ddply(df_to_use, .(Subfield), 
#                    function(subfield_df) {
#                      return(ddply(subfield_df, .(Direction), 
#                                   function(sub_df) {return(c(mean(sub_df[,"Percent"]), sd(sub_df[,"Percent"])))}))
#                    })
#     
#     
#     
#     
#     df_to_use$Session <- as.character(as.numeric(df_to_use$Session) %% 8)
#     df_to_use$Session[df_to_use$Session == "0"] <- "8"
#     
#     session_mean_sd_df <- 
#       ddply(df_to_use, .(Subfield), 
#             function(subfield_df) {
#               
#               return(ddply(subfield_df, .(Session), 
#                            function(sub_df) {
#                              sd_df <- ddply(sub_df, .(Direction), 
#                                             function(sub_df_inner) 
#                                             { est <- as.numeric(sub_df_inner[,"Percent"])
#                                             return(c(mean(est), sem(est)))})
#                              return(sd_df)
#                            }));
#             });
#     
#     #
#     colnames(sd_df) <- c("Subfield", "Direction", "Percent", "Sd")
#     colnames(session_mean_sd_df) <- c("Subfield", "Session", "Direction", "Percent", "Sd")
#     
#     sd_df <- sd_df[sd_df$Direction %in% c("Right", "Left"),]
#     df_to_use <- df_to_use[df_to_use$Direction%in% c("Right", "Left"),]
#     gr <- 
#       ggplot(sd_df, aes(x=Direction, y=Percent)) + 
#       geom_jitter(data=df_to_use,
#                   aes(x=Direction, y=Percent, group=interaction(Direction, Subfield)), position=position_jitter(0.2), color="gray20", size=1, alpha=0.1) + 
#       geom_boxplot(data=df_to_use, aes(group=interaction(Direction, Subfield), fill=Subfield, color=Subfield), width=0.5, size=0.75) + 
#       theme_light() +     
#       theme(panel.grid.major = element_blank(), 
#             panel.grid.minor = element_blank(),
#             axis.line = element_line(colour = "black"),
#             panel.border = element_blank(),
#             legend.position="NA",
#             panel.background = element_blank()) + 
#       scale_color_manual(values=c(adjustcolor("gray65", alpha=1), 
#                                   adjustcolor("royalblue4", alpha=1))) + 
#       
#       scale_fill_manual(values=c(adjustcolor("gray80", alpha=0.8), 
#                                  adjustcolor("royalblue2", alpha=0.8))) +
#       xlab("Running direction") +
#       ylab(ylabs[df_idx]) + 
#       ylim(0,1)
#     
#     
#     session_mean_sd_df_f <- session_mean_sd_df[session_mean_sd_df$Direction == "Both",]
#     df_to_use_f <- df_to_use[df_to_use$Direction == "Both",]
#     gr2 <- 
#       ggplot(session_mean_sd_df_f, aes(x=Session, y=Percent)) + 
#       geom_point(data=df_to_use_f, aes(x=Session, y=Percent, group=Subfield, color=Subfield), 
#                  size=3.5, alpha=0.2, position=position_dodge(0.4)) +
#       geom_line(aes(group=Subfield, color=Subfield), size=1.5, position=position_dodge(0.4)) +
#       
#       geom_errorbar(aes(ymin=Percent - Sd, ymax=Percent + Sd, group=Subfield, color=Subfield), 
#                     size=2, width=0.3, position=position_dodge(0.4)) +
#       geom_point(aes(group=Subfield, color=Subfield), size=5, position=position_dodge(0.4)) + 
#       theme_light() +     
#       theme(panel.grid.major = element_blank(), 
#             panel.grid.minor = element_blank(),
#             axis.line = element_line(colour = "black"),
#             panel.border = element_blank(),
#             panel.background = element_blank(),
#             legend.position = "NA") + 
#       xlab("Session") +
#       ylab(ylabs[df_idx]) + 
#       ylim(0,1.1)
#     
#     
#     bar_plots_list <- append(bar_plots_list,
#                              list(gr))
#     
#     bar_plots_list_2 <- append(bar_plots_list_2,
#                                list(gr2))
#   }
#   
#   bar_plots_list$nrow <- 1
#   first_row <- do.call(grid.arrange, bar_plots_list)
#   
#   estimated_res <- get_estimated_df(all_df, true_data_paths, sessions_to_use)
#   
#   estimation_list <- estimated_res$estimation_list
#   likelihood_list <- estimated_res$likelihood_list
#   measured_percent_list <- estimated_res$measured_percent_list
#   estimated_per_percent_list <- estimated_res$estimated_per_percent_list
#   session_metavar_list <- estimated_res$session_metavar_list
#   missing_paths <- estimated_res$session_metavar_list 
#   likelihood_df <- estimated_res$likelihood_df
#   estimated_per_percent_df <- estimated_res$estimated_per_percent_df
#   session_df <- estimated_res$session_df
#   estimation_vec <- estimated_res$estimation_vec
#   measured_vec <- estimated_res$measured_vec
#   estimated_df <- estimated_res$estimated_df
#   
#   
#   both_dir_res <- get_both_dir_df(estimated_df = estimated_df, all_df=all_df)
#   
#   both_dir_df <- both_dir_res$both_dir_df
#   mean_sd_df <- both_dir_res$mean_sd_df
#   session_mean_sd_df <- both_dir_res$session_mean_sd_df
#   both_dir_df_all <- both_dir_res$both_dir_df_all
#   session_mean_sd_df_all <- both_dir_res$session_mean_sd_df_all
#   
#   # session_mean_sd_df$Session <- factor(session_mean_sd_df$Session, levels=as.character(1:16))
#   # both_dir_df$Session <- factor(both_dir_df$Session, levels=as.character(1:16))
#   mean_sd_df <- mean_sd_df[mean_sd_df$Direction %in% c("Left", "Right"),]
#   both_dir_df <- both_dir_df[both_dir_df$Direction %in% c("Left", "Right"),]
#   
#   gest <- 
#     ggplot(mean_sd_df, aes(x=Direction, y=Estimated)) + 
#     geom_boxplot(data=both_dir_df, aes(group=interaction(Direction, Subfield), fill=Subfield, color=Subfield), width=0.5, size=1) + 
#     geom_jitter(data=both_dir_df,
#                 aes(x=Direction, y=Estimated, group=interaction(Direction, Subfield)), position=position_jitter(0.2), color="gray20", size=0.75, alpha=0.1) + 
#     scale_color_manual(values=c(adjustcolor("gray65", alpha=1), 
#                                 adjustcolor("royalblue4", alpha=1))) + 
#     
#     scale_fill_manual(values=c(adjustcolor("gray80", alpha=0.8), 
#                                adjustcolor("royalblue2", alpha=0.8))) +
#     theme_light() +     
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           axis.line = element_line(colour = "black"),
#           panel.border = element_blank(),
#           legend.position = "NA",
#           panel.background = element_blank()) + 
#     xlab("Running direction") +
#     ylab(ylabs[2]) + 
#     ylim(0,1.05)
#   
#   
#   ca1_b <- both_dir_df[both_dir_df$Subfield=="CA1",]
#   ca3_b <- both_dir_df[both_dir_df$Subfield=="CA3",]
#   # t.test(ca3_b[ca3_b$Direction == "Right","Estimated"], ca1_b[ca1_b$Direction == "Right","Estimated"]) 
#   # t.test(ca3_b[ca3_b$Direction == "Both","Estimated"], ca1_b[ca1_b$Direction == "Both","Estimated"]) 
#   # t.test(ca3_b[ca3_b$Direction == "Left","Estimated"], ca1_b[ca1_b$Direction == "Left","Estimated"]) 
#   
#   
#   session_mean_sd_df <- session_mean_sd_df[session_mean_sd_df$Direction == "Both",]
#   both_dir_df <- both_dir_df[both_dir_df$Direction == "Both",]
#   
#   # session_mean_sd_df$Session <- as.character(as.numeric(session_mean_sd_df$Session) %% 8)
#   # session_mean_sd_df$Session[session_mean_sd_df$Session == "0"] <- "8"
#   
#   weighted_avg_plots_list <- list()
#   for (pn in c(colnames(both_dir_df_all)[1:3])) {
#     
#     session_mean_work_df <- session_mean_sd_df_all[,c("Subfield", "Session", "Direction", sprintf("%s_%s", pn, c("mean", "sd")))]
#     work_df<- both_dir_df_all[,c("Subfield", "Session", "Direction", pn)]
#     colnames(session_mean_work_df) <- c("Subfield", "Session", "Direction", "Mean", "Sd")
#     colnames(work_df) <- c("Subfield", "Session", "Direction", "Mean")
#     g <- 
#       ggplot(session_mean_work_df, aes(x=Session, y=Mean)) + 
#       geom_point(data=work_df, aes(x=Session, y=Mean, group=Subfield, fill="gray75", color="gray75"), 
#                  size=.5, alpha=0.2, position=position_dodge(0.4)) +
#       geom_line(aes(group=Subfield, color=Subfield), size=1, position=position_dodge(0.4)) +
#       
#       geom_errorbar(aes(ymin=Mean - Sd, ymax=Mean + Sd, group=Subfield, color=Subfield), 
#                     size=.75, width=1, position=position_dodge(0.4)) +
#       geom_point(aes(group=Subfield, color=Subfield), size=2, position=position_dodge(0.4)) + 
#       theme_light() +     
#       theme(panel.grid.major = element_blank(), 
#             panel.grid.minor = element_blank(),
#             axis.line = element_line(colour = "black"),
#             panel.border = element_blank(),
#             panel.background = element_blank(),
#             legend.position="NA") + 
#       xlab("Session") +
#       ylim(0,1.1)
#     
#     weighted_avg_plots_list[[pn]] <- g
#   }
#   
#   gest2 <- 
#     ggplot(session_mean_sd_df, aes(x=Session, y=Estimated)) + 
#     geom_point(data=both_dir_df, aes(x=Session, y=Estimated, group=Subfield), 
#                fill="gray75", color="gray75",
#                size=.5, alpha=0.2, position=position_dodge(0.4)) +
#     geom_line(aes(group=Subfield, color=Subfield), size=1, position=position_dodge(0.4)) +
#     
#     geom_errorbar(aes(ymin=Estimated - Sd, ymax=Estimated + Sd, group=Subfield, color=Subfield), 
#                   size=0.75, width=1, position=position_dodge(0.4)) +
#     geom_point(aes(group=Subfield, color=Subfield), size=2, position=position_dodge(0.4)) + 
#     theme_light() +     
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           axis.line = element_line(colour = "black"),
#           panel.border = element_blank(),
#           panel.background = element_blank(),
#           legend.position="NA") + 
#     xlab("Session") +
#     ylab(ylabs[df_idx]) + 
#     ylim(0,1.1)
#   
#   gf <- grid.arrange(bar_plots_list[[1]], gest, nrow=1)
#   
#   gf2 <- grid.arrange(bar_plots_list_2[[1]], gest2 + geom_hline(yintercept = 1), nrow=1)
#   
#   gf3 <- grid.arrange(bar_plots_list_2[[2]], bar_plots_list_2[[3]], nrow=1)
#   
#   gf4 <- grid.arrange(weighted_avg_plots_list[["Active"]], weighted_avg_plots_list[["Of_Active"]], nrow=1)
#   gf5 <- grid.arrange(weighted_avg_plots_list[["Of_All"]], gest2 + geom_hline(yintercept = 1), nrow=1)
#   
#   pdf(file=sprintf("%s\\top_row.pdf", write_path),
#       height=(2.5),
#       width=(5))
#   plot(gf)
#   dev.off()  
#   
#   pdf(file=sprintf("%s\\naive_place_cells_row.pdf", write_path),
#       height=(2.5),
#       width=(5))
#   plot(gf3)
#   dev.off()  
#   
#   
#   pdf(file=sprintf("%s\\bottom_row.pdf", write_path),
#       height=(2.5),
#       width=(5))
#   plot(gf2)
#   dev.off()  
#   
#   pdf(file=sprintf("%s\\weighted_active_of_active.pdf", write_path),
#       height=(2.5),
#       width=(5))
#   plot(gf4)
#   dev.off()  
#   
#   pdf(file=sprintf("%s\\weighted_of_all_estimated.pdf", write_path),
#       height=(2.5),
#       width=(5))
#   plot(gf5)
#   dev.off()
#   
#   likdf_tmp <- likelihood_df
#   likelihood_df <- (t(apply(likelihood_df, 1, function(r) {return(r/sum(r))})))
#   likelihood_df <- as.data.frame(likelihood_df)
#   likelihood_sessions <- as.numeric(session_df[,1])
#   likelihood_sessions[likelihood_sessions > 8] <-  likelihood_sessions[likelihood_sessions > 8] - 8
#   likelihood_df <- cbind(likelihood_df, likelihood_sessions)
#   colnames(likelihood_df) <- c(as.character(seq(0.5,1,by=0.1)), "sessions")
#   
#   
#   plots_list <- list()
#   for (sf in c("CA1", "CA3")) { 
#     
#     op_likelihood_df <- likelihood_df[which(session_df[,4] == sf),]
#     print(dim(op_likelihood_df))
#     
#     melted_likelihood <- melt(op_likelihood_df, id.vars="sessions")
#     colnames(melted_likelihood) <- c("Session", "Simulated percent", "Normalized likelihood")
#     
#     likelihood_mean_sd_df <- ddply(melted_likelihood, .(`Simulated percent`), 
#                                    function(sub_df) {
#                                      ddply(sub_df, .(`Session`), 
#                                            function(sub_df_inner) {
#                                              mli <- sub_df_inner[,"Normalized likelihood"]
#                                              return(c(mean(mli), sem(mli)))})
#                                    })
#     
#     colnames(likelihood_mean_sd_df) <- c("Simulated percent", "Session", "Mean", "Sd")
#     
#     melted_likelihood$`Simulated percent` <- as.numeric(as.character(melted_likelihood$`Simulated percent`))
#     melted_likelihood$`Simulated percent` <- factor(melted_likelihood$`Simulated percent` * 100,
#                                                     as.character(seq(50,100, by=10)))
#     
#     likelihood_mean_sd_df$`Simulated percent` <- as.numeric(as.character(likelihood_mean_sd_df$`Simulated percent`))
#     likelihood_mean_sd_df$`Simulated percent` <- factor(likelihood_mean_sd_df$`Simulated percent` * 100, 
#                                                         as.character(seq(50,100, by=10)))
#     fdf <- likelihood_mean_sd_df[likelihood_mean_sd_df$Session %in% c(1:2,7:8),]
#     glikelihood <- 
#       ggplot(fdf, aes(x=`Simulated percent`, y=Mean, group=Session)) + 
#       geom_bar(stat="summary", aes(group=`Session`, fill=Session), color="NA", 
#                position = position_dodge(), width=0.5) + 
#       geom_errorbar(aes(ymin=Mean, ymax=Mean + Sd, group=`Session`),
#                     width=0.5, position=position_dodge()) +
#       # geom_jitter(data=melted_likelihood,
#       #             aes(x=`Simulated percent`, y=`Normalized likelihood`, group=`Session`, fill=Session), position=position_jitterdodge(0), color="gray40", size=3, alpha=0.5) + 
#       scale_fill_distiller(palette="RdYlBu") +
#       theme_light() +     
#       theme(panel.grid.major = element_blank(), 
#             panel.grid.minor = element_blank(),
#             axis.line = element_line(colour = "black"),
#             panel.border = element_blank(),
#             panel.background = element_blank(),
#             legend.position = "NA") + 
#       xlab("Simulated percent (%)") +
#       ylab("Normalized likelihood") + 
#       ylim(0,max(likelihood_mean_sd_df[,"Mean"]) + 0.02) +
#       ggtitle(sf)
#     plots_list <- append(plots_list, list(glikelihood))
#     
#   }
#   
#   plots_list$nrow = 1
#   glike_all <-   do.call(grid.arrange, plots_list)
#   
#   pdf(file=sprintf("%s\\likelihood_row.pdf", write_path),
#       height=(2.5),
#       width=(5))
#   plot(glike_all)
#   dev.off()  
#   
# }
# 
# figure_2_top_row <- function(path, session) {
#   write_path <- sprintf("%s\\figure_2\\",figures_path)
#   dir.create(write_path)  
#   
#   
#   dir.create(sprintf("%s\\big_histograms", write_path))
#   dir.create(sprintf("%s\\small_histograms", write_path))
#   dir.create(sprintf("%s\\estimations", write_path))
#   
#   for (np in all_data_paths) { 
#     for (ses in 1:16) {
#       
#   session_path <- sprintf("%s\\equalized\\session_%d", np, ses)
#   b <- get_spike_train_and_stim_trace_from_path(np, ses, equalize_prior=T)
#   
#   params_ex <- get_fit_params(session_path)
#   params <- params_ex$params
#   
#   print(params_ex)
#   
#   
#   
#   
#   graphs <- plot_gen_place_cell_histograms_new(average_width = params["average_width"],
#                                                sd_width = params["sd_width"],
#                                                noise_level = params["noise"],
#                                                double_peak_percent = params["double_peak_pct"],
#                                                true_spike_train = b[[1]],
#                                                stim_trace = b[[2]],
#                                                place_cell_percentage = params_ex$pct,
#                                                smooth_lambda = -1,
#                                                breaks_peaks=9,
#                                                breaks_SI=21,
#                                                pval_func = ks.test,
#                                                density_p = T,
#                                                title=T)
# 
#   
#   graphs_t <- lapply(graphs, function(pl) {pl + theme(plot.margin=unit(c(0,0,0,0), "cm"))})
#   small_graphs_t  <- lapply(graphs, function(pl) {pl + theme(plot.margin=unit(c(0,0,0,0), "cm")) + ggtitle("")})
#   graphs_t$nrow <- 2
#   small_graphs_t$nrow  <- 1
#   
#   hists_p <- do.call(arrangeGrob, graphs_t)
#   small_hists_p <- do.call(arrangeGrob, small_graphs_t)
#   
#   if (grepl("CA3", np)) {
#     mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', np))
#     mice_str <- substr(np, mice_str_index, mice_str_index+4) 
#   } else {
#     mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', np))
#     mice_str <- substr(np, mice_str_index, mice_str_index+3) 
#   }
#   
#   dir <- ifelse(grepl(np, "Right"), "R", "L")
#   
#   pdf(sprintf("%s\\big_histograms\\simulation_hists_%s_%d_%s.pdf",write_path, mice_str, ses, dir), 
#       height=4, width=4); 
#   plot(hists_p); 
#   dev.off() 
#   pdf(sprintf("%s\\small_histograms\\simulation_hists_%s_%d_%s.pdf",write_path, mice_str, ses, dir), 
#       height=1.5, width=6); 
#   plot(small_hists_p); 
#   dev.off() 
#     #pdf(sprintf("%s\\simulation_hists_%s.pdf",write_path, ext), height=3, width=3); plot(hists_p); dev.off() 
#     }
#   }
#   
#   for (np in all_data_paths) { 
#     for (ses in 1:16) {
#       
#   session_path <- sprintf("%s\\equalized\\session_%d", np, ses) 
#   likelihood_df <- get_fit_params(session_path, just_likelihood = T)
#   
#   measured_percent <- likelihood_df[1,"measured_pct"]
#   
#   likelihood_vec <- likelihood_df[,"likelihood"]
#   likelihood_df <- likelihood_df[,-which(colnames(likelihood_df) == "likelihood")] # Remove last
#   likelihood_df <- likelihood_df[,-which(colnames(likelihood_df) == "measured_pct")] # Remove last
#   
#   
#   # This is just to re-extract measured place cell percentage,
#   # In older versions of the analysis this was saved with the df as well
#   rownames(likelihood_df) <- rev(rownames(likelihood_df))
#   names(likelihood_vec) <- rownames(likelihood_df)
#   mdf <- melt(likelihood_df)
#   mdf$Var1 <- factor(mdf$Var1)
#   mdf$SimP <- factor(sprintf("%.2f",as.numeric(as.vector(mdf$Var1))),
#                      levels=sprintf("%.2f",seq(0.4,1,by=0.1)))
#   
#   mdf$MeasuredP <- mdf$value
#   
#   bar_df <- data.frame(Likelihood=likelihood_vec, 
#                        Percentage=factor(sprintf("%.2f",as.numeric(names(likelihood_vec))),
#                                          levels=sprintf("%.2f",seq(0.4,1,by=0.1))))
#   
#   bar_df$Likelihood <- bar_df$Likelihood / sum(bar_df$Likelihood)
#   
#   limits_of_plots <- c(min(mdf$value) - 0.05, max(mdf$value) + 0.1)
#   gbox <- 
#     ggplot(mdf) +
#     geom_violin(aes(y=MeasuredP, group=SimP, fill=SimP, x=SimP), size=1, alpha=0.5) + 
#     geom_jitter(aes(y=MeasuredP, group=SimP, x=SimP), size=0.05, position=position_jitter(0.15), alpha=0.25)  +
#     stat_summary(fun=mean, geom="point", size=1.5, color="red", aes(x=SimP, y=MeasuredP)) + 
#     scale_fill_brewer(palette="RdYlBu") +
#     theme_light() +
#     xlab("True simulation percentage (%)") +
#     ylab("Measured percentage (%)") +
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           axis.line = element_line(colour = "black"),
#           axis.text.x = element_text(angle = 45, vjust=0.5),
#           legend.position="NA",
#           panel.border = element_blank(),
#           panel.background = element_blank()) +
#     geom_hline(yintercept=(measured_percent), linetype="dashed", col="gray60", size=2) +
#     geom_text(y=(measured_percent + 0.01), x=2.5, label="Measured place cell percentage (actual data)",col="gray60",size=3) 
#   
#   gbar <- ggplot(bar_df, aes(x=Percentage, y=Likelihood)) +
#     geom_bar(stat="identity", fill="gray50", color="black", width=0.5) + 
#     theme_light() +
#     xlab("True simulation percentage (%)") +
#     ylab("Normalized Likelihood (AU)") +
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           axis.text.x = element_text(angle = 45, vjust=0.5),
#           axis.line = element_line(colour = "black"),
#           legend.position="NA",
#           panel.border = element_blank(),
#           panel.background = element_blank())
#   
#   margin <- theme(plot.margin=unit(c(0,0,0,0,0), "cm"))
#   
#   rest <- grid.arrange(gbox + margin, gbar + margin, nrow=1)
#   
#   if (grepl("CA3", np)) {
#     mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', np))
#     mice_str <- substr(np, mice_str_index, mice_str_index+4) 
#   } else {
#     mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', np))
#     mice_str <- substr(np, mice_str_index, mice_str_index+3) 
#   }
#   
#   dir <- ifelse(grepl("Right", np), "R", "L")
#   
#   pdf(sprintf("%s\\estimations\\box_and_bar_%s_%d_%s.pdf",write_path, mice_str, ses, dir), 
#       height=2.5, width=5); 
#   plot(rest); 
#   dev.off() 
#     }
#   }
#   
# }
# 
# 
# figure_2_middle_row <- function(path, session) {
#   
#   write_path <- sprintf("%s\\figure_2\\",figures_path)
#   dir.create(write_path)  
#   
#   session_path <- sprintf("%s\\equalized\\session_%d", path, session)
#   
#   load(file=sprintf("%s\\sim_vs_real\\sim_vs_real_by_sample.R", session_path))
#   simulation_vs_real_df <- final_result
#   names(simulation_vs_real_df)[[7]] <- "real"
#   
#   g_real_cyc <-  ggplot(simulation_vs_real_df$real, aes(x=Time, y=Cyclic)) + geom_line(size=2) + 
#     theme_light() +
#     ylim(c(0, 0.85)) +
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           axis.line = element_line(colour = "black"),
#           legend.position="NA",
#           panel.border = element_blank(),
#           panel.background = element_blank(),
#           plot.margin=unit(c(0,0,0,0), "cm")) +
#     labs(x="", y="")
#   #labs(x="Sample duration (minutes)", y = "Fraction of place cells (%)") +
#   #xlim(0, max(final_result[[7]][,"Time"]))
#   
#   
#   lev <- rev(as.numeric(names(simulation_vs_real_df)[-len(names(simulation_vs_real_df))]))
#   
#   sim_pct_colors = brewer.pal(n = len(lev),  name = "Spectral")
#   names(sim_pct_colors) <- lev
#   pval_df <- c()
#   plots_list <- list()
#   
#   for (sim_pct in lev) {
#     print(sim_pct)
#     simulation_df <- final_result[[as.character(sim_pct)]]
#     g_sim_cyc <- g_real_cyc# + ggtitle(sprintf("%d%%", as.numeric(sim_pct) * 100))
#     pval_vec <- c()
#     
#     for (i in 1:(ncol(simulation_df) / 2)) {
#       tmp_df_cyc <- data.frame(Fraction=simulation_df[, (i * 2)],
#                                Time=simulation_vs_real_df$real$Time)
#       
#       g_sim_cyc <- 
#         g_sim_cyc + geom_line(data=tmp_df_cyc,
#                               aes(x=Time, y=Fraction),
#                               color=sim_pct_colors[as.character(sim_pct)],
#                               size=2,
#                               alpha=0.5)
#       
#       
#       pval_vec <- c(pval_vec,
#                     wilcox.test(simulation_df[, (i * 2)],
#                                 simulation_vs_real_df$real$Cyclic)$p.value)
#       
#     }
#     
#     pval_df <- rbind(pval_df, pval_vec)
#     plots_list <- append(plots_list, list(g_sim_cyc))
#   }
#   
#   plots_list <- lapply(plots_list, function(pl) {pl + theme(plot.margin=unit(c(0,0,0,0), "cm"))})
#   plots_list$nrow <- 2
#   plots_list$left <-"Fraction of place cells (%)"
#   #plots_list$bottom <- "Sample duration (minutes)"
#   sim_real_p <- do.call(grid.arrange, plots_list)
#   pval_violin <- figure_2_pvalue_violin(new_paths[8:1],  use_JSD=T, pval_func = t.test,)
#   
#   grid.arrange(sim_real_p, pval_violin[[1]], nrow=1)
#   
#   pdf(sprintf("%s\\sim_vs_real_%s.pdf",write_path, ext), height=2.5, width=4.5); plot(sim_real_p); dev.off() 
#   pdf(sprintf("%s\\pval_violin_all.pdf",write_path), height=2.5, width=3.5); plot(pval_violin[[1]] + theme(#plot.margin=unit(c(0,0,0,0), "cm"), 
#     axis.text.x = element_text(angle = 0, vjust=1),
#     legend.position="NA")) + ylab(""); dev.off()
#   
#   pdf(sprintf("%s\\pval_violin_all_legend.pdf",write_path), height=2.5, width=3.5); plot(pval_violin[[1]] + theme(#plot.margin=unit(c(0,0,0,0), "cm"), 
#     axis.text.x = element_text(angle = 0, vjust=1))) + ylab(""); dev.off()
#   
#   pdf(sprintf("%s\\pval_violin_compare.pdf",write_path), height=2.5, width=3.5); plot(pval_violin[[2]] + theme(#plot.margin=unit(c(0,0,0,0), "cm"), 
#     axis.text.x = element_text(angle = 0, vjust=1),
#     legend.position="NA")) + ylab(""); dev.off()
#   
#   pdf(sprintf("%s\\pval_violin_compare_legend.pdf",write_path), height=2.5, width=3.5); plot(pval_violin[[2]] + theme(#plot.margin=unit(c(0,0,0,0), "cm"), 
#     axis.text.x = element_text(angle = 0, vjust=1))) + ylab(""); dev.off()
# }
# 
# figure_2_bottom_row <- function() {
#   
#   true_data_paths <- sprintf("%s\\%s", all_data_paths[c(1:18)], "equalized")
#   ses_ind <- c(1:16)
#   #by_activity <- plot_place_cell_pct_by_activity_plot(true_data_paths, cyclic_only=T, fit_sim = "fit_simulation_new_v_fit\\")
#   
#   by_duration <- plot_place_cell_pct_by_sample_duration_figure(ses_ind = ses_ind, 
#                                                                true_data_paths, 
#                                                                bin_size = 2.5, 
#                                                                ext = "",
#                                                                prefix = "JSD_2_",
#                                                                cyclic_only=T, 
#                                                                fit_sim = "fit_simulation\\", verbose=F)
#   by_pval <- plot_pval_by_sample_duration_figure(true_data_paths, ext = "", fit_sim = "fit_simulation\\", verbose=F, type1=F, ses_ind = ses_ind)
#   
#   #by_activity_real <- plot_place_cell_pct_by_activity_plot(true_data_paths, cyclic_only=T)
#   by_duration_real <- plot_place_cell_pct_by_sample_duration_figure(true_data_paths, bin_size = 2.5, 
#                                                                     ext = "", 
#                                                                     cyclic_only=T, 
#                                                                     fit_sim = "", 
#                                                                     verbose=F, 
#                                                                     ses_ind = ses_ind)
#   by_pval_real <- plot_pval_by_sample_duration_figure(true_data_paths, ext = "", fit_sim = "", verbose=F, type1 = F, ses_ind = ses_ind)
#   overlaid_duration <- by_duration_real[[1]] + geom_line(data=by_duration[[4]], linetype="dashed", aes(x=time, y=fraction), size=1, color="black") + geom_errorbar(size=0.8,data=by_duration[[4]], aes(ymin=fraction-sd, ymax=fraction+sd), color=adjustcolor("black", alpha=1), width=0.5)
#   
#   overlaid_pval <- by_pval_real[[1]] + geom_line(data=by_pval[[2]], aes(x=time, y=pvd), size=1, color="black") + 
#     geom_errorbar(size=0.3, data=by_pval[[2]], aes(ymin=pvd-pvdsd, ymax=pvd+pvdsd), color=adjustcolor("black", alpha=1), width=0.5)
#   
#   ap <- align_plots(overlaid_duration, overlaid_pval, align="v")
#   grid.arrange(ap[[1]], ap[[2]], nrow=2)
#   
#   colnames(by_pval[[3]]) <- 1:ncol(by_pval[[3]])
#   colnames(by_pval_real[[3]]) <- 1:ncol(by_pval_real[[3]])
#   melted_1 <- melt(by_pval[[3]])
#   melted_2 <- melt(by_pval_real[[3]])
#   colnames(melted_1) <- c("instance", "time", "pval")
#   colnames(melted_2) <- c("instance", "time", "pval")
#   
#   pval_df <- rbind(melted_1, melted_2)
#   pval_df <- as.data.frame(pval_df)
#   pval_df$group <- c(rep("real", nrow(melted_1)), rep("sim", nrow(melted_2)))
#   
#   
#   tmp_df <- t(by_activity[[4]])
#   tmp_df_real <- t(by_activity_real[[4]])
#   rownames(tmp_df) <- c()
#   rownames(tmp_df_real) <- c()
#   
#   melted_1 <- melt(tmp_df)
#   melted_2 <- melt(tmp_df_real)
#   colnames(melted_1) <- c("instance", "timebins", "fraction")
#   colnames(melted_2) <- c("instance", "timebins", "fraction")
#   
#   activity_df <- rbind(melted_1[,2:3], melted_2[,2:3])
#   activity_df <- as.data.frame(activity_df)
#   activity_df$group <- c(rep("real", nrow(melted_1)), rep("sim", nrow(melted_2)))
#   activity_df <- activity_df[!is.na(activity_df$fraction),]
#   
#   
#   tmp_df <- t(by_duration[[5]])
#   tmp_df_real <- t(by_duration_real[[5]])
#   rownames(tmp_df) <- c()
#   rownames(tmp_df_real) <- c()
#   
#   melted_1 <- melt(tmp_df[-1,])
#   melted_2 <- melt(tmp_df_real[-1,])
#   colnames(melted_1) <- c("instance", "time", "fraction")
#   colnames(melted_2) <- c("instance", "time", "fraction")
#   
#   duration_df <- rbind(melted_1[,2:3], melted_2[,2:3])
#   duration_df <- as.data.frame(duration_df)
#   duration_df$group <- c(rep("real", nrow(melted_1)), rep("sim", nrow(melted_2)))
#   duration_df <- duration_df[!is.na(duration_df$fraction),]
#   duration_df$fraction <- as.numeric(duration_df$fraction)
#   
#   
#   g <- 
#     grid.arrange(by_duration[[1]],
#                  by_activity[[1]],
#                  by_pval[[1]],
#                  nrow=1)
#   
#   #by_activity_real[[5]] + geom_line(data=by_activity[[3]], aes(x=activity, y=fraction)) + geom_errorbar(data=by_activity[[3]], aes(ymin=fraction-sem, ymax=fraction+sem)) + xlim(c(-1,135))
#   
#   overlaid_activity <- by_activity_real[[5]] + geom_line(data=by_activity[[3]], aes(x=activity, y=fraction)) + geom_errorbar(data=by_activity[[3]], aes(ymin=fraction-sem, ymax=fraction+sem), width=2.5, col="black")
#   overlaid_duration <- by_duration_real[[6]] + geom_line(data=by_duration[[4]], aes(x=time, y=fraction), size=1, color="black") + geom_errorbar(size=1,data=by_duration[[4]], aes(ymin=fraction-sd, ymax=fraction+sd), color=adjustcolor("black", alpha=1), width=0.5)
#   overlaid_pval <- by_pval_real[[1]] + geom_line(data=by_pval[[2]], aes(x=time, y=pvd), size=1, color="black") + 
#                                   geom_ribbon(data=by_pval[[2]], aes(ymin=pvd-pvdsd, ymax=pvd+pvdsd), fill=adjustcolor("black", alpha=0.2), color="NA")
#   
#   
#   goverlaid <- grid.arrange(overlaid_duration,
#                             overlaid_activity,
#                             overlaid_pval,
#                             nrow=1)
#   
#   #
#   for (dur in unique(duration_df$time)){
#     
#     ttest <- t.test(duration_df[duration_df$time == dur & duration_df$group == "real", "fraction"],
#                     duration_df[duration_df$time == dur & duration_df$group == "sim", "fraction"],
#                     paired=T)
#     wlk <- wilcox.test(duration_df[duration_df$time == dur & duration_df$group == "real", "fraction"],
#                        duration_df[duration_df$time == dur & duration_df$group == "sim", "fraction"],
#                        paired=T)
#     
#     print("_______________")
#     print(ttest$statistic)
#     print(wlk$statistic)
#     print("##################")
#     print(ttest$p.value)
#     print(wlk$p.value)
#     print("_______________")
#   }
#   
#   
#   
#   #
#   for (act in unique(activity_df$timebins)){
#     
#     ttest <- t.test(activity_df[activity_df$timebins == act & activity_df$group == "real", "fraction"],
#                     activity_df[activity_df$timebins == act & activity_df$group == "sim", "fraction"])
#     wlk <- wilcox.test(activity_df[activity_df$timebins == act & activity_df$group == "real", "fraction"],
#                        activity_df[activity_df$timebins == act & activity_df$group == "sim", "fraction"])
#     
#     print("_______________")
#     print(ttest$statistic)
#     print(wlk$statistic)
#     print("##################")
#     print(ttest$p.value)
#     print(wlk$p.value)
#     print("_______________")
#   }
#   
#   for (dur in unique(pval_df$time)){
#     
#     ttest <- t.test(pval_df[pval_df$time == dur & pval_df$group == "real", "pval"],
#                     pval_df[pval_df$time == dur & pval_df$group == "sim", "pval"],
#                     paired=T)
#     wlk <- wilcox.test(pval_df[pval_df$time == dur & pval_df$group == "real", "pval"],
#                        pval_df[pval_df$time == dur & pval_df$group == "sim", "pval"],
#                        paired=T)
#     
#     print("_______________")
#     print(ttest$statistic)
#     print(wlk$statistic)
#     print("##################")
#     print(ttest$p.value)
#     print(wlk$p.value)
#     print("_______________")
#   }
#   
#   pdf(sprintf("%s\\tmp_cosyne,pdf",write_path),height=1.9 * 4/3,width=1.9 * 4/3)
#   pdf(sprintf("%s\\bottom_row.pdf",write_path), height=1.9 * 4/3, width=1.9 * 4); plot(g); dev.off()
#   pdf(sprintf("%s\\bottom_row_overlaid.pdf",write_path), height=1.9 * 4/3, width=1.9 * 4); plot(goverlaid); dev.off()
#   
# }
# 
# get_spike_train_and_stim_trace_from_path <- function(path, session, old=F, equalize_frames=T, simulate=F) {
#   
#   
#   if (old) {
#     m <- readMat(paste(path, "\\stimulus_trace.mat", sep=""))
#     stim_trace <- m$stimulus.trace
#     
#     
#     m <- readMat(paste(path, spike_train_filename, sep=""))
#     spike_train <- m$spike.train
#   } else {  
#     m <- readMat(paste(path, stimulus_trace_filename, sep=""))
#     
#     if(len(grep("Left", path)) > 0) {
#       stim_trace <- m$position.left.trials
#     } else {
#       stim_trace <- m$position.right.trials
#     }
#     
#     
#     m <- readMat(paste(path, spike_train_filename, sep=""))
#     
#     if(len(grep("Left", path)) > 0) {
#       spike_train <- m$spikes.left.trials
#     } else {
#       spike_train <- m$spikes.right.trials
#     }
#   }
#   
#   tmp_path <- path
#   tmp_spike_train <- spike_train
#   tmp_stimulus_trace <- stim_trace
#   
#   
#   spike_train <- tmp_spike_train
#   stim_trace <- tmp_stimulus_trace
#   spike_train <- spike_train[[session]][[1]]
#   spike_train <- t(as.matrix(spike_train))
#   stim_trace <- stim_trace[[session]][[1]]
#   stim_trace <-  as.vector(stim_trace)
#   
#   
#   
#   if (!old && equalize_frames ) {
#     equalize_bins <- c(1:3,22:24)
#     indices <- 1:len(stim_trace)
#     
#     
#     final_ind <- rep(T, times=len(stim_trace))
#     
#     for (ebin in equalize_bins) {
#       print(ebin)
#       
#       num_to_draw <- gen_num_frames(stim_trace)
#       if (num_to_draw > len(which(stim_trace == ebin))) {
#         
#         next
#       }
#       
#       sampled_ind <- sample(which(stim_trace==ebin), num_to_draw)
#       tmp <- indices %in% sampled_ind | stim_trace[indices] != ebin       
#       final_ind <- final_ind & tmp
#     }
#     
#     spike_train <- spike_train[,final_ind]
#     stim_trace <- stim_trace[final_ind]
#   }
#   
#   
#   
#   if (simulate) {
#     fr <- rowMeans(spike_train) / dt
#     processed_real <- preprocess_spike_train(spike_train, stim_trace)
#     
#     true_cells_spike_train <- processed_real$working_cells_spike_train
#     true_firing_rate <- processed_real$working_firing_rate
#     true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
#     
#     print(sprintf("%s\\equalized\\session_%d", path, 2))
#     params <- get_fit_params_figures_from_path(sprintf("%s\\equalized\\session_%d", path, session) , T)
#     simulated_tuning_curve <-
#       generate_tuning_curves_cost(n = nrow(spike_train),
#                                   percentage = params$pct,
#                                   average_width = params$params["average_width"], 
#                                   sd_width = params$params["sd_width"],
#                                   fixed_fr=fr,
#                                   noise=params$params["noise"],
#                                   double_peak_pct = params$params["double_peak_pct"],
#                                   n_bins=20,
#                                   plot=F)
#     
#     pois_factor <- currate_spike_train_cost(simulated_tuning_curve, 
#                                             true_time_bins_per_cells,
#                                             stim_trace,
#                                             verbose = T)
#     
#     generated_spike_train <- 
#       generate_spike_trains_cost(tuning_curves = simulated_tuning_curve,
#                                  stim_trace = stim_trace,
#                                  factor=pois_factor,
#                                  fs=1)
#     
#     orig_st <- spike_train
#     orig_path <- path
#     spike_train <- generated_spike_train
#     return (list(spike_train, stim_trace, simulated_tuning_curve, orig_st))
#   }
#   
#   
#   return (list(spike_train, stim_trace))
# }
# 
# 
# 
# get_cost_compare_vectors_by_path <- function(path, true_spike_train, stim_trace, pval_func=wilcox.test, ret=F, smooth_lambda=-1, best=T, new_simulation=F, use_JSD=F, edge_free=F) {
#   
#   firing_rate <-rowMeans(true_spike_train) / dt
#   processed_real <- preprocess_spike_train(true_spike_train, stim_trace=stim_trace, verbose=F)
#   true_cells_spike_train <- processed_real$working_cells_spike_train
#   true_firing_rate <- processed_real$working_firing_rate
#   true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
#   
#   tmp <- compute_tuning(true_cells_spike_train, stim_trace)
#   stim_prob <- tmp[[1]]
#   true_tuning_curve <- tmp[[2]]
#   rm(tmp)
#   
#   
#   if (smooth_lambda == - 1){
#     tuning_for_peaks <- true_tuning_curve
#   } else {
#     tuning_for_peaks <- 
#       t(apply(true_tuning_curve, 1, function(trc) {
#         barp <- barplot(trc, plot=F)
#         smt <- smooth.spline(barp, trc, all.knots = T, lambda=smooth_lambda)
#         smt$y[smt$y < 0] <- 0
#         return(smt$y);
#       }))
#   }
#   
#   true_peaks <- unlist(apply(tuning_for_peaks, 1, 
#                              function(n){return(len(get_peaks(n, threshold_from_max = peaks_significance_threshold)))}))
#   true_active_spatial_bins_per_cell <- apply(true_tuning_curve, 1, function(n){sum(n>0)})
#   #print(true_active_spatial_bins_per_cell)
#   true_SI <- compute_SI(stim_prob, true_tuning_curve, true_firing_rate)
#   
#   params <- get_fit_params_figures_from_path(path, best = best, use_JSD = use_JSD, edge_free = edge_free)
#   SI_f <- c()
#   tb_f <- c()
#   peaks_f <- c()
#   sb_f <- c()
#   simulated_tuning_curve <-
#     generate_tuning_curves_cost(n = nrow(true_spike_train),
#                                 percentage = params$pct,
#                                 average_width = params$params["average_width"], 
#                                 sd_width = params$params["sd_width"],
#                                 fixed_fr=firing_rate,
#                                 noise=params$params["noise"],
#                                 double_peak_pct = params$params["double_peak_pct"],
#                                 n_bins = ifelse(edge_free, 24, 20),
#                                 plot=F)
#   pois_factor <- currate_spike_train_cost(simulated_tuning_curve, 
#                                           true_time_bins_per_cells,
#                                           stim_trace,
#                                           verbose = F)
#   
#   
#   for (i in 1:10) { 
#     simulated_tuning_curve <-
#       generate_tuning_curves_cost(n = nrow(true_spike_train),
#                                   percentage = params$pct,
#                                   average_width = params$params["average_width"], 
#                                   sd_width = params$params["sd_width"],
#                                   fixed_fr=firing_rate,
#                                   noise=params$params["noise"],
#                                   double_peak_pct = params$params["double_peak_pct"],
#                                   plot=F)
#     
#     
#     generated_spike_train <- 
#       generate_spike_trains_cost(tuning_curves = simulated_tuning_curve,
#                                  stim_trace = stim_trace,
#                                  factor=pois_factor,
#                                  fs=1)
#     
#     
#     processed_generated <- preprocess_spike_train(generated_spike_train, stim_trace, verbose=F)
#     generated_cells_active_st <- processed_generated$working_cells_spike_train
#     generated_firing_rate <- processed_generated$working_firing_rate
#     generated_time_bins <- processed_generated$working_time_bins_per_cells
#     
#     tmp <- compute_tuning(generated_cells_active_st, stim_trace)
#     gen_stim_prob <- tmp[[1]]
#     gen_tuning_curve <- tmp[[2]]
#     rm(tmp)
#     
#     if (smooth_lambda == - 1){
#       tuning_for_peaks <- gen_tuning_curve
#     } else {
#       tuning_for_peaks <- 
#         t(apply(gen_tuning_curve, 1, function(trc) {
#           barp <- barplot(trc, plot=F)
#           smt <- smooth.spline(barp, trc, all.knots = T, lambda=smooth_lambda)
#           smt$y[smt$y < 0] <- 0
#           return(smt$y);
#         }))
#     }
#     
#     generated_peaks <- unlist(apply(tuning_for_peaks, 1, 
#                                     function(n){return(len(get_peaks(n, threshold_from_max = peaks_significance_threshold)))}))
#     generated_active_spatial_bins_per_cell <- apply(gen_tuning_curve, 1, function(n){sum(n>0)})
#     generated_SI <- compute_SI(gen_stim_prob, gen_tuning_curve, generated_firing_rate)
#     
#     
#     if (ret) {
#       return(list(true_SI=true_SI[[1]],
#                   true_tb=true_time_bins_per_cells,
#                   true_peaks=true_peaks,
#                   true_sb=true_active_spatial_bins_per_cell,
#                   gen_SI=generated_SI[[1]],
#                   gen_tb=generated_time_bins,
#                   gen_peaks=generated_peaks,
#                   gen_sb=generated_active_spatial_bins_per_cell))
#     }
#     
#     SI=pval_func(true_SI[[1]], generated_SI[[1]])$p.value
#     tb=pval_func(true_time_bins_per_cells, generated_time_bins)$p.value
#     peaks=pval_func(true_peaks, generated_peaks)$p.value
#     sb=pval_func(true_active_spatial_bins_per_cell, generated_active_spatial_bins_per_cell)$p.value
#     
#     print(sprintf("SI %.3f, tb %.3f, peaks %.3f, sb %.3f",
#                   SI,
#                   tb,
#                   peaks,
#                   sb))
#     SI_f <- c(SI_f, SI)
#     tb_f <- c(tb_f, tb)
#     peaks_f <- c(peaks_f, peaks)
#     sb_f <- c(sb_f, sb)
#     
#   }
#   
#   return(c(mean(SI_f), mean(tb_f), mean(peaks_f), mean(sb_f)))
#   
# }
# 
# # 
# 
# 
# figure_2_pvalue_violin <- function(paths, sessions_to_use=c(5:8,13:16), verbose=F, use_JSD=F, pval_func=wilcox.test, edge_free=F) {
#   stimulus_trace_filename = "\\stim_trace.mat"
#   spike_train_filename = "\\spike_train.mat"
#   
#   equalized_paths <- sprintf("%s\\%s", paths, "equalized")
#   
#   most_likely <- list()
#   least_likely <- list()
#   
#   for (path_idx in 1:len(paths)) {
#     
#     
#     working_path <- paths[path_idx]
#     working_path <- paths[path_idx]
#     equalized_path <- equalized_paths[path_idx]
#     
#     m <- readMat(paste(working_path, stimulus_trace_filename, sep=""))
#     
#     if(len(grep("Left", working_path)) > 0) {
#       stim_trace <- m$position.left.trials
#     } else {
#       stim_trace <- m$position.right.trials
#     }
#     
#     
#     m <- readMat(paste(working_path, spike_train_filename, sep=""))
#     
#     if(len(grep("Left", working_path)) > 0) {
#       spike_train <- m$spikes.left.trials
#     } else {
#       spike_train <- m$spikes.right.trials
#     }
#     
#     
#     session_paths <- 
#       list.dirs(equalized_path, recursive=F)[grep("session", 
#                                                   list.dirs(equalized_path, 
#                                                             recursive=F))]
#     
#     
#     tmp_spike_train <- spike_train
#     tmp_stimulus_trace <- stim_trace
#     
#     # Run through all paths
#     for (tpath in session_paths) {
#       tmp_vec <- c()
#       tmp_vec_l <- c()
#       
#       for (jx in 1:5) { 
#         
#         
#         ses_idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
#         
#         # Ignore session
#         if (sum(ses_idx == sessions_to_use) == 0){
#           next
#         }
#         
#         spike_train <- tmp_spike_train
#         stim_trace <- tmp_stimulus_trace
#         spike_train <- spike_train[[ses_idx]][[1]]
#         spike_train <- t(as.matrix(spike_train))
#         stim_trace <- stim_trace[[ses_idx]][[1]]
#         stim_trace <-  as.vector(stim_trace)
#         
#         ###### equalize spike train
#         equalize_bins <- c(1:3,22:24)
#         indices <- 1:len(stim_trace)
#         
#         
#         final_ind <- rep(T, times=len(stim_trace))
#         
#         for (ebin in equalize_bins) {
#           print(ebin)
#           
#           num_to_draw <- gen_num_frames(stim_trace)
#           if (num_to_draw > len(which(stim_trace == ebin))) {
#             print(sprintf("No need to equalize bin %d", ebin))
#             next
#           }
#           
#           sampled_ind <- sample(which(stim_trace==ebin), num_to_draw)
#           tmp <- indices %in% sampled_ind | stim_trace[indices] != ebin       
#           final_ind <- final_ind & tmp
#         }
#         
#         spike_train <- spike_train[,final_ind]
#         stim_trace <- stim_trace[final_ind]
#         
#         if (edge_free) {
#           print("REMOVING EDGES!")
#           non_edges_ind <- which(stim_trace %in% c(3:22))
#           stim_trace <- stim_trace[non_edges_ind]
#           spike_train <- spike_train[,non_edges_ind]
#           stim_trace <- stim_trace - 2 # Refer to bin 3 as bin 1
#         }
#         
#         most <- get_cost_compare_vectors_by_path(tpath, spike_train, stim_trace, ret=F, smooth_lambda = 1e-5, pval_func = pval_func, use_JSD = use_JSD, edge_free=edge_free)
#         least <- get_cost_compare_vectors_by_path(tpath, spike_train, stim_trace, ret=F, smooth_lambda = 1e-5, pval_func = pval_func, best = F, use_JSD = use_JSD, edge_free=edge_free)
#         
#         tmp_vec <- rbind(tmp_vec, most)
#         tmp_vec_l <- rbind(tmp_vec_l, least)
#         print(tmp_vec)
#         
#       }
#       
#       if (sum(ses_idx == sessions_to_use) != 0){
#         most_likely <- append(most_likely, list(c(colMeans(tmp_vec), ses_idx)))
#         least_likely <- append(least_likely, list(c(colMeans(tmp_vec_l), ses_idx)))
#         
#       }
#       
#       
#       
#     }
#   }
#   
#   df_most <- do.call(rbind, most_likely)
#   df_least <- do.call(rbind, least_likely)
#   
#   melted_dfs <- list()
#   for (fdf in list(df_most[,1:4], df_least[,1:4])) {
#     colnames(fdf) <- c("SI", "Time bins", "Peaks", "Spatial bins")
#     fdf_2 <- fdf
#     #fdf_2 <- fdf[,c("SI", "Time bins", "Spatial bins")]
#     melted_df <- melt(fdf_2)
#     colnames(melted_df) <- c("#", "Comparsion", "Pval")
#     
#     melted_df$Ses <- melted_df$`#` %% 4 + 4
#     melted_df$Ses[melted_df$Ses == 4] <- melted_df$Ses[melted_df$Ses == 4] + 4
#     melted_df$Ses <- as.character(melted_df$Ses)
#     
#     melted_dfs <- append(melted_dfs, list(melted_df))
#     
#   }
#   
#   fdf <- rbind(melted_dfs[[1]], melted_dfs[[2]])
#   fdf$Fit = c(rep("Most likely", times=nrow(melted_dfs[[1]])),rep("Least likely", times=nrow(melted_dfs[[2]])))
#   fdf$Comparsion <- factor(fdf$Comparsion, levels=c("Peaks", "SI", "Spatial bins", "Time bins"))
#   g <-
#     ggplot(melted_dfs[[1]], aes(x=Comparsion, y=Pval)) +
#     geom_violin(alpha=0.1, size=0.75, position=position_dodge(), aes(x=Comparsion, y=Pval),
#                 color="darkmagenta", fill="darkmagenta") +
#     geom_jitter(position=position_jitterdodge(0.2), aes(col=`Ses`), size=0.5) +
#     geom_hline(yintercept=0.05, linetype="dashed", size=0.7) +
#     stat_summary(fun=mean, geom="point", size=2.5, color="red", aes(x=Comparsion)) + 
#     theme_light() +
#     scale_color_manual(values=c("gray80", "gray60", "gray40", "gray20"), "Session") +
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           axis.line = element_line(colour = "black"),
#           #legend.position="NA",
#           axis.text.x = element_text(angle = 45),
#           panel.border = element_blank(),
#           panel.background = element_blank())  +
#     ylab(TeX(r'($P_{value}$)'))
#   
#   g2 <- 
#     ggplot(fdf, aes(x=Comparsion, y=Pval, color=Fit, fill=Fit)) +
#     geom_jitter(data=fdf, position=position_jitterdodge(0.2), aes(x=Comparsion, y=Pval,col=`Ses`, group=Fit)) +
#     geom_hline(yintercept=0.05, linetype="dashed", size=1.25) + 
#     geom_boxplot(size=1, alpha=0.1) + 
#     #stat_summary(fun=mean, geom="point", size=5, color="red", aes(x=Comparsion, group=Fit)) + 
#     scale_color_manual(breaks=c("5","6","7","8","Most likely", "Least likely"), 
#                        values=c("gray80", "gray60", "gray40", "gray20", "darkmagenta", "goldenrod")) +
#     scale_fill_manual(breaks=c("Most likely", "Least likely"), 
#                       values=c("darkmagenta", "goldenrod")) +
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           axis.line = element_line(colour = "black"),
#           panel.border = element_blank(),
#           panel.background = element_blank())
#   
#   return(list(g, g2, df_most, df_least))
# }



plot_gen_place_cell_histograms_new <- 
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
           breaks_peaks=8,
           density_p=F,
           ttle_size=8,
           visualization_spline=-1,
           title=F)  {
    
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
                                  n_bins=24,
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
    
    ylab_p = ifelse(density_p, "Density", "Frequency")
    
    max_SI <- max(c(generated_SI[[1]], true_SI[[1]]))
    SI_breaks =seq(0, max_SI, length.out=breaks_SI)
    hgen <- hist(generated_SI[[1]], breaks=SI_breaks, plot=F)$counts
    htrue <- hist(true_SI[[1]], breaks=SI_breaks, plot=F)$counts
    
    if (visualization_spline != -1) {
      hgen <- smooth_spline_vector(hgen, visualization_spline)
      htrue <- smooth_spline_vector(htrue, visualization_spline)
    }
    
    
    SI_df <- data.frame(counts_gen=hgen,
                        counts_true=htrue,
                        SI=SI_breaks[2:breaks_SI])
    
    if (density_p) {
      SI_df$counts_gen <- SI_df$counts_gen / sum(SI_df$counts_gen)
      SI_df$counts_true <- SI_df$counts_true / sum(SI_df$counts_true)
      
    }
    gSI <- 
      ggplot(SI_df, aes(y=counts_gen, x=SI)) + 
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
            plot.title=element_text(size=ttle_size),
            panel.border = element_blank(),
            panel.background = element_blank()) +
      xlab("SI (bit/spikes)") +
      ylab(ylab_p) +
      geom_line(aes(y=counts_true), size=1, linetype="longdash")
    
    if (title){
      gSI <- gSI+
        ggtitle(sprintf("St=%.2f; Pval=%.4f",
                        pval_func(generated_SI[[1]], true_SI[[1]])$statistic,
                        pval_func(generated_SI[[1]], true_SI[[1]])$p.value))
    }
    
    
    # Compare time bins distribution distances
    max_time_bins <- max(generated_time_bins, true_time_bins_per_cells)
    time_bins_breaks = seq(0, max_time_bins, length.out=breaks_tb)
    hgen <- hist(generated_time_bins, breaks=time_bins_breaks, plot=F)$counts
    htrue <- hist(true_time_bins_per_cells, breaks=time_bins_breaks, plot=F)$counts
    
    if (visualization_spline != -1) {
      hgen <- smooth_spline_vector(hgen, visualization_spline)
      htrue <- smooth_spline_vector(htrue, visualization_spline)
    }
    
    time_bins_df <- data.frame(counts_gen=hgen,
                               counts_true=htrue,
                               Tbins=time_bins_breaks[2:breaks_tb])
    
    if (density_p) {
      time_bins_df$counts_gen <- time_bins_df$counts_gen / sum(time_bins_df$counts_gen)
      time_bins_df$counts_true <- time_bins_df$counts_true / sum(time_bins_df$counts_true)
      
    }
    
    gtimebins <- 
      ggplot(time_bins_df, aes(y=counts_gen, x=Tbins)) + 
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
            plot.title=element_text(size=ttle_size),
            panel.border = element_blank(),
            panel.background = element_blank()) +
      xlab("Time bins (#)") +
      ylab(ylab_p) +
      geom_line(aes(y=counts_true), size=1, linetype="longdash")
    
    if (title){
      gtimebins <- gtimebins+
        ggtitle(sprintf("St=%.2f; Pval=%.4f",
                        pval_func(generated_time_bins, true_time_bins_per_cells)$statistic,
                        pval_func(generated_time_bins, true_time_bins_per_cells)$p.value)) 
    }
    
    
    # Compare active spatial bins distribution distances
    max_spatial_bins <- max(generated_active_spatial_bins_per_cell, true_active_spatial_bins_per_cell)
    spat_bins_breaks =seq(0, max_spatial_bins, length.out=breaks_sb)
    hgen <- hist(generated_active_spatial_bins_per_cell, breaks=spat_bins_breaks, plot=F)$counts
    htrue <- hist(true_active_spatial_bins_per_cell, breaks=spat_bins_breaks, plot=F)$counts
    
    if (visualization_spline != -1) {
      hgen <- smooth_spline_vector(hgen, visualization_spline)
      htrue <- smooth_spline_vector(htrue, visualization_spline)
    }
    
    spat_bins_df <- data.frame(counts_gen=hgen,
                               counts_true=htrue,
                               Sbins=spat_bins_breaks[2:breaks_sb])
    
    if (density_p) {
      spat_bins_df$counts_gen <- spat_bins_df$counts_gen / sum(spat_bins_df$counts_gen)
      spat_bins_df$counts_true <- spat_bins_df$counts_true / sum(spat_bins_df$counts_true)
    }
    
    gspatbins <- 
      ggplot(spat_bins_df, aes(y=counts_gen, x=Sbins)) + 
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
            plot.title=element_text(size=ttle_size),
            panel.border = element_blank(),
            panel.background = element_blank()) +
      xlab("Spatial bins (#)") +
      ylab(ylab_p) + 
      geom_line(aes(y=counts_true), size=1, linetype="longdash")
    
    if (title){
      gspatbins <- gspatbins+
        ggtitle(sprintf("St=%.2f; Pval=%.4f",
                        pval_func(generated_active_spatial_bins_per_cell, true_active_spatial_bins_per_cell)$statistic,
                        pval_func(generated_active_spatial_bins_per_cell, true_active_spatial_bins_per_cell)$p.value)) 
    }
    
    # Compare peaks bins distribution distances
    max_peaks <- max(generated_peaks, true_peaks)
    peaks_breaks = seq(0, max_peaks, length.out=breaks_peaks)
    #peaks_breaks = seq(0, 12, by=2)
    print(generated_peaks)
    print(peaks_breaks)
    print(max_peaks)
    hgen <- hist(generated_peaks, breaks=peaks_breaks, plot=F)$counts
    htrue <- hist(true_peaks,      breaks=peaks_breaks, plot=F)$counts
    
    if (visualization_spline != -1) {
      hgen <- smooth_spline_vector(hgen, visualization_spline)
      htrue <- smooth_spline_vector(htrue, visualization_spline)
    }
    
    peaks_df <- data.frame(counts_gen=c(0,hgen),
                           counts_true=c(0,htrue),
                           Peaks=c(0.5,peaks_breaks[2:breaks_peaks]))
    
    if (density_p) {
      peaks_df$counts_gen <- peaks_df$counts_gen / sum(peaks_df$counts_gen)
      peaks_df$counts_true <- peaks_df$counts_true / sum(peaks_df$counts_true)
    }
    
    
    
    gpeaks <- 
      ggplot(peaks_df, aes(y=counts_gen, x=Peaks)) + 
      geom_bar(stat="identity", 
               color="black", 
               position="identity", 
               fill="cornflowerblue",
               alpha=0.8,
               width=1) + 
      theme_light() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position="NA",
            plot.title=element_text(size=ttle_size),
            panel.border = element_blank(),
            panel.background = element_blank()) +
      xlab("Peaks (#)") +
      ylab(ylab_p) + 
      #xlim(c(0,max_peaks)) +
      #scale_x_continuous(breaks = floor(seq(0, max_peaks, length.out=))) +
      geom_line(aes(y=counts_true), size=1, linetype="longdash")
    
    if (title){
      gpeaks <- gpeaks+
        ggtitle(sprintf("St=%.2f; Pval=%.4f",
                        pval_func(generated_peaks, true_peaks)$statistic,
                        pval_func(generated_peaks, true_peaks)$p.value))  
    }
    
    
    return(list(gSI, gtimebins, gspatbins, gpeaks))
    
  }









figure_5 <- function(paths,
                     estimated_df_folder_cyclic="KS_simulations_likelihood_c",
                     estimated_df_folder_rand="KS_simulations_likelihood_r")  {
  write_path <- sprintf("%s\\figure_5\\",figures_path)
  dir.create(write_path)
  
  true_data_paths <- sprintf("%s\\%s", paths, "equalized")
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
  
  
  estimated_res_cyclic <- get_estimated_df(all_df, true_data_paths, sessions_to_use, estimated_df_folder=estimated_df_folder_cyclic)
  estimated_res_rand <- get_estimated_df(all_df, true_data_paths, sessions_to_use, estimated_df_folder=estimated_df_folder_rand)
  
  
  
  
  both_dir_cyc <- get_both_dir_df(estimated_df=estimated_res_cyclic$estimated_df, 
                                  all_df=all_df)
  both_dir_rand <- get_both_dir_df(estimated_df=estimated_res_rand$estimated_df, 
                                   all_df=all_df)
  
  
  
  sessions <- unique(estimated_res_rand$estimated_df$Session)
  directions <- unique(estimated_res_rand$estimated_df$Direction)
  mic <- unique(estimated_res_rand$estimated_df$Mice)
  
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
  
  gest <- 
    ggplot(mean_sd_df, aes(x=Direction, y=Estimated)) + 
    geom_boxplot(data=both_dir_df, aes(group=interaction(Direction, Subfield, Shuffle), fill=Subfield, color=Subfield), width=0.5, size=1) + 
    geom_jitter(data=both_dir_df,
                aes(x=Direction, y=Estimated, group=interaction(Direction, Subfield, Shuffle)), position=position_jitter(0.2), color="gray20", size=0.75, alpha=0.1) + 
    scale_color_manual(values=c(adjustcolor("gray65", alpha=1), 
                                adjustcolor("royalblue4", alpha=1))) + 
    
    scale_fill_manual(values=c(adjustcolor("gray80", alpha=0.8), 
                               adjustcolor("royalblue2", alpha=0.8))) +
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          legend.position = "NA",
          panel.background = element_blank()) + 
    xlab("Running direction") +
    ylab(ylabs[2]) + 
    ylim(0,1.05)
  
  
  ca1_b <- both_dir_df[both_dir_df$Subfield=="CA1",]
  ca3_b <- both_dir_df[both_dir_df$Subfield=="CA3",]
  # t.test(ca3_b[ca3_b$Direction == "Right","Estimated"], ca1_b[ca1_b$Direction == "Right","Estimated"]) 
  # t.test(ca3_b[ca3_b$Direction == "Both","Estimated"], ca1_b[ca1_b$Direction == "Both","Estimated"]) 
  # t.test(ca3_b[ca3_b$Direction == "Left","Estimated"], ca1_b[ca1_b$Direction == "Left","Estimated"]) 
  
  
  session_mean_sd_df <- session_mean_sd_df[session_mean_sd_df$Direction == "Both",]
  both_dir_df <- both_dir_df[both_dir_df$Direction == "Both",]
  
  # session_mean_sd_df$Session <- as.character(as.numeric(session_mean_sd_df$Session) %% 8)
  # session_mean_sd_df$Session[session_mean_sd_df$Session == "0"] <- "8"
  
  
  g_cyc_rand_ca1_ca3 <- 
    ggplot(session_mean_sd_df, aes(x=Session, y=Estimated)) + 
    geom_point(data=both_dir_df, aes(x=Session, 
                                     y=Estimated, 
                                     group=interaction(Subfield, Shuffle), 
                                     color=interaction(Subfield, Shuffle)), 
               size=3.5, 
               alpha=0.2, 
               position=position_dodge(0.4)) +
    geom_line(aes(group=interaction(Subfield, Shuffle), 
                  color=interaction(Subfield, Shuffle)), 
              size=1.5, 
              position=position_dodge(0.4)) +
    
    geom_errorbar(aes(ymin=Estimated - Sd, ymax=Estimated + Sd, 
                      group=interaction(Subfield, Shuffle), 
                      color=interaction(Subfield, Shuffle)), 
                  size=2, 
                  width=0.3, 
                  position=position_dodge(0.4)) +
    geom_point(aes(group=interaction(Subfield, Shuffle), 
                   color=interaction(Subfield, Shuffle)), 
               size=5, 
               position=position_dodge(0.4)) + 
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position="NA") + 
    geom_hline(yintercept = 1,
               linetype="dashed",
               color="black") +
    xlab("Session") +
    ylab(ylabs[df_idx]) + 
    ylim(0,1.1)
  
  
  
  
  delta_both_dir <- get_both_dir_df(estimated_df=delta_df, 
                                    all_df=all_df)
  
  
  session_mean_sd_df <- delta_both_dir$session_mean_sd_df
  both_dir_df <- delta_both_dir$both_dir_df
  session_mean_sd_df <- session_mean_sd_df[session_mean_sd_df$Direction == "Both",]
  both_dir_df <- both_dir_df[both_dir_df$Direction == "Both",]
  
  # session_mean_sd_df$Session <- as.character(as.numeric(session_mean_sd_df$Session) %% 8)
  # session_mean_sd_df$Session[session_mean_sd_df$Session == "0"] <- "8"
  
  
  gdelta <- 
    ggplot(session_mean_sd_df, aes(x=Session, y=Estimated)) + 
    # geom_point(data=both_dir_df, aes(x=Session, 
    #                                  y=Estimated, 
    #                                  group=interaction(Subfield), 
    #                                  color=interaction(Subfield)), 
    #            size=3.5, 
    #            alpha=0.2, 
    #            position=position_dodge(0.4)) +
    geom_line(aes(group=interaction(Subfield), 
                  color=interaction(Subfield)), 
              size=1.5, 
              position=position_dodge(0.4)) +
    
    geom_errorbar(aes(ymin=Estimated - Sd, ymax=Estimated + Sd, 
                      group=interaction(Subfield), 
                      color=interaction(Subfield)), 
                  size=2, 
                  width=0.3, 
                  position=position_dodge(0.4)) +
    geom_point(aes(group=interaction(Subfield),
                   color=interaction(Subfield)),
               size=5,
               position=position_dodge(0.4)) +
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position="NA") + 
    xlab("Session") +
    ylab("Difference in estimated place cell fraction (%)") + 
    ylim(0,0.4)
  
  
  cost_df <- get_predicted_params_cost_df(new_paths[1:18], sessions_to_use = c(1:2,9:10))
  
  cdf <- 
    data.frame(Y=c(cost_df$all_min, (cost_df$all_min_1 + cost_df$all_min_09) / 2),
               X=c(rep("Estimated percentage", times=len(cost_df$all_min)), 
                   rep("0.9-1.0", times=len(cost_df$all_min_1))))
  
  pdf(file=sprintf("%s\\Estimated_vs_100_cost.pdf", write_path),
      height=(2.5),
      width=(2.5))
  plot(my_boxplot(cdf, xl="MOS place fraction", yl="Minimal cost"))
  dev.off()  
  
  gall <- arrangeGrob(g_cyc_rand_ca1_ca3, gdelta, nrow=1)
  
  pdf(file=sprintf("%s\\delta_shuffles.pdf", write_path),
      height=(2.5),
      width=(5))
  plot(gall)
  dev.off()  
  
  
  pdf(file=sprintf("%s\\diff_box.pdf", write_path),
      height=(3.5),
      width=(5))
  plot(my_boxplot(mdf, yl="Delta measured place cell fration (Cyclic - Random %)", "", xaxis_mod = element_text(angle = 45, vjust = 1, hjust=1)))
  dev.off()  
  
  likdf_tmp <- likelihood_df
  likelihood_df <- (t(apply(likelihood_df, 1, function(r) {return(r/sum(r))})))
  likelihood_df <- as.data.frame(likelihood_df)
  likelihood_sessions <- as.numeric(session_df[,1])
  likelihood_sessions[likelihood_sessions > 8] <-  likelihood_sessions[likelihood_sessions > 8] - 8
  likelihood_df <- cbind(likelihood_df, likelihood_sessions)
  colnames(likelihood_df) <- c(as.character(seq(0.5,1,by=0.1)), "sessions")
  
  
  plots_list <- list()
  for (sf in c("CA1", "CA3")) { 
    
    op_likelihood_df <- likelihood_df[which(session_df[,4] == sf),]
    print(dim(op_likelihood_df))
    
    melted_likelihood <- melt(op_likelihood_df, id.vars="sessions")
    colnames(melted_likelihood) <- c("Session", "Simulated percent", "Normalized likelihood")
    
    likelihood_mean_sd_df <- ddply(melted_likelihood, .(`Simulated percent`), 
                                   function(sub_df) {
                                     ddply(sub_df, .(`Session`), 
                                           function(sub_df_inner) {
                                             mli <- sub_df_inner[,"Normalized likelihood"]
                                             return(c(mean(mli), sem(mli)))})
                                   })
    
    colnames(likelihood_mean_sd_df) <- c("Simulated percent", "Session", "Mean", "Sd")
    
    melted_likelihood$`Simulated percent` <- as.numeric(as.character(melted_likelihood$`Simulated percent`))
    melted_likelihood$`Simulated percent` <- factor(melted_likelihood$`Simulated percent` * 100,
                                                    as.character(seq(50,100, by=10)))
    
    likelihood_mean_sd_df$`Simulated percent` <- as.numeric(as.character(likelihood_mean_sd_df$`Simulated percent`))
    likelihood_mean_sd_df$`Simulated percent` <- factor(likelihood_mean_sd_df$`Simulated percent` * 100, 
                                                        as.character(seq(50,100, by=10)))
    fdf <- likelihood_mean_sd_df[likelihood_mean_sd_df$Session %in% c(1:2,7:8),]
    glikelihood <- 
      ggplot(fdf, aes(x=`Simulated percent`, y=Mean, group=Session)) + 
      geom_bar(stat="summary", aes(group=`Session`, fill=Session), color="NA", 
               position = position_dodge(), width=0.5) + 
      geom_errorbar(aes(ymin=Mean, ymax=Mean + Sd, group=`Session`),
                    width=0.5, position=position_dodge()) +
      # geom_jitter(data=melted_likelihood,
      #             aes(x=`Simulated percent`, y=`Normalized likelihood`, group=`Session`, fill=Session), position=position_jitterdodge(0), color="gray40", size=3, alpha=0.5) + 
      scale_fill_distiller(palette="RdYlBu") +
      theme_light() +     
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = "NA") + 
      xlab("Simulated percent (%)") +
      ylab("Normalized likelihood") + 
      ylim(0,max(likelihood_mean_sd_df[,"Mean"]) + 0.02) +
      ggtitle(sf)
    plots_list <- append(plots_list, list(glikelihood))
    
  }
  
  plots_list$nrow = 1
  glike_all <-   do.call(grid.arrange, plots_list)
  
  pdf(file=sprintf("%s\\likelihood_row.pdf", write_path),
      height=(2.5),
      width=(5))
  plot(glike_all)
  dev.off()  
  
}



get_all_df <- function(true_data_paths, 
                       sessions_to_use, 
                       file_name="percent_df_cyclic_SI_eq_prior.R", 
                       multiple_pvals=F,
                       verbose=F) {
  
  
  unique_paths_both <- 
    unique(unlist(lapply(str_split(true_data_paths, "Matrices"), function(sr) {return(sr[[1]])})))
  
  percent_df_list <- lapply(unique_paths_both, 
                            function(p) {
                              if (verbose) {
                                print(sprintf("Loading %s\\%s", 
                                              p,
                                              file_name))
                              }
                              final_df_list = loadRData(sprintf("%s\\%s", p, file_name))
                              return(final_df_list)
                            })
  
  
  if (multiple_pvals) {
    iterations <- names(percent_df_list[[1]][[16]])
    all_df <- list()
    
    for (iter in iterations) {
      all_df[[as.character(iter)]] <- c()
    }
  } else {
    iterations = .5
    all_df <- c()
  }
  
  
  
  for (iter in iterations) {
    for(i in 1:len(percent_df_list)) {
      measured_percent_df <- percent_df_list[[i]]
      measurements_df <- c()
      
      if (multiple_pvals) {
        cols_to_iterate <- 1:ncol(measured_percent_df[[1]][[iter]])
      } else {
        cols_to_iterate <- 1:ncol(measured_percent_df[[1]])
      }
      
      
      for (col_idx in cols_to_iterate) {
        
        
        # some dataframes were mistakenly saved in an aggregation
        # Meaning instead of saving 5 reps of 16 sessions,
        # Each session was saved and then the aggregated session 
        # Etc 1, 1-2, 1-3 ... 1-16, X 5 times instead of just 1-16 X 5 times
        # Thus, is this is the case
        if (len(measured_percent_df) > 5) {
          full_indices <- seq(16,len(measured_percent_df), by=16)
        } else {
          full_indices <- 1:5
        }
        
        col_sessions <- 
          do.call(rbind,
                  lapply(full_indices,
                         function(idx) {
                           if (multiple_pvals) {
                            df_rep <- measured_percent_df[[idx]][[iter]]
                           } else {
                             df_rep <- measured_percent_df[[idx]]
                           }
                           return(df_rep[sessions_to_use, col_idx])
                         }))
        
        measurements_df <- cbind(measurements_df,
                                 colMeans(col_sessions))
      }
      
      measurements_df <- cbind(measurements_df,
                               sessions_to_use)
      
      measurements_df <- as.data.frame(measurements_df)
      
      split_str <- ifelse(grepl("CA3", unique_paths_both[i]), "CA3\\\\", "CA1\\\\")
      subfield <- ifelse(grepl("CA3", unique_paths_both[i]), "CA3", "CA1")
      measurements_df <- cbind(measurements_df,
                               rep(str_split(str_split(unique_paths_both[[i]], split_str)[[1]][2], "\\\\")[[1]][1],
                                   times=len(sessions_to_use)))
      measurements_df <- cbind(measurements_df, rep(subfield, times=len(sessions_to_use)))
      
      if (multiple_pvals) {
        all_df[[iter]] <- rbind(all_df[[iter]],
                        measurements_df)        
      } else {
      all_df <- rbind(all_df,
                      measurements_df)
      }
  
      }
    
    if (multiple_pvals) {
        colnames(all_df[[iter]]) <- c(colnames(measured_percent_df[[1]][[iter]]),
                                      "Session", "Mice", "Subfield")
    } else {
        colnames(all_df) <- c(colnames(measured_percent_df[[1]]),
                              "Session", "Mice", "Subfield")
    }
  }
  
  if (multiple_pvals) {
    names(all_df) <- iterations
  }
  
  return(all_df)
  
}

get_tunning_corr_df <- function(true_data_paths, 
                                sessions_to_use, 
                                ext="",
                                central_tendency=mean,
                                verbose=F,
                                all_cells=F,
                                file_name="non_place_tunning_corr_within.Rda",
                                extract_corr_func=function(df) {central_tendency(colMeans(df[,colSums(df) != 0]))}) {
  
  tuning_df <- data.frame(Corr=c(),
                          Session=c(),
                          Direction=c(),
                          Mice=c(),
                          Subfield=c())
  
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
      
      if (len((grep("tunning_corr_within.Rda", list.files(sprintf("%s%s", tpath, ext))))) == 0) {
        print(tpath)
        missing_paths <- c(missing_paths, tpath)
        print("Missing")
        next
      }
      
      if (verbose) {
        print(sprintf("Loading %s", sprintf("%s%s\\%s", tpath, ext, file_name)))
      }
      load(sprintf("%s%s\\%s", tpath, ext, file_name), verbose=verbose)
      
      
      exf <- function(df) {central_tendency(df[!is.na(df)])}
      corr_value <- exf(res)
      
      if (grepl("CA3", working_path)) {
        mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', tpath))
        mice_str <- substr(tpath, mice_str_index, mice_str_index+4) 
      } else {
        mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', tpath))
        mice_str <- substr(tpath, mice_str_index, mice_str_index+3) 
      }
      
      subfield <- ifelse(grepl("CA3", working_path), "CA3", "CA1")

      
      if (all_cells) { 
        tuning_df <- rbind(tuning_df,
                           data.frame(res, 
                                      rep(idx, len(res)),
                                      rep(ifelse(grepl("Left", tpath), "Left", "Right"), len(res)), 
                                      rep(mice_str, len(res)), 
                                      rep(subfield, len(res))))
      } else {

        tuning_df <- rbind(tuning_df,
                           data.frame(corr_value, idx, ifelse(grepl("Left", tpath), "Left", "Right"), mice_str, subfield))

      }
    }
  }
  
  return(tuning_df)
}


get_estimated_df <- function(all_df, 
                             true_data_paths, 
                             sessions_to_use, 
                             estimated_df_folder="KS_simulations_likelihood_c",
                             edge_free=F) {
  estimation_list <- list()
  likelihood_list <- list()
  measured_percent_list <- list()
  estimated_per_percent_list <- list()
  session_metavar_list <- list()
  missing_paths <- c()
  
  all_paths <- c();
  
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
      
      if (len((grep("pct_estimate_df.R", list.files(sprintf("%s\\%s%s", 
                                                            tpath, 
                                                            ifelse(edge_free, "edge_free\\", ""),
                                                            estimated_df_folder))))) == 0) {
        print(tpath)
        missing_paths <- c(missing_paths, tpath)
        print("Missing")
        next
      }
      
      
      
      all_paths <- c(all_paths, tpath)
      
      #final_df <- get_fit_params_figures_from_path(tpath, just_likelihood = T)
      load(sprintf("%s\\%s%s\\%s", 
                   tpath, 
                   ifelse(edge_free, "edge_free\\", ""),
                   estimated_df_folder, 
                   "pct_estimate_df.R"))
      #print(colnames(final_df))
      
      measured_p <- -1
      if ("measured_pct" %in% colnames(final_df)) {
        print("Got measured percent")
        measured_p <- final_df[1,"measured_pct"]
        final_df <- final_df[,-ncol(final_df)]
      }
      
      
      likelihood <- final_df[,ncol(final_df)] # Put last column on different vector
      final_df <- final_df[,-ncol(final_df)] # remove last column
      
      if (len(likelihood) == 7) {
        likelihood = likelihood[-1]
        final_df <- final_df[-1,]
      }
      
      # In case we saved it in reverse (mistake in some datasets)
      # that 100% were saved as 50%
      if (mean(final_df[1,]) > mean(final_df[nrow(final_df),])) {
        print("Reversing")
        final_df_t <- apply(final_df, 2, rev)
        likelihood_t <- rev(likelihood)
        
        names(likelihood_t) <- names(likelihood)
        rownames(final_df_t) <- rownames(final_df)
        
        final_df <- final_df_t
        likelihood <- likelihood_t
      }
      
      if (measured_p == -1) {
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
        
        measured_p <- as.numeric(names(which.max(percents_tab)))
      }
      
      measured_percent_list <- append(measured_percent_list, list(measured_p))
      estimation_list <- append(estimation_list, list(as.numeric(names(which.max(likelihood)))))
      likelihood_list <- append(likelihood_list, list(likelihood))
      estimated_per_percent_list <- append(estimated_per_percent_list, list(rowMeans(final_df)))
      
      split_str <- ifelse(grepl("CA3", working_path), "CA3\\\\", "CA1\\\\")
      subfield <- ifelse(grepl("CA3", working_path), "CA3", "CA1")
      
      #mice_str <- str_split(str_split(tpath, split_str)[[1]][2], "\\\\Matrices")[[1]]
      
      if (grepl("CA3", working_path)) {
        mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', tpath))
        mice_str <- substr(tpath, mice_str_index, mice_str_index+4) 
      } else {
        mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', tpath))
        mice_str <- substr(tpath, mice_str_index, mice_str_index+3) 
      }
      
      session_metavar_list <- append(session_metavar_list, list(c(idx, mice_str, ifelse(grepl("Left", tpath), "Left", "Right"), subfield)))
      
    }
  }
  
  likelihood_df <- do.call(rbind, likelihood_list)
  estimated_per_percent_df <- do.call(rbind, estimated_per_percent_list)
  session_df <- do.call(rbind, session_metavar_list)
  
  estimation_vec <- unlist(estimation_list)
  measured_vec <- unlist(measured_percent_list)
  
  estimated_df <- data.frame(percent_estimated=estimation_vec)
  estimated_df <- cbind(estimated_df, session_df)
  colnames(estimated_df) <- c("Estimated", "Session", "Mice", "Direction", "Subfield")
  
  
  return(list(estimation_list=estimation_list,
              likelihood_list=likelihood_list,
              measured_percent_list=measured_percent_list,
              estimated_per_percent_list=estimated_per_percent_list,
              session_metavar_list=session_metavar_list,
              missing_paths=missing_paths, 
              likelihood_df=likelihood_df,
              estimated_per_percent_df=estimated_per_percent_df,
              session_df=session_df,
              estimation_vec=estimation_vec,
              measured_vec=measured_vec,
              estimated_df=estimated_df,
              all_paths=all_paths))
}

get_both_dir_df <- function(all_df, estimated_df, dev_func=sd) {
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
                                           c(weighted_estimate_both_dir, ses, mice, "Both", left[1,"Subfield"]))
                         }
                         
                         return(sub_df)
                       })
  
  
  both_dir_df_all <- ddply(all_df, .(Mice), 
                       function(sub_df) {
                         mice = sub_df[1,"Mice"]
                         subfield <- sub_df[1, "Subfield"]
                         
                         fdf <- data.frame()
                         for (ses in sub_df$Session) {
                           
                           ses_df  <- sub_df[sub_df$Session == ses,]
                           right_p <- ses_df[,"Right active"]
                           left_p <- ses_df[,"Left active"]
                           
                           
                           weighted_pc_out_all <- 
                             (ses_df["Right_out_of_all"] * right_p + 
                              ses_df["Left_out_of_all"] * left_p) / (left_p + right_p)
                           
                           weighted_pc_out_active <- 
                             (ses_df["Right"] * right_p + 
                              ses_df["Left"] * left_p) / (left_p + right_p)
                           
                           weighted_active <- 
                              (right_p + left_p) / 2
                           
                           fdf <- rbind(fdf,
                                        data.frame(as.numeric(weighted_pc_out_all),
                                          as.numeric(weighted_pc_out_active),
                                          as.numeric(weighted_active), 
                                          ses, mice, "Both", subfield))
                         }
                         
                         
                         
                         fdf <- as.data.frame(fdf)
                         
                         colnames(fdf) <- c("Of_All", "Of_Active", "Active", "Session",
                                            "Mice", "Direction", "Subfield")
                         return(fdf)
                       })
  # 
  # estimated_per_percent_df <- cbind(estimated_per_percent_df, measured_vec)
  # estimated_per_percent_df <- as.data.frame(estimated_per_percent_df)
  # estimated_per_percent_df <- cbind(estimated_per_percent_df, session_df)
  
  
  
  # melted_likelihood <- melt(t(apply(likelihood_df[,-1], 1, function(r) {return(r/sum(r))})))
  # colnames(melted_likelihood) <- c("#", "Simulated percent", "Normalized likelihood")
  
  both_dir_df <- as.data.frame(both_dir_df)
  both_dir_df[,"Estimated"] <- as.numeric(both_dir_df[,"Estimated"])
  
  
  
  both_dir_df$TwoEnvSession <- both_dir_df$Session
  both_dir_df$Session <- as.character(as.numeric(both_dir_df$Session) %% 8)
  both_dir_df$Session[both_dir_df$Session == "0"] <- "8"
  
  both_dir_df_all$Session <- as.character(as.numeric(both_dir_df_all$Session) %% 8)
  both_dir_df_all$Session[both_dir_df_all$Session == "0"] <- "8"
  
  mean_sd_df <- ddply(both_dir_df, .(Subfield), 
                      function(subfield_df) { 
                        return(ddply(subfield_df, .(Direction), 
                                     function(sub_df) {return(c(mean(sub_df[,"Estimated"]),dev_func(sub_df[,"Estimated"])))}))
                      });
  
  session_mean_sd_df <- 
    ddply(both_dir_df, .(Subfield), 
          function(subfield_df) {
            
            return(ddply(subfield_df, .(Session), 
                         function(sub_df) {
                           sd_df <- ddply(sub_df, .(Direction), 
                                          function(sub_df_inner) 
                                          { est <- as.numeric(sub_df_inner[,"Estimated"])
                                          return(c(mean(est), dev_func(est)))})
                           return(sd_df)
                         }));
          });
  
  
  session_mean_sd_df_all <- 
    ddply(both_dir_df_all, .(Subfield), 
          function(subfield_df) {
            
            return(ddply(subfield_df, .(Session), 
                         function(sub_df) {
                           sd_df <- ddply(sub_df, .(Direction), 
                                          function(sub_df_inner) 
                                          { of_all <- as.numeric(sub_df_inner[,"Of_All"])
                                            of_active <- as.numeric(sub_df_inner[,"Of_Active"])
                                            active <- as.numeric(sub_df_inner[,"Active"])
                                          return(c(mean(of_all), dev_func(of_all),
                                                   mean(of_active), dev_func(of_active),
                                                   mean(active), dev_func(active)))
                                          })
                           return(sd_df)
                         }));
          });
  
  colnames(mean_sd_df) <- c("Subfield", "Direction", "Estimated", "Sd")
  colnames(session_mean_sd_df) <- c("Subfield", "Session", "Direction", "Estimated", "Sd")
  colnames(session_mean_sd_df_all) <- c("Subfield", "Session", "Direction", 
                                        "Of_All_mean", "Of_All_sd",
                                        "Of_Active_mean", "Of_Active_sd",
                                        "Active_mean", "Active_sd")
  
  return(list(both_dir_df=both_dir_df,
              mean_sd_df=mean_sd_df,
              session_mean_sd_df=session_mean_sd_df,
              both_dir_df_all=both_dir_df_all,
              session_mean_sd_df_all=session_mean_sd_df_all))
}



get_predicted_params_cost_df <- function(paths, sessions_to_use=c(1:16),range=seq(0.5,1,by=0.1), prediction_path="KS_simulations_estimate", estimated_df_folder="KS_simulations_likelihood_c") {
  
  true_data_paths <- sprintf("%s\\%s", paths, "equalized")
  
  
  all_df <- get_all_df(true_data_paths, sessions_to_use)
  estimated <- get_estimated_df(all_df, true_data_paths, sessions_to_use, estimated_df_folder=estimated_df_folder_cyclic)  
  
  
  all_min_09 <- c()
  all_min_1 <- c()
  all_mean <- c()
  all_min <- c()
  
  for (idx in 1:len(estimated$all_paths)) {
    curr_path <- estimated$all_paths[idx] 
    
    estimated_percent <- estimated$estimation_vec[idx]
    load(sprintf("%s\\%s\\param_fit_%.3f\\pred_df.R", curr_path, prediction_path, estimated_percent))
    predicted <- t(predicted_df[,5:24])
    
    
    all_min <- c(all_min, min(colMeans(predicted))) 
    all_mean <- c(all_mean, mean(colMeans(predicted))) 
    
    load(sprintf("%s\\%s\\param_fit_%.3f\\pred_df.R", curr_path, prediction_path, 1))
    predicted_1 <- t(predicted_df[,5:24])
    all_min_1 <- c(all_min_1, min(colMeans(predicted_1))) 
    
    load(sprintf("%s\\%s\\param_fit_%.3f\\pred_df.R", curr_path, prediction_path, 0.9))
    predicted_09 <- t(predicted_df[,5:24])
    all_min_09 <- c(all_min_09, min(colMeans(predicted_09))) 
    
  }
  
  return(list(all_min=all_min,
              all_min_1=all_min_1,
              all_mean=all_mean,
              all_min_09=all_min_09))
  
}




my_boxplot <- function(df, xl, yl, xaxis_mod=element_text()) {
  
  gbox <- 
    ggplot(df, aes(x=X, y=Y)) + 
    geom_jitter(data=df, aes(x=X, y=Y), position=position_jitter(0.2), color="gray20", size=0.75, alpha=0.1) + 
    geom_boxplot(data=df, width=0.5, size=1, 
                 color=adjustcolor("gray65", alpha=1), 
                 fill=adjustcolor("gray80", alpha=0.8)) +
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          legend.position = "NA",
          panel.background = element_blank(),
          axis.text.x=xaxis_mod) + 
    xlab(xl) +
    ylab(yl) 
  
  return(gbox)
}

figure_6_tuning_corr <- function() {
 
  true_data_paths <- sprintf("%s\\%s", all_data_paths[1:18], "equalized")
  ct_func = mean
  extraction_func = function(df) {central_tendency(df[!is.na(df)])}
  file_name = "place_tunning_corr_within.Rda"
  tuning_df <- get_tunning_corr_df(true_data_paths, 
                                   1:16, 
                                   ext="\\fit_simulation", 
                                   central_tendency = ct_func, 
                                   extract_corr_func = extraction_func, 
                                   file_name=file_name,
                                   all_cells=F)
  
  
  colnames(tuning_df) <- c("Corr", "Session", "Direction",  "Mice", "Subfield")
  
  tuning_df$Session <- as.character(as.numeric(tuning_df$Session) %% 8)
  tuning_df$Session[tuning_df$Session == "0"] <- "8"
  
  min_corr_1 <- min(tuning_df$Corr)
  max_corr_1 <- max(tuning_df$Corr)
  
  session_mean_sd_df <- 
    ddply(tuning_df, .(Subfield), 
          function(subfield_df) {
            
            return(ddply(subfield_df, .(Session), 
                         function(sub_df) {
                           corr_vec <- as.numeric(sub_df[,"Corr"])

                           res <- c(mean(corr_vec, na.rm=T), sem(corr_vec))
                           return(res)
                         }));
          });
  
  colnames(session_mean_sd_df) <- c("Subfield", "Session",  "Corr", "Sd")
  
  gtuning <- 
    ggplot(session_mean_sd_df, aes(x=Session, y=Corr)) + 
    geom_point(data=tuning_df, aes(x=Session, y=Corr, group=Subfield, color=Subfield), 
               size=.5, alpha=0.2, position=position_dodge(0.4)) +
    geom_line(aes(group=interaction(c(Subfield)), 
                  color=interaction(c(Subfield))), size=1.5, position=position_dodge(0.4)) +
    
    geom_errorbar(aes(ymin=Corr - Sd, ymax=Corr + Sd, group=Subfield, color=Subfield), 
                  size=2, width=0.3, position=position_dodge(0.4)) +
    geom_point(aes(group=Subfield, color=Subfield), size=5, position=position_dodge(0.4)) + 
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "NA") + 
    xlab("Session") + 
    ggtitle("Simulated")
  
  
  tuning_df <- get_tunning_corr_df(true_data_paths, 
                                   1:16, 
                                   ext="", 
                                   central_tendency = ct_func, 
                                   extract_corr_func = extraction_func, 
                                   file_name=file_name,
                                   all_cells = F)
  
  
  colnames(tuning_df) <- c("Corr", "Session", "Direction",  "Mice", "Subfield")
  
  tuning_df$Session <- as.character(as.numeric(tuning_df$Session) %% 8)
  tuning_df$Session[tuning_df$Session == "0"] <- "8"
  
  min_corr_2 <- min(tuning_df$Corr)
  max_corr_2 <- max(tuning_df$Corr)
  
  session_mean_sd_df <- 
    ddply(tuning_df, .(Subfield), 
          function(subfield_df) {
            
            return(ddply(subfield_df, .(Session), 
                         function(sub_df) {
                           corr_vec <- as.numeric(sub_df[,"Corr"])
                           res <- c(mean(corr_vec, na.rm=T), sem(corr_vec))
                           return(res)
                         }));
          });
  
  colnames(session_mean_sd_df) <- c("Subfield", "Session",  "Corr", "Sd")
  
  gtuning2 <- 
    ggplot(session_mean_sd_df, aes(x=Session, y=Corr)) + 
    geom_point(data=tuning_df, aes(x=Session, y=Corr, group=Subfield, color=Subfield), 
               size=.5, alpha=0.2, position=position_dodge(0.4)) +
    geom_line(aes(group=interaction(c(Subfield)), 
                  color=interaction(c(Subfield))), size=1.5, position=position_dodge(0.4)) +
    
    geom_errorbar(aes(ymin=Corr - Sd, ymax=Corr + Sd, group=Subfield, color=Subfield), 
                  size=2, width=0.3, position=position_dodge(0.4)) +
    geom_point(aes(group=Subfield, color=Subfield), size=5, position=position_dodge(0.4)) + 
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "NA") + 
    xlab("Session") +
    ggtitle("Real")
  
  grid.arrange(gtuning + ylim(c(min(1.1 * c(min_corr_1, min_corr_2)), max(1.1 * c(max_corr_1, max_corr_2)))), 
               gtuning2 + ylim(c(min(1.1 * c(min_corr_1, min_corr_2)), max(1.1 * c(max_corr_1, max_corr_2)))), nrow=1)
  
  

  
}


color_plot <- function(g, x_axes) {
  gf <-  g + 
    scale_fill_manual(values=
                        c(adjustcolor(c("royalblue2", "gray80"), alpha=0.8),
                          adjustcolor(c("royalblue2", "gray80"), alpha=0.8),
                          adjustcolor(c("royalblue2", "gray80"), alpha=0.8)),
                      breaks=c("Both.CA3", "Both.CA1",
                               "Left.CA3", "Left.CA1",
                               "Right.CA3", "Right.CA1")) +
    scale_color_manual(values=
                         c("royalblue3","gray65",
                                       "royalblue3","gray65",
                                       "royalblue3","gray65"),
                                       breaks=c("Both.CA3", "Both.CA1",
                                                "Left.CA3", "Left.CA1",
                                                "Right.CA3", "Right.CA1")) + 
    base_plot_theme + 
    theme(text=element_text(size=15)) +
    xlab(x_axes) +
    ylab("Fraction") + 
    ylim(0,1.05)
  
  return(gf)
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
        tpath_f <- sprintf("%s/%s", tpath, fit_sim)
        tpath2 <- sprintf("%s%s%s%s.Rda",
                          tpath_f,
                          prefix,
                          "place_cells_by_sample_duration_cor",
                          ext)
        tpath4 <- sprintf("%s/%s", tpath3, "duration_pval_matrices\\")
        
        if (len((grep(sprintf("%splace_cells_by_sample_duration_cor%s.Rda", prefix, ext), list.files(tpath_f)))) == 0) {
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
        
        load(tpath2, verbose=T)
        
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
                         
                         # df <- df_list[[i]]
                         # new_df <- data.frame(bins=seq(1, 20, length.out = 8),
                         #                      rand=df[,1],
                         #                      cyclic=df[,2])
                         # 
                         # rownames(new_df) <- seq(1, 20, length.out = 8)
                         
                         df <- df_list[[i]]
                         new_df <- data.frame(bins=seq(1, 20, length.out = 15),
                                              pct=df[,1])
                         
                         rownames(new_df) <- seq(1, 20, length.out = 15)
                         return(new_df) 
                       })
    
    
    merged_rand = df_list2[[1]][,c("bins", "rand")]
    merged_cyclic = df_list2[[1]][,c("bins", "cyclic")]
    
    
    merged_cor = df_list2[[1]][,c("bins", "pct")]
    
    # Merge dataframes
    for (i in 2:len(df_list2)) {
      # merged_rand <- merge(merged_rand,
      #                      df_list2[[i]][,c("bins", "rand")],
      #                      by="bins",
      #                      all=T)
      # 
      # merged_cyclic <- merge(merged_cyclic,
      #                        df_list2[[i]][,c("bins", "cyclic")],
      #                        by="bins",
      #                        all=T)        
      
      merged_cor <- merge(merged_cor,df_list2[[i]][,c("bins", "pct")],
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
    
    
    
    
    means_cor <- rowMeans(merged_cor[,-1], na.rm=T)
    means_rand <- rowMeans(merged_rand[,-1], na.rm=T)
    sem_rand <- apply(merged_rand[,-1], 1, sd, na.rm=T)
    sem_cor <- apply(merged_cor[,-1], 1, sd, na.rm=T)
    #labels_rand <- sapply(merged_rand[,1], function(b) {as.numeric(str_split(str_split(b, "\\(")[[1]][2], ",")[[1]][1])})
    labels_rand <- merged_rand[,1]
    
    means_cyc <- rowMeans(merged_cyclic[,-1], na.rm=T)
    sem_cyc <- apply(merged_cyclic[,-1], 1, sd, na.rm=T)
    #labels_cyc <- sapply( merged_cyclic[,1], function(b) {as.numeric(str_split(str_split(b, "\\(")[[1]][2], ",")[[1]][1])})
    labels_cyc <- merged_cyclic[,1]
    
    labels_cor <- merged_cor[,1]
    
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
    
    df_cor <- data.frame(fraction=means_cor,
                         time=as.numeric(labels_cor),
                         sd=sem_cor)
    
    plot_df <- df_cor
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
    
    
    plot_df <- df_cor
    
    gcor <- 
      ggplot(plot_df, aes(x=time, y=fraction)) + 
      geom_ribbon(aes(ymin=fraction-sd, ymax=fraction+sd), color="NA", alpha=0.1, size=1) +
      geom_line(size=1.2) +
      theme_light() +
      ylim(c(0.25, 0.7))+
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



plot_pval_by_sample_duration_figure <- 
  function(working_path_list, 
           ses_ind=c(5:8, 13:16), 
           ext="",
           fit_sim="",
           prefix="",
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
        tpath2 <- sprintf("%s\\%s%s%s%s.Rda",
                          tpath3,
                          fit_sim,
                          prefix,
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
                         
                         f_vec <- t(apply(all_cells, 1, function(r) {r}))
                         #f_vec <- all_cells
                         f_vec <- f_vec[!is.nan(f_vec[,6]),]
                         f_vec <- apply(f_vec, 2, median)
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


plot_neur_raster <- function(stim_trace, spike_train, idx, coltoc="red", cx=1) {
  plot(stim_trace, type="l", col="gray60")
  points(which(spike_train[idx,] != 0 ), stim_trace[which(spike_train[idx,] != 0 )], col=coltoc,
         pch=19, cex=cx)
}

