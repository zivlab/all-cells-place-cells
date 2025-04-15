library(ggrepel)

figure_1_literature_scatter <- function() {
wpath <- sprintf("%s//place_cell_fractions_literature_v2.csv", base_path)
df <- read.csv(wpath)

write_path <- sprintf("%s\\figure_1\\",figures_path)
dir.create(write_path)

sizes=c(big=3.5,
        big_1=3,
        medium=2.75,
        medium_1=2.5,
        medium_2=2.25,
        small=2)

# gfractions <- 
# ggplot(df, aes(x=Index, y=Fraction, label=Paper)) + 
#   geom_point(aes(shape=interaction(Animal, Subfield, sep=" "), color=Recording), size=2) +
#   geom_text_repel(hjust=0.1) +   theme_classic() + ylim(c(0.2,1)) + 
#   theme(legend.title=element_blank(), axis.text.x=element_blank()) + 
#   xlab("Year") + ylab("Reported fraction of place cells")


xaxis_labels <- 
  str_replace(str_replace(str_replace(as.character(levels(cut(df$Year, c(1988, df$Year[c(5,10,15)], 2023)))), ",", "-"), "]", ""), "\\(", "")
names(xaxis_labels) <- seq(5,20, by=5)
# gfractions_no_leg <- 
#   ggplot() + 
#   geom_point(data=df, aes(x=Index, 
#                          y=Fraction,
#                          shape=interaction(Animal, sep=" "), 
#                          color=Recording), size=2) +
#   geom_text_repel(hjust=0.1, data=df, aes(x=Index, y=Fraction, label=Paper)) +   theme_classic() + ylim(c(0.2,1)) + 
#   xlab("Year") + ylab("Reported fraction of place cells") +
#   scale_x_discrete(labels=xaxis_labels) + 
#   theme(legend.title=element_blank(), 
#         legend.position = "NA", 
#         text=element_text(size=15))

gfrac_no_legend <- 
ggplot(data=df, aes(x=factor(Index), 
                    y=Fraction, 
                    label=Paper)) + 
  geom_point(aes(shape=Animal, 
                  color=Recording), 
             size=2) +
  geom_text_repel(hjust=0.1) +   
  theme_classic() + 
  ylim(c(0.2,1)) + 
  xlab("Year") + 
  ylab("Reported fraction of place cells") + 
  scale_x_discrete(breaks=c(5,10,15,20), 
                   labels=xaxis_labels) + 
  scale_color_manual(values=c("#F05A28", "#652D90")) + 
  theme(legend.title=element_blank(), 
        legend.position = "NA", 
        text=element_text(size=15))
  


g_frac <- 
  ggplot(data=df, aes(x=Year,#factor(Index), 
                      y=Fraction, 
                      label=Paper)) + 
  geom_point(aes(shape=Animal, 
                 color=Recording), 
             size=2) +
  geom_text_repel(hjust=0.1) +   
  theme_classic() + 
  ylim(c(0.2,1)) + 
  xlab("Year published") + 
  ylab("Reported fraction of place cells") + 
  #scale_x_discrete(breaks=c(5,10,15,20), 
   #                labels=xaxis_labels) + 
  scale_color_manual(values=c("#F05A28", "#652D90")) + 
  theme(legend.title=element_blank(), 
        text=element_text(size=15, color="black"),
        axis.text=element_text(color="black"),
        legend.position="None")
    


for (size_name in names(sizes)) {
  size = sizes[[size_name]]
  dir.create(sprintf("%s\\%s", write_path, size_name))
  pdf(file=sprintf("%s\\%s\\new_literature_fractions_no_legend.pdf",
                   write_path,
                   size_name),
      height=size,
      width=size * 1.9)
  
  plot(gfrac_no_legend)
  dev.off()
  
  pdf(file=sprintf("%s\\%s\\new_literature_fractions.pdf",
                   write_path,
                   size_name),
      height=size,
      width=size * 2)
  
  plot(g_frac)
  dev.off()
}

}




figure_1_npc_ill <- function() {
  write_path <- sprintf("%s\\figure_1\\",figures_path)
  dir.create(write_path)
  
  sizes=c(big=3.5,
          big_1=3,
          medium=2.75,
          medium_1=2.5,
          medium_2=2.25,
          small=2)

npc <- c(rlnorm(500000, sd=0.6, mean=-1.5),rlnorm(500000, sd=0.2, mean=0.15))
dens <- density(npc)
dens2 <- density(npc, n =100)
min_points <- which(diff(as.numeric(diff(dens$y) < 0)) == -1)
min_points_2 <- which(diff(as.numeric(diff(dens2$y) < 0)) == -1)

dens_df <- data.frame(x=dens2$x, y=dens2$y)

npc_df <- dens_df[1:min_points_2[1],]
npc_df <- rbind(npc_df, c(npc_df[nrow(npc_df),1], 0))
pc_df <- dens_df[(min_points_2[1]):nrow(dens_df),]
pc_df <- rbind(pc_df, c(pc_df[nrow(pc_df),1], 0))
pc_df <- rbind(c(pc_df[1,1], 0), pc_df)


# gbimodal <- 
# ggplot(dens_df, aes(x=x,y=y)) + 
#   theme_classic() + 
#   geom_polygon(data=npc_df, aes(x=x,y=y), fill="red") +
#   geom_polygon(data=pc_df, aes(x=x,y=y), fill="royalblue4") + 
#   geom_vline(xintercept=dens_df[min_points[1],"x"], 
#            linetype="dashed",
#            linewidth=1) + 
#   geom_bar(linewidth=1, stat=stat_count()) + 
#   xlab("Spatial Information") + 
#   xlim(c(min(npc_df),3)) +
#   ylab("Density")

gbimodal2 <- 
ggplot(dens_df, aes(x=x,y=y)) + 
  theme_classic() + 
  #geom_bar(stat="identity") +
  geom_polygon(data=npc_df, aes(x=x,y=y), fill="#CE3736", alpha=0.9) +
  geom_polygon(data=pc_df, aes(x=x,y=y), fill="#3F5DAB", alpha=0.9) + 
  geom_vline(xintercept=dens_df[min_points_2[1],"x"], 
             linetype="dashed",
             linewidth=1) + 
  xlab("Spatial Information") + 
  xlim(c(min(npc_df),2.4)) +
  ylab("Density") + 
  theme(text=element_text(size=15)) + 
  scale_y_continuous(expand=c(0,0))

true_pc_dens <- rlnorm(50000, mean=0.2, sd=0.2)

true_pc_df <- data.frame(x=true_pc_dens)#[which(true_pc_dens < 2.5)])




gtrue <- 
ggplot(true_pc_df, aes(x=x)) + 
  geom_density(aes(x=x), fill="#02A0A0", alpha=0.8, color=NA) + 
  xlab("Spatial Information") + 
  xlim(c(.4,2.5)) +
  ylab("") + 
  xlab("") + 
  theme_classic() + 
  theme(text=element_text(size=12)) +
  scale_y_continuous(expand=c(0,0))

for (size_name in names(sizes)){
  size <- sizes[[size_name]]
  pdf(file=sprintf("%s\\%s\\bimodal.pdf",
                   write_path,
                   size_name),
      height=size,
      width=size)
  
  plot(gbimodal2)
  dev.off()
  
  pdf(file=sprintf("%s\\%s\\true.pdf",
                   write_path,
                   size_name),
      height=size * .4,
      width=size * .4)
  
  plot(gtrue)
  dev.off()  

  }

}



figure_1_place_cell_calc_demo <- function()
{

a <- get_spike_train_and_stim_trace_from_path(all_data_paths[3],7)
spike_train <- a[[1]]
stim_trace <- a[[2]]
fr <- rowMeans(spike_train) / dt
processed_real <- preprocess_spike_train(spike_train, stim_trace)

true_cells_spike_train <- processed_real$working_cells_spike_train
true_firing_rate <- processed_real$working_firing_rate
true_time_bins_per_cells <- processed_real$working_time_bins_per_cells

#params <- get_fit_params(path, estimation_path="simulations_2_new_new", likelihood_path = "simulations_3_new", pct_range=seq(0.5,1,by=0.1))
params <- get_fit_params(sprintf("%s\\equalized\\session_%d", all_data_paths[3], 7))

print(params)
simulated_tuning_curve <-
  generate_tuning_curves_cost(n = nrow(spike_train),
                              percentage = params$pct,
                              average_width = params$params["average_width"], 
                              sd_width = params$params["sd_width"],
                              fixed_fr=fr,
                              noise=params$params["noise"],
                              double_peak_pct = params$params["double_peak_pct"],
                              n_bins=24,
                              plot=F)

pois_factor <- currate_spike_train_cost(simulated_tuning_curve, 
                                        true_time_bins_per_cells,
                                        stim_trace,
                                        verbose = T,
                                        jump=0.05)

generated_spike_train <- 
  generate_spike_trains_cost(tuning_curves = simulated_tuning_curve,
                             stim_trace = stim_trace,
                             factor=pois_factor,
                             fs=1)

spike_train <- generated_spike_train

fr <- rowMeans(spike_train) / dt
processed_real <- preprocess_spike_train(spike_train, stim_trace)

true_cells_spike_train <- processed_real$working_cells_spike_train
true_firing_rate <- processed_real$working_firing_rate
true_time_bins_per_cells <- processed_real$working_time_bins_per_cells

idx <- processed_real$ind[sample(order(true_firing_rate, decreasing=T)[1], 1)]
firing_ind <- which(spike_train[idx,] > 0)

run_df <- data.frame(Time=((1:len(stim_trace) / len(stim_trace)) * 20), 
                     Position=stim_trace * 4) 
run_df$Position[which(run_df$Position == 4)] <- 0



shuffle <- generate_shuffled_spike_train(spike_train = spike_train,
                                         idx = idx,
                                         num_of_shuffles = 1000,
                                         shuffle_type = "c")


# spike train in the size of Nshuffles x T 
shuffled_spike_train <- t(array_reshape(shuffle, c(ncol(spike_train),1000)))
shuffled_idx <- sample(1:1000, 1)
shuffled_firing_ind <- which(shuffled_spike_train[shuffled_idx,] > 0)
rm(shuffle)


cell_df <- data.frame(Positions=stim_trace[firing_ind] * 4,
                      Times=run_df$Time[firing_ind],
                      Shuffled_times=run_df$Time[shuffled_firing_ind],
                      Shuffled_positions=stim_trace[shuffled_firing_ind] * 4)

shuffled_spike_train <- rbind(spike_train[idx,],
                              shuffled_spike_train)

shuffled_rate <- rowMeans(shuffled_spike_train) / dt
tmp <- compute_tuning(shuffled_spike_train, stim_trace)
shuffled_stim_prob <- tmp[[1]]
shuffled_tuning <- tmp[[2]]

# Sanity check
# assert(all(shuffled_stim_prob == stim_prob))
# assert(rowSums())

# Compute SI for all shuffles
shuffled_SI <- compute_SI(shuffled_stim_prob,
                          shuffled_tuning,
                          shuffled_rate)



SI <- shuffled_SI[[1]][2:1001]


real_SI <- shuffled_SI[[1]][1]
tuning_curve <- (unlist(lapply(sort(unique(stim_trace)), function (sb) {mean(spike_train[idx,which(stim_trace == sb)]) / dt}))) 
barp <- barplot(tuning_curve)
smoothed_tuning <- smooth.spline(barp, tuning_curve, all.knots = T, lambda=1e-4)
smoothed_tuning$y[smoothed_tuning$y < 0] <- 0
pos <- sort(unique(stim_trace)) * 4; pos[pos == 4] <- 0
smoothed_df <- data.frame(FR=smoothed_tuning$y,
                          Positions=pos)
grundf <-
  ggplot() + 
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
  geom_point(data=cell_df, aes(x=Times,y=Positions),
             fill="NA",
             color="#DA1C5C",
             stroke=0,
             alpha=0.8,
             size=2) + 
  coord_flip() +
  ggtitle(sprintf("SI: %f", real_SI))


grundf_shuffled <-
  ggplot() + 
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
  geom_point(data=cell_df, aes(x=Shuffled_times,y=Shuffled_positions),
             fill="NA",
             color="#DA1C5C",
             stroke=0,
             alpha=0.8,
             size=2) + 
  coord_flip() +
  ggtitle(sprintf("SI: %f", SI[shuffled_idx]))



si_df <- data.frame(si=SI)

gSI <- 
  ggplot(si_df, aes(x=si)) + 
  geom_histogram(bins=100, fill="#6D91CB", color="black", size=.05) +
  theme_classic() +
  xlab("Spatial Information") + 
  ylab("Count") +
  geom_vline(xintercept=,
             linetype="dashed",
             linewidth=1) + 
  ggtitle("Shuffles distribution")
  



}



gaus <- function(x, mn=0, sd=1) {(1/(sd * sqrt(2*pi))) * exp((-.5)*(((x - mn)**2)/(sd ** 2)))}
mdf <- data.frame(x=1:800, g2=gaus(x_reig, sd=.6, mn=1.7), g1=gaus(x_reig))

#gexamp <- 
ggplot(mdf) + 
  geom_line(aes(x=x,y=g1+g2)) + 
  geom_vline(xintercept=466) + 
  geom_line(aes(x=x,y=g1), linetype="dashed") + 
  geom_line(aes(x=x,y=g2), linetype="dashed") + 
  base_plot_theme + scale_y_discrete(expand=c(0,0))



for (size_name in names(sizes)){
  size <- sizes[[size_name]]
  pdf(file=sprintf("%s\\%s\\example_var.pdf",
                   write_path,
                   size_name),
      height=size,
      width=size)
  
  plot(gexamp)
  dev.off()

  
}
