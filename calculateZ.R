library(edf)
library(plyr)
library(ggplot2)

calculate_Zcore <- function(input_data, time_window=1, threshold=7, sample_name="EEG_test", output_dir="/Users/fanyingtang/Downloads/", channel="EEG1") {
  
  # time_window is in secs (1sec)
  # input_data is EEG file
  # default cutoff is 7
  
  edf_file <- read.edf(input_data)
  EEG_signal <- edf_file$signal[[channel]]$data
  
  # binsize the signal, using samples over 1sec as a bin
  # FFT: "signal over time" into "frequency distribution"
  HZ <- edf_file$header.signal[[channel]]$samplingrate #number of sampling per second
  bin=HZ*time_window
  sample_num <- length(EEG_signal)/bin
  threshold_true=threshold*time_window
  EEG_signal_split <- split(EEG_signal, ceiling(seq_along(EEG_signal)/bin))
  
  # Do fast fourier transform
  EEG_signal_fft <- lapply(EEG_signal_split, fft)
  EEG_signal_fft_abs=lapply(EEG_signal_fft, function(x) abs(x))
  
  # use 7 as the threshold. Determine fast or slow fraction. Calculate zscore=(slow-fast)/(slow+fast)
  # zscore > 0 sleep; zscore < 0 awake
  slow_fast = lapply(EEG_signal_fft_abs, function(x) return(c(slow=sum(x[1:(threshold_true+1)]),fast=sum(x[(threshold_true+2):(bin)]))))
  Zscore = lapply(slow_fast, function(x) (x[1]-x[2])/(x[1]+x[2]))
  Zscore_vec <- as.numeric(Zscore)
  
  # get the zscore plot. output zscore and FFT results in csv files
  time_axis <- seq(1,sample_num)
  lm(Zscore_vec ~ time_axis)
  pdf(paste0(paste0(output_dir,sample_name,".pdf")))
  plot(time_axis, Zscore_vec, pch=16, cex=0.3, xlab="Time (secs)", ylab="Zscore")
  lm(Zscore_vec ~ time_axis)
  dev.off()
  write.table(Zscore_vec, paste0(output_dir, sample_name, ".csv"), quote=F)
  write.table(as.data.frame(EEG_signal_fft), paste0(output_dir, sample_name, "_FFT.csv"), quote=F,sep=",",row.names=F)
}

calculate_Zcore(input_data = "/Users/fanyingtang/Downloads/EEG_sample_data.edf")
