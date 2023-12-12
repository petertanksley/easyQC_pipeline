#=====Install packages=========================================================
# library(pacman)
# p_load(forestplot,
#        plotrix)
# install.packages('https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_23.8.tar.gz', repos = NULL, type = 'source')


#=====Load packages============================================================
library(pacman)
p_load(EasyQC,
       data.table,
       forestplot,
       plotrix)

EasyQC("../temp/ecf_files/easyqc_ext2.0_afr_alcp.ecf")


