tmp = installed.packages()

installedpackages = as.vector(tmp[is.na(tmp[,"Priority"]), 1])
save(installedpackages, file="/home/yongbao/scrna-kulkarni-macrophages/misc/installed_packages.rda")

load("~/Desktop/installed_packages.rda")

for (count in 1:length(installedpackages)) install.packages(installedpackages[count])