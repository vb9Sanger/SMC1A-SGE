install.packages("beeswarm")

library(beeswarm)
library(ggplot2)

color3 <- exon3_pilot_LFC[["pch"]]

#plot for exon3 - all clinVar
beeswarm(
  # Plot the petal lengths as a function of the species of iris flower
  exon3_pilot_LFC[["LFC"]] ~ exon3_pilot_LFC[["Variant"]],
  # Set the x- and y-axis labels
  xlab = "Variant Classification", ylab = "LFC (Day15/Day4)",
  # Set the labels of the groups
  labels = c("B/LB", "VUS", "P/LP"),
  pwcol = color3,
  pch = 19,
  ylim= c(-3.5,3.5), cex.lab =1.5, cex.axis = 1.5)
abline(h=0, lty=2)
legend("topright", legend = c("Synonymous", "Missense", "LoF"), title = "Consequence", pch = 19, col = c("chartreuse4", "blue3", "red4"))

#######################

color4 <- exon4_pilot_LFC[["pch"]]

#plot for exon4 - all ClinVar
beeswarm(
  # Plot the petal lengths as a function of the species of iris flower
  exon4_pilot_LFC[["LFC"]] ~ exon4_pilot_LFC[["Variant"]],
  # Set the x- and y-axis labels
  xlab = "Variant Classification", ylab = "LFC (Day15/Day4)",
  # Set the labels of the groups
  labels = c("B/LB", "VUS", "P/LP"),
  pwcol = color4,
  pch = 19,
  #col = c("chartreuse4", "blue3", "red4")
  ylim= c(-3.5,3.5), cex.lab =1.5, cex.axis = 1.5)
abline(h=0, lty=2)
legend("topright", legend = c("Synonymous", "Missense", "LoF"), title = "Consequence", pch = 19, col = c("chartreuse4", "blue3", "red4"))
