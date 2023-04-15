library("related")
setwd("~/Desktop")
input <- readgenotypedata("KBTOGenotypes_Related.txt") #reads in file with all loci and genotypes 
simdata <- familysim(input$freqs , 600) 
errors <- c(0.0408, 0.7303, 0.0140, 0.2215, 0.1099, 0.3269, 0.000, 0.3355, 0.3073, 0.1075, 0.1128, 0.6325, 0.000, 0.2232, 0.0154, 0.6325, 0.4000, 0.0159) #null allele rates calculated from other analyses
output <- coancestry(simdata, error.rates=errors, lynchrd=1, quellergt=1) #running to find point estimates of relatedness for simulated data, incorporating null allele into simulated runs 
simrel <- cleanuprvals(output$relatedness, 600) #reduces to set of 600 pairs per dyad type of interest 
lynchrdpo <- simrel [1:600 , 8] #lines 8-15 are selecting the range of rows and columns that correspond to the appropriate relatedness value and estimator.
lynchrdfs <- simrel[(600 + 1) : (2 * 600), 8]
lynchrdhs <- simrel[((2 * 600) + 1) : (3 * 600), 8] 
lynchrdur <- simrel[((3 * 600) + 1) : (4 * 600), 8]
quellergtpo <- simrel [1:600 , 10]
quellergtfs <- simrel[(600 + 1) : (2 * 600), 10]
quellergths <- simrel[((2 * 600) + 1) : (3 * 600), 10] 
quellergtur <- simrel[((3 * 600) + 1) : (4 * 600), 10]
t.test(lynchrdpo, quellergtpo) #lines 16-19 compare sampling variance between both estimators for a given relationship
t.test(lynchrdfs, quellergtfs)
t.test(lynchrdhs, quellergths)
t.test(lynchrdur, quellergtur)
lynchrd <- rep("LRM", 600) #lines 20-23 create a list of labels for the different estimators, with each repeated the appropriate number of times
quellergt <- rep("QGM", 600)
estimator2 <- c(lynchrd, quellergt)
Estimator <- rep(estimator2 , 2)
po <- rep("Parent-Offspring", (2 * 600))  #lines 24-28 creates a list of labels for the different relatedness types
fs <- rep("Full-Sibs", (2 * 600))
hs <- rep("Half-Sibs", (2 * 600)) 
ur <- rep("Unrelated", (2 * 600))
relationship <- c(po, fs, hs, ur)
relatednesspo <- c(lynchrdpo, quellergtpo) #lines 29-33 combine the different values for each estimator based on relatedness type, as lists.
relatednessfs <- c(lynchrdfs, quellergtfs)
relatednesshs <- c(lynchrdhs, quellergths)
relatednessur <- c(lynchrdur, quellergtur)
Relatedness_Value <- c(relatednesspo, relatednessfs, relatednesshs, relatednessur)
combineddata <- as.data.frame(cbind(Estimator , relationship , Relatedness_Value)) #combines data 
combineddata$Relatedness_Value <- as.numeric(as.character(combineddata$Relatedness_Value)) #plots data 
ggplot(combineddata , aes(x = Estimator , y = Relatedness_Value), ylim = c(-0.5, 1.0)) +
geom_boxplot() + facet_wrap(~ relationship)

urval <- rep(0, 600) #here, we are determining with Pearson's correlation how well expectations for the different order relatedness accords with observation
hsval <- rep (0.25 , 600)
fsval <- rep (0.5 , 600)
poval <- rep (0.5 , 600)
relvals <- c(poval , fsval , hsval , urval)
cor(relvals , simrel[, 8])
cor(relvals , simrel[, 10]) 
write.csv(combineddata, "/Users/jeffreylernercoleman/Desktop/combineddata.csv", row.names=TRUE) #writes combined data frame to CSV 
