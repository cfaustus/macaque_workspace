#####
## SI plotting scenarios
## 2018 june 07
#####

library(scales)

scenarios = read.csv('output/bayesian_scenarios.csv', header = TRUE)

par(las = 1)
hist(scenarios$t1.0, col= alpha('darkgoldenrod4', 0.4),
     xlim=c(0,1), ylim=c(0,500), border = F,
     breaks=seq(0,0.5,0.01),
     ylab = 'frequency',
     xlab = 'posterior estimates of the probability that a macaque \nhas a variant alpha globin phenotype',
     main = '',
     cex.lab = 1,
     cex.axis = 1)
hist(scenarios$t1.1A, col=alpha('darkgoldenrod4',0.4), add=T, border = T,
     breaks=seq(0,1,0.01))
hist(scenarios$t2.0, col=alpha('firebrick3',0.4), add=T, border = F,
     breaks=seq(0,1,0.01))
hist(scenarios$t2.1A, col=alpha('firebrick3',0.4), add=T, border = T,
     breaks=seq(0,1,0.01))

legend(0.3, 500,
       legend = c('presence of virulent malaria unlikely','presence of virulent malaria highly likely'),
       pch = 19,
       bty='n',
       col= c('darkgoldenrod4','firebrick3'))
legend(0.3, 420,
       legend = c('0.02 cutoff (main text)','0.04 cutoff'),
       pch = 22,
       bty = 'n',
       col= c('lightgray','black'))


par(las = 1)
hist(scenarios$t1.2A, col= alpha('darkgoldenrod4', 0.4),
     xlim=c(0,1), ylim=c(0,500), border = F,
     breaks=seq(0,0.5,0.01),
     ylab = 'frequency',
     xlab = 'posterior estimates of the probability that a macaque \nhas a variant alpha globin phenotype',
     main = '',
     cex.lab = 1,
     cex.axis = 1)
hist(scenarios$t1.2B, col=alpha('darkgoldenrod4',0.4), add=T, border = T,
     breaks=seq(0,1,0.01))
hist(scenarios$t2.2A, col=alpha('firebrick3',0.4), add=T, border = F,
     breaks=seq(0,1,0.01))
hist(scenarios$t2.2B, col=alpha('firebrick3',0.4), add=T, border = T,
     breaks=seq(0,1,0.01))

legend(0.3, 500,
       legend = c('presence of virulent malaria unlikely','presence of virulent malaria highly likely'),
       pch = 19,
       bty='n',
       col= c('darkgoldenrod4','firebrick3'))
legend(0.3, 420,
       legend = c('unknowns nonvirulent','unknowns virulent'),
       pch = 22,
       bty = 'n',
       col= c('lightgray','black'))


par(las = 1)
hist(scenarios$t1.0, col= alpha('darkgoldenrod4', 0.4),
     xlim=c(0,1), ylim=c(0,500), border = F,
     breaks=seq(0,0.5,0.01),
     ylab = 'frequency',
     xlab = 'posterior estimates of the probability that a macaque \nhas a variant alpha globin phenotype',
     main = '',
     cex.lab = 1,
     cex.axis = 1)
hist(scenarios$t1.3, col=alpha('darkgoldenrod4',0.4), add=T, border = T,
     breaks=seq(0,1,0.01))
hist(scenarios$t2.0, col=alpha('firebrick3',0.4), add=T, border = F,
     breaks=seq(0,1,0.01))
hist(scenarios$t2.3, col=alpha('firebrick3',0.4), add=T, border = T,
     breaks=seq(0,1,0.01))

legend(0.3, 500,
       legend = c('presence of virulent malaria unlikely','presence of virulent malaria highly likely'),
       pch = 19,
       bty='n',
       col= c('darkgoldenrod4','firebrick3'))
legend(0.3, 420,
       legend = c('main text','malaria proxy of neighbor'),
       pch = 22,
       bty = 'n',
       col= c('lightgray','black'))

par(las = 1)
hist(scenarios$t1.0, col= alpha('darkgoldenrod4', 0.4),
     xlim=c(0,1), ylim=c(0,500), border = F,
     breaks=seq(0,0.5,0.01),
     ylab = 'frequency',
     xlab = 'posterior estimates of the probability that a macaque \nhas a variant alpha globin phenotype',
     main = '',
     cex.lab = 1,
     cex.axis = 1)
hist(scenarios$t1.4, col=alpha('darkgoldenrod4',0.4), add=T, border = T,
     breaks=seq(0,1,0.01))
hist(scenarios$t2.0, col=alpha('firebrick3',0.4), add=T, border = F,
     breaks=seq(0,1,0.01))
hist(scenarios$t2.4, col=alpha('firebrick3',0.4), add=T, border = T,
     breaks=seq(0,1,0.01))

legend(0.3, 500,
       legend = c('presence of virulent malaria unlikely',
                  'presence of virulent malaria highly likely'),
       pch = 19,
       bty='n',
       col= c('darkgoldenrod4','firebrick3'))
legend(0.3, 420,
       legend = c('only major variants','including minor variants'),
       pch = 22,
       bty = 'n',
       col= c('lightgray','black'))

