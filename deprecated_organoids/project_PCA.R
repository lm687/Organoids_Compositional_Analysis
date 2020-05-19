## 5 samples and 3 variables
x <- MCMCpack::rdirichlet(5, c(1,1,1,1))
clr <- function(x) log(x)-mean(log(x))
x_clr <- t(apply(x, 1, clr))
stopifnot(all(dim(x) == dim(x_clr)))
princomp_x <- princomp(x_clr)
princomp_x$loadings
princomp_x$scores

princomp_x$center
princomp_x$scale

apply(x_clr, 2, mean)
princomp_x$center
apply(sweep(x_clr, 2, princomp_x$center, '+'), 2, mean)

apply(sweep(x_clr, 2, princomp_x$center, '-'), 2, mean)

##################################################

par(mfrow=c(1,1))
dev.off()
x_clr_centered <- x_clr
x_clr_centered <- sweep(x_clr, 2, princomp_x$center, '-')
apply(x_clr_centered, 2, mean)

projection <- ( x_clr_centered %*% (princomp_x$loadings))

projection <- projection[1,]
projection_to_plot <- t(projection[1:2])
plot(princomp_x$scores[,1:2], pch=19, col=c('red', 'black', 'black', 'black', 'black'),
     xlim=c(min(c(princomp_x$scores[,1], projection_to_plot[,1])),max(c(princomp_x$scores[,1], projection_to_plot[,1]))),
     ylim=c(min(c(princomp_x$scores[,2], projection_to_plot[,2])),max(c(princomp_x$scores[,2], projection_to_plot[,2]))))
points(projection_to_plot, col='blue')

## comp1 is correct but comp2 isn't

#projection_to_plot
## still need to take into account scale

#points(predict(princomp_x, newdata = x_clr[1:2,])[,1:2],
#       col='blue', pch=4)

