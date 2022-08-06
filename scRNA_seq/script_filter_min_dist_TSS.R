
#Script to find peaks with genes closest to TSS.

df <- data.frame(
A=c("A", "A", "A", "B", "B", "B", "C", "C", "C"),
x=c(5, 10, 2, 20, 13, 40, 25, 55, 35),
y=rnorm(9)
)

library(plyr)
ddply(df, .(A), function(z) {
    z[z$x == min(z$x), ][1, ]
})


