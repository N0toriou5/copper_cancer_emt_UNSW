# Modified function to separate labels on scatter plot
library(wordcloud)
textplot3<-function (x, y, words, cex = 1, pch = 16, pointcolor = "#FFFFFF00", 
                     new = FALSE, show.lines = TRUE, line.col="black",rstep=0.5,...) 
{
    if (new) {
        plot(x, y, type = "n", ...)
    }
    lay <- wordlayout(x, y, words, cex, rstep=rstep, ...)
    if (show.lines) {
        for (i in seq_len(length(x))) {
            xl <- lay[i, 1]
            yl <- lay[i, 2]
            w <- lay[i, 3]
            h <- lay[i, 4]
            if (x[i] < xl || x[i] > xl + w || y[i] < yl || y[i] > 
                yl + h) {
                points(x[i], y[i], pch = pch, col = pointcolor, 
                       cex = 0.5)
                nx <- xl + 0.5 * w
                ny <- yl + 0.5 * h
                lines(c(x[i], nx), c(y[i], ny), col = line.col)
            }
        }
    }
    text(lay[, 1] + 0.5 * lay[, 3], lay[, 2] + 0.5 * lay[, 4], 
         words, cex = cex, ...)
}