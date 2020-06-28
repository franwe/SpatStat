build.convex.grid <- function (x, y, npts) {
   library(splancs) 
      ch <- chull(x, y) # index for pts on convex hull 
   ch <- c(ch, ch[1]) # Add first point on to the end.
   border <- cbind(x[ch], y[ch])  # This works as a splancs poly
   # Now fill it with grid points
   xy.grid <- gridpts(border, npts)
   return(xy.grid)
}