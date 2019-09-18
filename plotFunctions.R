
plot_ecdfErr = function (x, main = NULL, sub = NULL, xlab = deparse(substitute(x)),
 ax = seq(0,1,0.001),colLine = "gray59",colErr = "gray67",addAll = F)
{

    force(xlab)
    if (is.null(main))
        main <- paste0("ecdf(", deparse(substitute(x)), ") + 95% K.S. bands")

    ci.col = colErr

    n <- length(x)
    if (is.null(sub))
        sub <- paste("n = ", n)
    ec <- ecdf(x)
    xx <- get("x", envir = environment(ec))
    yy <- get("y", envir = environment(ec))
    D <- KSd(n)
    yyu <- pmin(yy + D, 1)
    yyl <- pmax(yy - D, 0)
    ecu <- stepfun(xx, c(yyu, 1))
    ecl <- stepfun(xx, c(yyl, yyl[n]))
    if(!addAll)
      plot(ax,ec(ax), main = main, sub = sub, xlab = xlab,type ="l",col = colLine,ylim = c(0,1),xlim = c(0,1),lwd = 2)
    else
      lines(ax,ec(ax),lwd = 2.5)

    plot(ecu, add = TRUE, verticals = TRUE, do.points = FALSE,
        col.hor = ci.col, col.vert = ci.col,lty = 4)
    plot(ecl, add = TRUE, verticals = TRUE, do.points = FALSE,
        col.hor = ci.col, col.vert = ci.col,lty = 4)
}
