.First.lib <- function(lib, pkg) {
	if(!require(methods)) {
		stop("'methods' package is required.")
		require(methods)
	}
	packageStartupMessage(sprintf("Package %s (%s) loaded.
To cite, see citation(\"%s\")\n", pkg, packageDescription(pkg)$Version, pkg))
}

setClass(Class = "HWEdata", representation(data = "matrix", size = "numeric", vec = "numeric"))

setMethod("initialize", "HWEdata",
		function(.Object, data = matrix(NA), size = numeric(0), vec = numeric(0)) {
			tmp <- data
			R <- length(tmp)
			r <- (-1 + sqrt(1 + 8*R))/2
			data <- matrix(NA, nrow = r, ncol = r)
			for (i in 1:r) {
				data[i, 1:i] <- tmp[((i - 1)*i/2 + 1):(i*(i + 1)/2)]
			}
			.Object@data <- data
			.Object@size <- r
			.Object@vec <- tmp
			.Object
		}
)

setClass(Class = "HWEintr", representation(bf = "numeric", npp = "numeric", draws = "numeric", data = "matrix"))

setMethod("initialize", "HWEintr",
		function(.Object, bf = numeric(0), npp = numeric(0), draws = numeric(0), data = matrix(NA)) {
			.Object@bf <- bf
			.Object@npp <- npp
			.Object@draws <- draws
			.Object@data <- data
			.Object
		}
)

setMethod("summary",
		"HWEdata",
		function(object, plot = FALSE, ...) {
			if (plot) {
				plot(object)
			}
			out <- list(n = sum(object@vec), size = nrow(object@data), data = object@data)
			out
		}
)

setMethod("plot",
		signature(x = "HWEdata", y = "missing"),
		function(x, ndens = 5, xlab = "x", ylab = "density", nbin = 10, histogram = FALSE, bands = FALSE, confid = .95, start = 1, ...) {
			y <- x@data
			r <- dim(y)[1]
			n <- sum(y, na.rm = TRUE)
			p <- y/n
			lbl <- paste("A", 1:r, sep = "")
			p_a <- (rowSums(y, na.rm = TRUE) + colSums(y, na.rm = TRUE))/(2*n)
			p_hwe <- matrix(NA, nrow = r, ncol =r)
			for (i in 1:r) {
				for (j in 1:i) {
					p_hwe[i, j] <- 2*p_a[i]*p_a[j] 
				}
				p_hwe[i, i] <- p_a[i]^2
			}
			side <- sqrt(p)
			side_hwe <- sqrt(p_hwe)
			old.par <- par(no.readonly = TRUE)
			par(mai = c(.82, .82, 1.02, .42))
			plot(1:r, 1:r, type = "n", xlim = c(0, (r+1)), ylim = c(0, (r+1)), xlab = "", ylab = "First chromosome", main = "",
					lab = c(rep((r + 1), 2), 7), xaxt = "n", yaxt = "n")
			clr <- as.numeric(col2rgb("brown"))
			color <- rgb(clr[1], clr[2], clr[3], alpha=180, maxColorValue=255)
			clr <- as.numeric(col2rgb("orange"))
			color_hwe <- rgb(clr[1], clr[2], clr[3], alpha=180, maxColorValue=255)
			axis(side = 3, at = 1:r, labels = lbl, tick = FALSE)
			axis(side = 2, at = 1:r, labels = rev(lbl), tick = FALSE, las = 1)
			for (i in 1:r) {
				abline(v = (i - .5), lty = 3, col = gray(.5), lwd = .35)
				abline(h = (i - .5), lty = 3, col = gray(.5), lwd = .35)
				for (j in 1:r) {
					rect(xleft = (i - side[j, i]/2), ybottom = ((r - j + 1) - side[j, i]/2), xright = (i + side[j, i]/2), ytop = ((r - j +1) + side[j, i]/2),
							col = color, border = NA, density = NA)
					rect(xleft = (i - side_hwe[j, i]/2), ybottom = ((r - j + 1) - side_hwe[j, i]/2), xright = (i + side_hwe[j, i]/2), ytop = ((r - j +1) +
							side_hwe[j, i]/2), col = color_hwe, border = NA, density = NA)
					text(j, (r - i + 1), y[i, j], cex = .55, col = gray(.4))
				}
				segments(x0 = (i - .5), y0 = ((r - i + 1) + .5), x1 = (i + .5), y1 = ((r - i + 1) + .5), lwd = .5)
				segments(x0 = (i + .5), y0 = ((r - i + 1) + .5), x1 = (i + .5), y1 = ((r - i + 1) - .5), lwd = .5)
			}
			mtext(text = "Second chromosome", side = 3, line = ((par()$mar)[3] - 2.5))
			abline(v = ((r + 1) - .5), lty = 3, col = gray(.5), lwd = .35)
			abline(h = ((r + 1) - .5), lty = 3, col = gray(.5), lwd = .35)
			segments(x0 = .5, y0 = (r + .5), x1 = .5, y1 = .5, lwd = .5)
			segments(x0 = .5, y0 = .5, x1 = (r + .5), y1 = .5, lwd = .5)
			mtext(text = "Brown squares areas: sample genotype proportions", side = 1, line = ((par()$mar)[1] - 3.7), cex = .7)
			mtext(text = "Orange squares areas: sample Hardy-Weinberg proportions", side = 1, line = ((par()$mar)[1] - 2.7), cex = .7)
			mtext(text = "Numbers within squares: sample genotype counts", side = 1, line = ((par()$mar)[1] - 1.7), cex = .7)
			par(old.par)
		}
)

setMethod("summary",
		"HWEintr",
		function(object, plot = FALSE, ...) {
			if (plot) {
				plot(object)
			}
			out <- list(bf = object@bf, npp = object@npp, mc.draws = summary.default(object@draws))
			out
		}
)

setMethod("plot",
		signature(x = "HWEintr", y = "missing"),
		function(x, xlab = "Iteration", ylab = "Monte Carlo Average of Bayes Factor", ...) {
			M <- length(x@draws)
			d <- x@draws 
			res <- cumsum(d)/(1:M)
			sdres <- sqrt(cumsum((d - res)^2))/(1:M)
			plot(1:M, res, xlab = xlab, ylab = ylab, main = "", col = gray(0), type = "n", ylim = mean(res) + 20*c(-sdres[M], sdres[M]))
			abline(h = x@bf, col = gray(.5), lty = "longdash")
			lines(1:M, res, xlab = xlab, ylab = ylab, main = "", col = gray(0), lwd = 2)
			lines(res + 2*sdres, col = "gold", lwd = 2)
			lines(res - 2*sdres, col = "gold", lwd = 2)
		}
)
