#' Example data generator for DA_CV
#'
#' Creates a simple raster stack and random sample points to demonstrate
#' how \code{DA_CV()} can be run. This is only for testing/demo purposes.
#'
#' @return A list with `samples` (sf object), `target` (SpatRaster),
#'   `env_stack` (SpatRaster), and the output of \code{DA_CV()}.
#' @export
#'
#' @examples
#' \dontrun{
#' demo <- example_DA_CV()
#' print(demo$out)
#' }
example_DA_CV <- function() {
	library(sf)
	library(terra)

	# Create a simple raster as "target"
	target <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
	terra::values(target) <- matrix(runif(100, 50, 150), nrow = 10, ncol = 10)

	# Create two environmental predictor rasters
	env1 <- terra::rast(target)
	terra::values(env1) <- runif(terra::ncell(env1))
	env2 <- terra::rast(target)
	terra::values(env2) <- runif(terra::ncell(env2))
	env_stack <- c(target, env1, env2)
	names(env_stack) <- c("response", "predictor_1", "predictor_2") # ensure RF variable names match

	# Generate some random sample points
	set.seed(42)
	pts <- data.frame(
		x = runif(10, 0, 10),
		y = runif(10, 0, 10)
	)
	samples <- sf::st_as_sf(pts, coords = c("x", "y"), crs = terra::crs(target))
	samples <- terra::extract(env_stack, samples, bind = TRUE) |>
		st_as_sf()

	# Run DA_CV
	out <- DA_CV(
		samples = samples,
		predictors = env_stack,
		response = "response",
		folds_k = 5,
		autoc_threshold = 0.2,
		cate_num = 5,
		seed = 10
	)

	return(list(samples = samples, target = target, env_stack = env_stack, out = out$DA_folds))
}
