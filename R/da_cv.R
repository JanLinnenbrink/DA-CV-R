#' Domain Adaptive Cross-Validation (DA-CV)
#'
#' Implements prediction-domain adaptive cross-validation (DA-CV). This method combines random-domain cross-validation (RDM-CV)
#' for "similar" prediction regions and spatial cross-validation (SP-CV) for
#' "dissimilar" prediction regions, based on dissimilarity quantification via
#' adversarial validation.
#'
#' @param samples An \code{sf} object. Must contain sample locations (POINT geometry)
#'   and attributes including target variable and (optionally) sample IDs.
#' @param target A \code{SpatRaster}. Target raster (e.g., biomass).
#' @param env_stack A \code{SpatRaster}. Stack of environmental predictor rasters.
#' @param nodata_value Numeric. Value in \code{target} that indicates missing data.
#' @param folds_k Integer. Number of folds for cross-validation.
#' @param cate_num Integer. Number of categorical covariates.
#' @param autoc_threshold Numeric. Spatial autocorrelation threshold for SP-CV.
#'
#' @return A list containing a data frame with DA-CV fold assignments, the dissimilarity raster, as well as the threshold used to classify this raster.
#' Additionally writes probability and category rasters, and a CSV with fold assignments.
#'
#' @examples
#' library(sf)
#' library(terra)
#'
#' # Create a simple raster as "target"
#' target <- rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
#' values(target) <- matrix(runif(100, 50, 150), nrow = 10, ncol = 10)
#'
#' # Create two simple environmental predictor rasters
#' env1 <- rast(target); values(env1) <- runif(ncell(env1))
#' env2 <- rast(target); values(env2) <- runif(ncell(env2))
#' env_stack <- c(env1, env2)
#'
#' # Generate some random sample points inside the raster extent
#' set.seed(42)
#' pts <- data.frame(
#'   x = runif(10, 0, 10),
#'   y = runif(10, 0, 10),
#'   ID = 1:10
#' )
#' samples <- st_as_sf(pts, coords = c("x", "y"), crs = crs(target))
#'
#' # Run DA-CV with small folds for demo
#' # (requires RDM_CV() and spatial_plus_cv() to be implemented in the package)
#' out <- DA_CV(
#'   samples = samples,
#'   target = target,
#'   env_stack = env_stack,
#'   nodata_value = NA,
#'   folds_k = 2,
#'   cate_num = 2,
#'   autoc_threshold = 0.3
#' )
#' print(out)
#'
#' @export
DA_CV <- function(samples, target, env_stack, nodata_value, folds_k, cate_num, autoc_threshold) {
	start_time <- Sys.time()

	# ---- Check inputs ----
	stopifnot(inherits(samples, "sf"))
	stopifnot(inherits(target, "SpatRaster"))
	stopifnot(inherits(env_stack, "SpatRaster"))

	row_length <- nrow(target)
	col_length <- ncol(target)

	samples <- samples[, c(names(env_stack), names(target))]

	# ---- Step 1: Dissimilarity quantification ----
	# Find no-data cells in target
	nodata_cells <- which(terra::values(target) == nodata_value)

	# Extract environmental variables at sample locations
	sample_envs <- terra::extract(env_stack, samples, ID = FALSE)
	sample_coords <- sf::st_coordinates(samples)
	n_samples <- nrow(sample_envs)

	# Select random prediction cells (same number as samples), excluding no-data
	all_cells <- seq_len(terra::ncell(target))
	valid_cells <- setdiff(all_cells, nodata_cells)
	set.seed(123)
	pred_cells <- sample(valid_cells, n_samples, replace = FALSE)
	pred_envs <- terra::extract(env_stack, pred_cells)

	# Prepare adversarial validation dataset
	X <- rbind(sample_envs, pred_envs)
	y <- c(rep(1, n_samples), rep(0, n_samples))
	trainDat <- data.frame(y = y, X)
	trainDat <- trainDat[complete.cases(trainDat), ]

	rf <- ranger::ranger(formula = y ~ ., trainDat, num.trees = 500, probability = TRUE)

	pred <- stats::predict(rf, data = X, type = "response")
	probs <- pred$predictions[, 2]

	roc_obj <- pROC::roc(y, probs)
	auc_val <- pROC::auc(roc_obj)
	diss_value <- ifelse(auc_val <= 0.5, 0, round((auc_val - 0.5) * 2, 2))
	threshold <- diss_value * 0.5
	message("Dissimilarity value = ", round(diss_value * 100), "%, threshold = ", threshold)

	# ---- Step 2: Apply classifier to all grid cells ----
	all_envs <- terra::as.data.frame(env_stack, na.rm = FALSE)
	all_envs <- as.data.frame(all_envs) # ensure type matches X
	stopifnot(identical(names(X), names(all_envs))) # sanity check

	all_pred <- stats::predict(rf, all_envs, type = "response")
	all_probs <- all_pred$predictions[, 2]

	prob_raster <- target
	terra::values(prob_raster) <- all_probs

	cate_raster <- prob_raster
	terra::values(cate_raster) <- ifelse(terra::values(prob_raster) >= threshold, 2, 1)

	sim_count <- sum(terra::values(cate_raster) == 2, na.rm = TRUE)
	dissim_count <- sum(terra::values(cate_raster) == 1, na.rm = TRUE)
	sim_ratio <- sim_count / (sim_count + dissim_count)
	dissim_ratio <- 1 - sim_ratio

	samples_category <- terra::extract(cate_raster, terra::vect(samples), ID = FALSE, raw = TRUE)
	samples$category <- samples_category[, 1]

	print(table(samples$category))

	# ---- Step 3: Run CV on subsets ----
	similar_samples <- samples[samples$category == 1, ]
	dissim_samples <- samples[samples$category == 2, ]
	if (nrow(similar_samples) > 0) {
		RDM_folds <- RDM_CV(similar_samples, folds_k)
	}

	if (nrow(dissim_samples) > 0) {
		SP_folds <- spatial_plus_cv(
			samples = dissim_samples,
			response_name = names(target),
			cate_col_start = 0,
			cate_col_end = 0,
			k = folds_k,
			sp_threshold = autoc_threshold,
			method = "SP"
		)
	}

	# ---- Step 4: Combine results ----
	folds_df <- data.frame(
		ID = if ("ID" %in% names(samples)) samples$ID else seq_len(nrow(samples)),
		fold = NA_integer_,
		x = sample_coords[, 1],
		y = sample_coords[, 2],
		category = samples$category
	)

	# Assign folds based on category
	if (nrow(similar_samples) > 0) {
		folds_df$fold[folds_df$category == 1] <- RDM_folds$fold
	}

	if (nrow(dissim_samples) > 0) {
		folds_df$fold[folds_df$category == 2] <- SP_folds$fold
	}

	# Drop the category column if you don’t want it in the final output
	folds_df$category <- NULL

	# Convert back to sf
	DA_folds <- sf::st_as_sf(folds_df, coords = c("x", "y")) |>
		sf::st_set_crs(sf::st_crs(samples))

	runtime <- Sys.time() - start_time
	message("DA-CV completed in ", round(runtime, 1), " seconds.")

	return(list("DA_folds" = DA_folds, "dissimilarity_map" = prob_raster, "threshold" = threshold))
}
