pts <- sf::st_as_sf(data.frame(ID = 1:10, x = runif(10), y = runif(10)), 
                coords = c("x", "y"), crs = 4326)
env <- data.frame(ID = 1:10, var1 = rnorm(10), var2 = runif(10), target = runif(10))
pts <- cbind(pts, env)
response_name <- "target"
folds <- spatial_plus_cv(samples = pts, response_name = response_name,
                        cate_col_start = 0, cate_col_end = 0,
                        k = 3, sp_threshold = 1, method = "SP")