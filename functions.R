
get_angle <- function(xy) {
  ta <- base::atan2(xy[, 2], xy[, 1])
  ta <- ta - (pi/2)
  ta
}


kappa_to_cos_ta <- function(x) {
  checkmate::assert_number(x, lower = 0)
  c("cos_ta_" = x )
}

scale_to_sl <- function(x) {
  checkmate::assert_number(x, lower = .Machine$double.eps)
  c("sl_" = -1 / x )
}

shape_to_log_sl <- function(x) {
  checkmate::assert_number(x, lower = .Machine$double.eps)
  c("log_sl_" = x - 1)
}

sl_to_scale <- function(x) {
  c("scale" = -1/x)
}

log_sl_to_shape <- function(x) {
  c("shape" = x + 1)
}

cos_ta_to_cappa <- function(x) {
  c("kappa" = x)
}

kernel_setup <- function(template, max.dist = 100, position = c(0, 0)) {
  
  checkmate::assert_class(template, "RasterStack")
  checkmate::assert_number(max.dist, lower = 0)
  checkmate::assert_numeric(position, len = 2)
  
  p <- sf::st_sf(geom = sf::st_sfc(sf::st_point(position))) %>% 
    sf::st_buffer(dist = max.dist)
  
  # 2. Rasterize buffer
  r1 <- raster::rasterize(p, template)
  
  # 3. Get xy from buffer
  xy <- raster::rasterToPoints(r1)
  
  # 4. Substract start, so we are centered around 0,0
  xy[, 1] <- xy[, 1] - position[1]
  xy[, 2] <- xy[, 2] - position[2]
  xy[, 1:2]
  
  k<- tibble(
    x = xy[, 1],
    y = xy[, 2],
    ta_ = get_angle(xy), 
    sl_ = sqrt(x^2 + y^2),
    log_sl_ = log(sl_))
  k
}

kernel_shift <- function(kernel, position) {
  checkmate::assert_tibble(kernel)
  checkmate::assert_numeric(position, len = 2)
  
  kernel$x <- kernel$x + position[1]
  kernel$y <- kernel$y + position[2]
  kernel
}

kernel_rotate <- function(kernel, direction) {
  checkmate::assert_tibble(kernel)
  checkmate::assert_number(direction)
  
  kernel$ta_ <- (kernel$ta_ + direction) %% (2 * pi)
  kernel$cos_ta_ <- cos(kernel$ta_)
  kernel
}

kernel_add_covars <- function(kernel, spatial.covars, position, temporal.covars = NULL) {
  checkmate::assert_tibble(kernel)
  checkmate::assert_class(spatial.covars, "RasterStack")
  checkmate::assert_data_frame(temporal.covars, null.ok = TRUE)
  
  cells <- cellFromXY(spatial.covars, cbind(kernel$x, kernel$y))
  sc1 <- as.data.frame(spatial.covars[cells])
  
  cc <- cellFromXY(spatial.covars, cbind(position[1], position[2]))
  ii <- which(cc == cells)[1]
  sc2 <- sc1[ii, , drop = FALSE]
  names(sc1) <- paste0(names(sc1), "_end")
  names(sc2) <- paste0(names(sc2), "_start")
  
  
  dplyr::bind_cols(kernel[, c("x", "y", "sl_", "ta_", "log_sl_", "cos_ta_")],
                   sc1, sc2, if (!is.null(temporal.covars)) temporal.covars)
  
}

kernel_finish <- function(kernel, formula, coefficients, return.raster = FALSE) {
  
  checkmate::assert_tibble(kernel)
  checkmate::assert_formula(formula)
  checkmate::assert_vector(coefficients, names = "unique")
  checkmate::assert_logical(return.raster)
  
  if (attr(terms(formula), "intercept")) {
    wx <- as.formula(paste("~ 0 + ", as.character(formula)[2]))
  }
  design_matrix <- model.matrix(wx, kernel) 
  
  # Multiply with coefficients
  kernel$dk <- exp((design_matrix %*% coefficients[colnames(design_matrix)])[, 1])
  kernel$dk <- kernel$dk / sum(kernel$dk, na.rm = TRUE)
  
  if(return.raster) {
    raster::rasterFromXYZ(kernel[, c("x", "y", "dk")])
  } else {
    kernel
  }
}

dispersal_kernel <- function(
  wx, coefficients, start, 
  spatial.covars, temporal.covars = NULL, 
  direction = 0, max.dist = 100, return.raster = FALSE) {
  
  dk0 <- kernel_setup(spatial.covars, max.dist = max.dist, position = start)
  dk1 <- kernel_shift(dk0, position = start) 
  dk2 <- kernel_rotate(dk1, direction = direction)
  dk3 <- kernel_add_covars(dk2, spatial.covars = spatial.covars, 
                           position = start,
                           temporal.covars = temporal.covars)
  dk4 <- kernel_finish(dk3, formula = wx, coefficients = coefficients, 
                       return.raster = return.raster)
  
}


  
simulate_track <- function(
  wx, coefficients, start, spatial.covars,
  direction = 0, temporal.covars = NULL, max.dist = 100, n = 10, 
  as.track = TRUE, start.time = ymd_hms("2020-01-01 00:00:00"), 
  delta.time = hours(2)) {
  
  ##
  if (FALSE) {
  
wx = formula
coefficients = coefs
start = c(5000, 5000)
spatial.covars = lscp
               max.dist = 500
               n = 10
  }
###
  
  checkmate::assert_number(n, lower = 1)
  checkmate::assert_logical(as.track)
  checkmate::assert_class(delta.time, "Period")
  
  k <- kernel_setup(spatial.covars, max.dist = max.dist, position = start)
  
  n.px <- nrow(k)
  dir <- direction
  res <- tibble::tibble(x = NA_real_, y = rep(NA_real_, n))
  start1 <- start
  
  xmn <- raster::xmin(spatial.covars)
  xmx <- raster::xmax(spatial.covars)
  ymn <- raster::ymin(spatial.covars)
  ymx <- raster::xmax(spatial.covars)
  
  for (i in 1:n) {
    dk <- kernel_shift(k, position = start1) 
    dk <- kernel_rotate(dk, direction = dir) 
    dk <- kernel_add_covars(dk,
      spatial.covars = spatial.covars, position = start1,
      temporal.covars = if(!is.null(temporal.covars)) 
        temporal.covars[i, , drop = FALSE] else NULL)
    
    dk <- kernel_finish(dk, formula = wx, coefficients = coefficients, 
                        return.raster = FALSE)
    nxt.cell <- wrswoR::sample_int_expj(n.px, size = 1, dk$dk)
    dir <- dir + dk$ta_[nxt.cell]
    start1 <- c(dk$x[nxt.cell], dk$y[nxt.cell])
    res[i, ] <- as.list(start1)
   
    if ((start1[1] - max.dist) < xmn | (start1[1] + max.dist) > xmx |
        (start1[2] - max.dist) < ymn | (start1[2] + max.dist) > ymx) {
      break()
    }
  }
  
  trk <- bind_rows(tibble(x = start[1], y = start[2]), res[1:i, ])
  
  if (as.track) {
    trk$ts <- start.time + 0:i * delta.time
    make_track(trk, x, y, ts)
  } else {
    trk
  }
  
}


library(MASS)

blur <- function(trk,k) {
  n=dim(trk)[1]
  SE<-sqrt(k)
  trk$x_<-round(trk$x_+rnorm(n,0,SE))
  trk$y_<-round(trk$y_+rnorm(n,0,SE))
  return(trk)
}

fit<-function(trk,lscp){
  m<-trk %>% steps %>% random_steps() %>%
    extract_covariates(lscp) %>%
    mutate(log_sl_ = log(sl_), cos_ta_ = cos(ta_))%>% 
   fit_clogit(case_ ~ var + sl_ + log_sl_ + cos(ta_) + strata(step_id_))
  return(c(summary(m)$coef[1,1],summary(m)$coef[1,3]))
}


simex<-function(df){
  fit = lm(variable ~ error+I(error^2), data = df)
  new_df<-data.frame(error=0)
  p<-predict(fit, newdata = new_df, interval = "confidence", type = "response")
  return(p)
}

landscape<-function(autocorr,nug){
  lscp <- NLMR::nlm_gaussianfield(150, 150,nug=nug, resolution = 10,autocorr_range=autocorr,user_seed = 4,rescale = FALSE)
  lscp <- stack(lscp)
  names(lscp) <- "var"
  return(lscp)
}

prepare_track<-function(track){
  trk<-track %>% steps
  trk<-subset(trk, abs(trk$t1_-trk$t2_)<40)
  trk<- trk %>% random_steps %>% amt::extract_covariates(ras)%>% mutate(log_sl_ = log(sl_), cos_ta_ = cos(ta_))
  trk[sapply(trk, is.infinite)] <- NA
  nlcd<-read_csv(here("/Users/clara/Documents/Master Thesis/Code/Crane/nlcd_legend.csv"))
  #rename variable to make join simpler
  colnames(nlcd)[2]<-'nlcd'
  colnames(trk)[13]<-'nlcd'
  #this joins the actual nlcd habitat name and also a consolidated category from 20->8 levels
  trk<-left_join(trk, nlcd[,1:3])
  trk$category<-factor(trk$category)
  trk$category<-relevel(trk$category,"wetlands")
  #replace shrubland by forest because shrubland is never picked
  trk[trk$category=="shrubland",]$category="forest"
  trk$category<-droplevels(trk$category,"shrubland")
  return(trk)
}


find_variance<-function(df){ 
  df$var<-df$sd^2
  mean<-df %>% group_by(error) %>% summarise(mean(var))
  dimmean<-dim(mean)[1]
  second <- c()
  for (j in 1:dimmean){
    store<-filter(df,df$error==mean$error[j])
    s<-0
    for (i in 1:50){
      s<-s+(mean$`mean(var)`[j]-store[i,4])^2
    }
    second<-c(second,s/(49))
    
  }
  mean$var<-mean$`mean(var)`-second
  fit = lm(var ~ error+I(error^2), data = mean)
  new_df<-data.frame(error=0)
  final<-predict(fit, newdata = new_df, interval = "confidence", type = "response")
  return(final[1])
}

ConfidenceInt<-function(beta,sd){
  return(c(beta-1.96*sd,beta+1.96*sd))
}
