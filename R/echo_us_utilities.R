# Trickery to handle working directories
if(isTRUE(file.info("data")$isdir)){
  # working directory contains data subdirectory
  us_properties <- read.csv("data/acoustic_properties_thigh.csv")
} else if(isTRUE(file.info("../data")$isdir)){
  # parent of working directory contains data
  us_properties <- read.csv("../data/acoustic_properties_thigh.csv")
} else {
  stop(paste("Needed file, data/acoustic_properties_thigh.csv, wasn't found.",
             "Your working directory is ", getwd(), "."))
}

# Given a vector of tissue IDs (integers between 3 and 11) and a pixel size, dx,
# return the times and amplitudes of reflected responses to a unit impulse.
#
#' @param tissue_ids a vector of integers in [3, 11]
#' @param dx pixel size in mm.
#' @param us_properties a data frame of ultrasound properties
#' @return a data frame with columns: ID, time, amplitude
impulseResponse <- function(tissue_ids, dx){
  n <- length(tissue_ids)
  # Persiflage to retain our conventional tissue/ID association
  temp <- us_properties$ID
  idrow <- sapply(tissue_ids, function(id)which(id == temp))
  # Round-trip times (in milliseconds) from each of the n-2 interior boundaries
  speeds <- us_properties[idrow, "c"]
  # The factor of 10^6 comes from 10^3 ms/s and 10^3 mm/mm. Speeds are in m/s.
  tof <- 2*1e+6*dx*(1:(n-1)/speeds[1:(n-1)])
  # Reflection coefficients from each of the boundaries
  Z1 <- us_properties[idrow[-n], "Z"]
  Z2 <- us_properties[idrow[-1], "Z"]
  rcs <- ((Z1-Z2)/(Z1+Z2))^2
  # Recursively calculate reflected and transmitted amplitudes at each boundary.
  # NOTE: This is the way I THINK it works. Confident, but not certain. If, for
  # example, forward transmission were a beam--which it can be but doesn't seem
  # to be in our case--then forward losses would be much less than inverse square.
  famp <- 1.0 # initial forward amplitude
  refl <- numeric(n-1) # array to hold reflected amplitudes
  for(i in 1:(n-1)){
    # inverse square loss before reflection
    famp <- famp*(i/(i+1))^2
    # reflected
    refl[i] <- rcs[i]*famp
    # transmitted
    famp <- famp - refl[i]
  }
  # Recursively calculate the amplitude of each reflection on its return trip.
  amp <- numeric(n-1)
  for(i in 1:(n-1)){
    bamp <- refl[i] # reflected amplitude at the ith internal boundary
    for(j in i:1){
      # inverse square loss before reflection
      bamp <- bamp*(j/(j+1))^2
      # backwards transmission losses at an internal boundary (if present)
      # NOTE: neglecting any other effect from reflection of the backward impulse.
      if(j > 1){
        bamp <- bamp*(1-refl[j-1])
      }
    }
    amp[i] <- bamp
  }
  # Return results in a data frame
  data.frame(ID=tissue_ids[-n], tof=tof, amplitude=amp)
}