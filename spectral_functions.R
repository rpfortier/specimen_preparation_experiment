### VIP.R: Implementation of VIP (variable importance in projection)(*) for the
### `pls' package.
### $Id: VIP.R,v 1.2 2007/07/30 09:17:36 bhm Exp $

### Copyright © 2006,2007 Bjørn-Helge Mevik
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License version 2 as
### published by the Free Software Foundation.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.

### A copy of the GPL text is available here:
### http://www.gnu.org/licenses/gpl-2.0.txt

### Contact info:
### Bjørn-Helge Mevik
### bhx6@mevik.net
### Rødtvetvien 20
### N-0955 Oslo
### Norway

### (*) As described in Chong, Il-Gyo & Jun, Chi-Hyuck, 2005, Performance of
### some variable selection methods when multicollinearity is present,
### Chemometrics and Intelligent Laboratory Systems 78, 103--112.

## VIP returns all VIP values for all variables and all number of components,
## as a ncomp x nvars matrix.
VIP <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}


## VIPjh returns the VIP of variable j with h components
VIPjh <- function(object, j, h) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  b <- c(object$Yloadings)[1:h]
  T <- object$scores[,1:h, drop = FALSE]
  SS <- b^2 * colSums(T^2)
  W <- object$loading.weights[,1:h, drop = FALSE]
  Wnorm2 <- colSums(W^2)
  sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS))
}

###########################################################################################
#### The following function is custom and computes the first and second derivatives for hyperspectral data
compute_hyperspectral_derivatives <- function(data, 
                                              wavelength_range = c(450, 2400),
                                              metadata_cols = 1:28) {
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("Input data must be a data frame")
  }
  
  # Identify columns that start with "Wavelength_"
  wavelength_columns <- grep("^Wavelength_", names(data), value = TRUE)
  
  if (length(wavelength_columns) == 0) {
    stop("No columns found starting with 'Wavelength_'. Please check your data format.")
  }
  
  # Convert the selected wavelength columns to numeric
  data[wavelength_columns] <- lapply(data[wavelength_columns], as.numeric)
  
  # Extract wavelengths from column names
  wavelengths <- as.numeric(gsub("Wavelength_", "", wavelength_columns))
  
  # Filter to specified wavelength range
  wavelength_filter <- wavelengths >= wavelength_range[1] & wavelengths <= wavelength_range[2]
  spectral_cols <- wavelength_columns[wavelength_filter]
  wavelengths_restricted <- wavelengths[wavelength_filter]
  
  if (length(spectral_cols) == 0) {
    stop(paste0("No wavelength columns found in range ", 
                wavelength_range[1], "-", wavelength_range[2], " nm"))
  }
  
  message(paste0("Processing ", length(spectral_cols), 
                 " wavelength columns (", min(wavelengths_restricted), 
                 "-", max(wavelengths_restricted), " nm)"))
  
  # Convert spectral columns to matrix
  spectra_matrix <- as.matrix(data[, spectral_cols])
  
  # Compute 1st derivative
  # Take 1st derivative across wavelengths (columns)
  # Using t(diff(t(matrix))) to apply diff across columns
  spectra_deriv <- t(diff(t(spectra_matrix), differences = 1))
  message(paste0("First derivative: ", ncol(spectra_deriv), " columns"))
  
  # Compute 2nd derivative
  spectra_deriv2 <- t(diff(t(spectra_matrix), differences = 2))
  message(paste0("Second derivative: ", ncol(spectra_deriv2), " columns"))
  
  # Adjust column names for first derivative
  # Remove first wavelength since diff reduces by 1
  new_wavelengths <- wavelengths_restricted[-1]
  colnames(spectra_deriv) <- paste0("Deriv_Wavelength_", new_wavelengths)
  
  # Adjust column names for second derivative
  # Remove first two wavelengths since diff with differences=2 reduces by 2
  new_wavelengths2 <- wavelengths_restricted[-c(1, 2)]
  colnames(spectra_deriv2) <- paste0("Deriv2_Wavelength_", new_wavelengths2)
  
  # Concatenate first and second derivatives
  spectra_deriv_total <- cbind(spectra_deriv, spectra_deriv2)
  
  # Convert to data frame
  spectra_deriv_df <- as.data.frame(spectra_deriv_total)
  
  # Add metadata back
  # Handle both numeric and character metadata_cols
  if (is.numeric(metadata_cols)) {
    metadata <- data[, metadata_cols, drop = FALSE]
  } else {
    metadata <- data[, metadata_cols, drop = FALSE]
  }
  
  hyperspectral_deriv_data <- cbind(metadata, spectra_deriv_df)
  
  message(paste0("Output dimensions: ", nrow(hyperspectral_deriv_data), 
                 " rows × ", ncol(hyperspectral_deriv_data), " columns"))
  
  return(hyperspectral_deriv_data)
}

########################################################################################
### The remaining functions are taken from Kothari et al. (2023)
### https://github.com/ShanKothari/pressed-leaf-models/blob/main/00%20useful_functions.R
###

## root mean squared deviation
RMSD<-function(measured,predicted){
  not.na<-which(!is.na(measured) & !is.na(predicted))
  return(sqrt(sum((measured-predicted)^2,na.rm=T)/(length(not.na)-1)))
}

## percent RMSD (based on data quantiles)
## set min and max to 0 and 1 for range as denominator
## or to 0.25 and 0.75 for IQR as denominator
percentRMSD<-function(measured,predicted,min,max,na.rm=T){
  RMSD_data<-RMSD(measured,predicted)
  range<-unname(quantile(measured,probs=max,na.rm=na.rm)-quantile(measured,probs=min,na.rm=na.rm))
  return(RMSD_data/range)
}

## applying coefficients to validation spectra
apply.coefs<-function(coef.list,val.spec,intercept=T){
  if(sum(lapply(coef.list,length)==ncol(val.spec)+intercept) < length(coef.list)){
    stop("some coefficients have the wrong length")
  }
  
  coef.matrix<-matrix(unlist(coef.list),
                      nrow=length(coef.list),
                      byrow=T)
  
  if(intercept==T){
    pred.matrix<-t(t(as.matrix(val.spec) %*% t(coef.matrix[,-1]))+coef.matrix[,1])
  } else {
    pred.matrix<-as.matrix(val.spec) %*% t(coef.matrix)
  }
}


