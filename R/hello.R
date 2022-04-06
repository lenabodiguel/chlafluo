#' @Calculation of parameters from RLC process_psi (https://github.com/bmjesus/psifluo/blob/master/R/process_psi.R)
#' @rETR, NPQ, YNPQ, Eff, Qp, SigmaPSII, Fv', Fv/Fm, Tau, Ctau1
#' @title function to fit several PI models
#' @description  Preview of RLC curves and associated parameters
#' @param x rlc object from the process_psi() function
#' @param starting_values list of starting values to fit the rETR model, needs to be changed to auto values
#' @export
#' @keywords external

rlc_parameters <- function(x, strain, growth_light, protocol, duration, used){
  #x = psifluo::process_psi(file_name =  , calib_file = , plot_matrix = TRUE, str_model = "tau1_subset", subset_time = 300)

  #light
  par<-x$data$measuring_steps$par
  #number of light steps
  n_light<-length(x$data$measuring_steps$par)
  #measurement type (light or after 2 seconds of dark, i.e., light, dark)
  type<-c(rep("light",n_light),rep("dark",n_light))

  #PSII quantum efficiency values
  eff<-x$sti_parameters$psII_eff_sti
  #relative ETR
  retr <- par * eff

  #NPQ
  #(Fm-Fm')/Fm'
  # fm_light<-x$sti_parameters$fm_sti[1]
  # fm_dark<-x$sti_parameters$fm_sti[n_light + 1]
  fm_light<-max(x$sti_parameters$fm_sti[1:11])
  fm_dark<-max(x$sti_parameters$fm_sti[(n_light + 1):(n_light * 2)])
  fmp_light<-x$sti_parameters$fm_sti[1:n_light]
  fmp_dark<-x$sti_parameters$fm_sti[(1 + n_light):(n_light * 2)]
  npq_light<-(fm_light-fmp_light)/fmp_light
  npq_dark<-(fm_dark-fmp_dark)/fmp_dark

  npq<-c( npq_light, npq_dark )

  #same as before but using the maximum Fm observed in the series instead of
  #the same value
  #NPQm
  fmm_light<-max(na.omit(fmp_light))
  fmm_dark<-max(na.omit(fmp_dark))
  npq_light_m<-(fmm_light-fmp_light)/fmp_light
  npq_dark_m<-(fmm_dark-fmp_dark)/fmp_dark

  npq_m<-c( npq_light_m, npq_dark_m )


  #YNPQ
  #(F/Fm')-(F/Fm)
  f_light<-x$sti_parameters$fo_sti[1:n_light]
  f_dark<-x$sti_parameters$fo_sti[(1 + n_light):(n_light * 2)]
  ynpq_light<-(f_light/fmp_light)-(f_light/fm_light)
  ynpq_dark<-(f_dark/fmp_dark)-(f_dark/fm_dark)

  ynpq<-c( ynpq_light, ynpq_dark )


  #same as before but using the maximum Fm observed in the series instead of
  #the same value
  #YNPQm

  ynpqm_light<-(f_light/fmp_light)-(f_light/fmm_light)
  ynpqm_dark<-(f_dark/fmp_dark)-(f_dark/fmm_dark)

  ynpqm<-c( ynpqm_light, ynpqm_dark )


  #Fv/Fm
  #(Fm-F0)/Fm
  # fmp_light<-x$sti_parameters$fm_sti[1:n_light]
  # fmp_dark<-x$sti_parameters$fm_sti[(1 + n_light):(n_light * 2)]

  Fv_Fm_light<-(fmp_light-f_light)/fmp_light
  Fv_Fm_dark<-(fmp_dark-f_dark)/fmp_dark
  Fv_Fm <- c(Fv_Fm_light, Fv_Fm_dark)

  # Fv_Fm_light<-eff[1:n_light]
  # Fv_Fm_dark<-eff[(1 + n_light):(n_light * 2)]
  # Fv_Fm <- c(Fv_Fm_eff_light, Fv_Fm_eff_dark)


  #Qp
  #(Fm'-F)/(Fm'-F'0)
  qp<-(fmp_light-f_light)/(fmp_light-x$sti_parameters$fo_sti[12])

  #sigma
  sigma_p_light<-x$sti_parameters$sigma_sti[1:n_light]
  sigma_light<-sigma_p_light[1]

  sigma_p_dark<-x$sti_parameters$sigma_sti[(1 + n_light):(n_light * 2)]
  sigma_dark<-sigma_p_dark[1]

  sigma<-c(sigma_p_light,sigma_p_dark)

  sigma_se_light<-x$sti_parameters$sigma_se_sti[1:n_light]
  sigma_se_dark<-x$sti_parameters$sigma_se_sti[(1 + n_light):(n_light * 2)]
  sigma_se<-c(sigma_se_light, sigma_se_dark)


  #Fvp
  #Fm-F0
  fvp_light <- x$sti_parameters$fm_sti[1:n_light] - x$sti_parameters$fo_sti[1:n_light]
  fvp_dark <- x$sti_parameters$fm_sti[(1 + n_light):(n_light * 2)] - x$sti_parameters$fo_sti[(1 + n_light):(n_light * 2)]

  fvp<-c(fvp_light,fvp_dark)

  # fmp_light<-x$sti_parameters$fm_sti[1:n_light]
  # fmp_dark<-x$sti_parameters$fm_sti[(1 + n_light):(n_light * 2)]
  # f_light<-x$sti_parameters$fo_sti[1:n_light]
  # f_dark<-x$sti_parameters$fo_sti[(1 + n_light):(n_light * 2)]

  #tau
  tau_light<-x$str_parameters$tau1_str[1:n_light]
  tau_dark<-x$str_parameters$tau1_str[(1 + n_light):(n_light * 2)]

  tau<-c(tau_light,tau_dark)

  tau_se_light<-x$str_parameters$tau1_se_str[1:n_light]
  tau_se_dark<-x$str_parameters$tau1_se_str[(1 + n_light):(n_light * 2)]
  tau_se<-c(tau_se_light, tau_se_dark)

  #rho
  rho_light<-x$sti_parameters$rho_sti[1:n_light]
  rho_dark<-x$sti_parameters$rho_sti[(1 + n_light):(n_light * 2)]
  rho<-c(rho_light,rho_dark)

  rho_se_light<-x$sti_parameters$rho_se_sti[1:n_light]
  rho_se_dark<-x$sti_parameters$rho_se_sti[(1 + n_light):(n_light * 2)]
  rho_se<-c(rho_se_light, rho_se_dark)


  #Ctau1
  #Estimation of the PSII closure
  #Ctau1 = 1/(1+(sigmaPSII*PAR*Tau1))
  #SigmaPSII in m2.photons
  #PAR in photons.m-2.s-1
  #Tau1 in s
  Ctau1 <- 1/(1+(sigma_p_light*par*tau_light))

  #Strain
  Strain <- c(strain, strain,strain,strain,strain,strain,strain,strain,strain,strain,strain,strain,strain,strain,strain,strain,strain,strain,strain,strain,strain,strain)

  #growth light in uE
  GL_uE <- c(growth_light, growth_light, growth_light, growth_light, growth_light, growth_light, growth_light, growth_light, growth_light, growth_light, growth_light, growth_light, growth_light, growth_light, growth_light, growth_light, growth_light, growth_light, growth_light, growth_light, growth_light, growth_light)

  #protocol
  Protocol<-c(protocol, protocol, protocol, protocol, protocol, protocol, protocol, protocol, protocol, protocol, protocol, protocol, protocol, protocol, protocol, protocol, protocol, protocol, protocol, protocol, protocol, protocol)
  #duration in second
  Duration_s <- c(duration, duration, duration, duration, duration, duration, duration, duration, duration, duration, duration, duration, duration, duration, duration, duration, duration, duration, duration, duration, duration, duration)
  #sample use
  if (used==0) {
    use <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  }
  else if (used==1) {
    use <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
  }
  else {
    use <- c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
  }

  output <- data.frame(Strain, GL_uE, Protocol, Duration_s, use, type, par, retr, npq, npq_m, ynpq, ynpqm, eff, qp, sigma, fvp, Fv_Fm, tau, Ctau1, sigma_se, tau_se, rho_se,  stringsAsFactors=FALSE)
  return(output)
}
