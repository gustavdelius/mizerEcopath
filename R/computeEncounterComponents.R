#' Set up intake rate and metabolic loss parameters for a tuned non-interacting
#' model. Take the existing tuned encounter rates and compute the correct
#' internal component values, and sets them based on a supplied assimilation
#' efficiency, critical feeding level,  and feeding level. Defaults for these
#' are set following Andersen (2019).
#'
#'
#' @param params A MizerParams object
#' @param fc Critical feeding level. Default is 0.2
#' @param f Feeding level. Default is 0.6
#' @return A MizerParams object with the adjusted max intake rate (h) and
#' metablolic loss (k) parameters
#' @export
computeEncounterComponents<-function(params, fc = 0.2, f = 0.6) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params

    if (any(sp$h != Inf)) {
        stop("initial h values must be Inf")
    }
    if (any(sp$ks != 0)) {
        stop("initial ks values must be 0")
    }

    Eiw<-getExtEncounter(params)

    for(i in seq(nrow(sp))){
        a<-sp[i,]$alpha
        Eriw<-Eiw[i,]*a

        kiw<-a*(1-fc)*Eiw[i,]

        Eiw[i,]<-(Eriw+kiw)/(a*(1-f))

        hiw<-(Eiw[i,]/f)-Eiw[i,]

        ki<-kiw[1]/(w(params)[1]^sp[i,]$n)
        hi<-hiw[1]/(w(params)[1]^sp[i,]$n)
        sp[i,]$ks<-ki
        sp[i,]$h<-hi
    }
    params@species_params<-sp
    ext_encounter(params)<-Eiw

    return(params)
}
