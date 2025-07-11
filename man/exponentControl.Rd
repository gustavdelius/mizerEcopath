% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exponentControl.R
\name{exponentControl}
\alias{exponentControl}
\alias{exponentControlUI}
\title{Controlling the allometric exponents}
\usage{
exponentControl(input, output, session, params, params_old, flags, ...)

exponentControlUI(params, input)
}
\arguments{
\item{input}{Reactive holding the inputs}

\item{output}{Reactive holding the outputs}

\item{session}{Shiny session}

\item{params}{The MizerParams object currently being tuned.}

\item{params_old}{Reactive value holding non-updated MizerParams object}

\item{flags}{Environment holding flags to skip certain observers}

\item{...}{Unused}
}
\value{
A tagList with sliders for the exponents
}
\description{
Controlling the allometric exponents
}
