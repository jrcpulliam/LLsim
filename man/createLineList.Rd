% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/llFunction.R
\name{createLineList}
\alias{createLineList}
\title{Convert simulation output to line list}
\usage{
createLineList(
  cases,
  obs_prob,
  mean_report_delay,
  mean_confirm_delay,
  zero_date,
  obs_date = as.character(Sys.Date())
)
}
\arguments{
\item{cases}{A data frame generated by \code{simpleSim} with rows
representing all cases that occur in the outbreak simulation, having at
minimum the following columns: \code{id} (unique individual ID),}

\item{obs_prob}{The individual-level probabilty that a case is reported
(numeric)}

\item{mean_report_delay}{The average number of days between disease onset and
reporting to a health facility (numeric)}

\item{mean_confirm_delay}{The average number of days between reporting to a
health facility and confirmation of a case (numeric)}

\item{zero_date}{Character string representing the date on which the initial
case was infected (YYYY-MM-DD)}

\item{obs_date}{Character string representing the date on which observation
of the state of the line list takes place (YYYY-MM-DD; default = system
date)}
}
\value{
Returns a data frame the provides the line list, as observed on
  \code{obs_date}

  TODO: Add more detail
}
\description{
\code{convertLineList} takes simulation output on cases, converts numeric
time to dates, and adds an observation layer to output observed cases on a
specified date.
}
\details{
TODO: Add more details later
}
