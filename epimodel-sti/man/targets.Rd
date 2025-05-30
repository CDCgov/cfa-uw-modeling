% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nw_prep.R
\name{targets}
\alias{targets}
\alias{calc_joint_nodefactor}
\alias{calc_single_attr_nodefactor}
\alias{calc_edges}
\alias{calc_absdiff}
\alias{calc_concurrent}
\alias{calc_cross_network}
\alias{calc_nodecov_age}
\alias{check_targets}
\alias{check_conditions}
\title{Target Stats Calculation Helpers}
\usage{
calc_joint_nodefactor(params, attrs, joint_attrs, joint_name, rel)

calc_single_attr_nodefactor(params, attr_name, joint_attrs, nf_joint_counts)

calc_edges(nf_joint_counts)

calc_absdiff(params, rel, count_type, edges)

calc_concurrent(params, rel, num)

calc_cross_network(params, rel)

calc_nodecov_age(
  params,
  rel,
  attr_name,
  edges,
  level = NULL,
  joint_attrs,
  nf_joint_counts,
  attr_squared
)

check_targets(edges, final_targets, count_type, threshold = 0.01)

check_conditions(nw, params, rel, count_type, attr_name, joint_attrs)
}
\arguments{
\item{params}{the parameter list object with parameters information for
all networks estimated from data}

\item{attrs}{which two attributes?}

\item{joint_attrs}{string of length 2, the two attrs used to calculate joint distribution
of estimates calculated from empirical data}

\item{joint_name}{Name of joint attribute as found in parameter input, calculated in calc_targets()}

\item{rel}{string, which relationship/network sub-lists to pull parameters from.}

\item{attr_name}{string of length 1, used for nodefactor and nodematch targets. If NULL, produces
targets based on the two joint_attrs specifed}

\item{nf_joint_counts}{output from calc_joint_nodefactor()}

\item{count_type}{string, which ERGM term to generate targets for. Currently can only be
in the form "edges", "nodefactor", "nodematch", "absdiff_sqrt_age", or "concurrent".}

\item{edges}{output from calc_edges()}

\item{level}{additional statification for nodedov target calculation function}

\item{attr_squared}{for nodecov target calculation, use squared version of attribute? (usually, age)}

\item{final_targets}{vector, final ergm term targets to be checked before output}

\item{threshold}{default = 0.01, proportion of expected activity based on edges that
calulated target is allowed within (+/- threshold)}

\item{nw}{the network object outputted from "generate_init_network"
(usually already in environment during workflow)}
}
\description{
Small helper functions used as part of calc_targets()
}
