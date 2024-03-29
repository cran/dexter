% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dexter.R
\name{start_new_project}
\alias{start_new_project}
\title{Start a new project}
\usage{
start_new_project(rules, db_name = "dexter.db", person_properties = NULL)
}
\arguments{
\item{rules}{A data frame with columns \code{item_id}, \code{response}, and \code{item_score}.
The order is not important but spelling is. Any other columns will be ignored.}

\item{db_name}{A string specifying a filename
for a new sqlite database to be created. If this name does not
contain a path, the file will be created in the work
directory. Any existing file with the same name will be overwritten. For an in-memory database
you can use the string \code{":memory:"}. A connection object is also allowed.}

\item{person_properties}{An optional list of person properties. Names should correspond to person_properties intended to be used in the project.
Values are used as default (missing) values. The datatype will also be inferred from the values.
Known person_properties will be automatically imported when adding response data with \code{\link{add_booklet}}.}
}
\value{
a database connection object.
}
\description{
Imports a complete set of scoring rules and starts a new project (database)
}
\details{
This package only works with closed items (e.g. likert, MC or possibly short answer)
it does not score any open items.
The first step to creating a project is to import an exhaustive list of all items and
all admissible responses, along with the score that any of the latter will be given.
Responses may be integers or strings but they will always be treated as strings.
Scores must be integers, and the minimum score for an item must be 0.
When inputting data, all responses not specified in the rules can optionally be treated as
missing and ultimately scored 0, but it is good style to include the missing
responses in the list. NA values will be treated as the string "NA"'.
}
\examples{
\donttest{
head(verbAggrRules)
db_name = tempfile(fileext='.db')
db = start_new_project(verbAggrRules, db_name, 
                       person_properties = list(gender = "unknown"))
}

}
