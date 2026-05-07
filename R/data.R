#' @title A_example dataset
#'
#' @description
#' An example dataset containing artificial personal data.
#'
#' @format
#' A `data.frame` with 10 records. Each row represents one record, with the following columns:
#' `name`, `surname`, and `city`. Some records can be matched with records in the
#' `B_example` dataset.
#'
#' @docType data
#' @keywords datasets
#' @name A_example
#' @rdname A_example
#' @examples
#' data("A_example")
#' A_example
"A_example"

#' @title B_example dataset
#'
#' @description
#' An example dataset containing artificial personal data.
#'
#' @format
#' A `data.frame` with 12 records. Each row represents one record, with the following columns:
#' `name`, `surname`, and `city`. Some records can be matched with records in the
#' `A_example` dataset.
#'
#' @docType data
#' @keywords datasets
#' @name B_example
#' @rdname B_example
#' @examples
#' data("B_example")
#' B_example
"B_example"

#' @title Fictional census data
#'
#' @description
#' This dataset was created by Paula McLeod, Dick Heasman and Ian Forbes, ONS,
#' for the ESSnet DI on-the-job training course, Southampton, 25-28 January 2011.
#' It contains fictional data representing some observations from a decennial Census.
#'
#' @format A `data.table` with 25343 records. Each row represents one record, with the following columns:
#' \itemize{
#' \item{`person_id` -- a unique number for each person, consisting of postcode, house number and person number,}
#' \item{`pername1` -- forename,}
#' \item{`pername2` -- surname,}
#' \item{`sex` -- gender (M/F),}
#' \item{`dob_day` -- day of birth,}
#' \item{`dob_mon` -- month of birth,}
#' \item{`dob_year` -- year of birth,}
#' \item{`hse_num` -- house number, a numeric label for each house within a street,}
#' \item{`enumcap` -- an address consisting of house number and street name,}
#' \item{`enumpc` -- postcode,}
#' \item{`str_nam` -- street name of person's household's street,}
#' \item{`cap_add` -- full address, consisting of house number, street name and postcode,}
#' \item{`census_id` -- person ID with "CENS" added in front.}
#' }
#'
#' @references
#' McLeod, P., Heasman, D., Forbes, I. (2011). Simulated data for the ESSnet DI on-the-job training course,
#' Southampton, 25-28 January 2011.
#' \url{https://wayback.archive-it.org/12090/20231221144450/https://cros-legacy.ec.europa.eu/content/job-training_en}
#'
#' @docType data
#' @keywords datasets
#' @name census
#' @rdname census
#' @examples
#'
#' data("census")
#' head(census)
#'
"census"

#' @title Fictional customer data
#'
#' @description
#' This dataset was created by Paula McLeod, Dick Heasman and Ian Forbes, ONS,
#' for the ESSnet DI on-the-job training course, Southampton, 25-28 January 2011.
#' It contains fictional observations from Customer Information System,
#' which is combined administrative data from the tax and benefit systems.
#'
#' @format A `data.table` with 24613 records. Each row represents one record, with the following columns:
#' \itemize{
#' \item{`person_id` -- a unique number for each person, consisting of postcode, house number and person number,}
#' \item{`pername1` -- forename,}
#' \item{`pername2` -- surname,}
#' \item{`sex` -- gender (M/F),}
#' \item{`dob_day` -- day of birth,}
#' \item{`dob_mon` -- month of birth,}
#' \item{`dob_year` -- year of birth,}
#' \item{`enumcap` -- an address consisting of house number and street name,}
#' \item{`enumpc` -- postcode,}
#' \item{`cis_id` -- person ID with "CIS" added in front.}
#' }
#'
#' @references
#' McLeod, P., Heasman, D., Forbes, I. (2011). Simulated data for the ESSnet DI on-the-job training course,
#' Southampton, 25-28 January 2011.
#' \url{https://wayback.archive-it.org/12090/20231221144450/https://cros-legacy.ec.europa.eu/content/job-training_en}
#'
#' @docType data
#' @keywords datasets
#' @name cis
#' @rdname cis
#' @examples
#'
#' data("cis")
#' head(cis)
#'
"cis"
