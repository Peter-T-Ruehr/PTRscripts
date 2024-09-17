#' Add leading zero to a vector of numbers.
#' @param numbers A vector of `numeric` or `character` values.
#' @param length A single `numeric` value defining the desired length of each
#' `character` of the resulting vector.
#' @return A vector of `characters`.
#' @export
#' @examples
#' # add leading zeros to numbers
#' add_leading_zeros(numbers = c(1, 2, 3),
#'                    length = 4)
#'
#' # add leading zeros to characters
#' add_leading_zeros(numbers = c("01", "2", "003"),
#'                    length = 4)
add_leading_zeros <- function(numbers, length){
  require(stringr)
  str_pad(numbers, length, pad = "0")
}


#' Create empty tibble with column names and column types.
#' @param nrow `Numeric` value defining number of desired rows. `Default = 0`.
#' @param names A vector of `characters` for desired column names.
#' @param type A `character` value for the desired type of all columns.
#' Takes `number`, 'character`, or `factor`. `Default = number`.
#' @export
#' @examples
#' empty_tibble(nrow = 0,
#'               names = c("col1", "col2", "col3"),
#'               type = "character")
empty_tibble <- function(nrow = 0,
                         names,
                         type = "number"){
  tibble <- as_tibble(setNames(data.frame(matrix(nrow = nrow, ncol = length(names))), names))
  if(type == "character") {tibble %>% mutate_all(as.character)}
  else if(type == "factor") {tibble %>% mutate_all(as.factor)}
  else {tibble %>% mutate_all(as.numeric)}
}


#' Show loop progress in percent.
#' @param current `Numeric` value of current loop iteration.
#' @param end `Numeric` value of loop end value.
#' @export
#' @examples
#' n = 100
#' for(i in 1:n){
#'   print_progress(i, n)
#' }
print_progress <- function(current,
                           end){
  cat("\r", paste0(round(current/end*100,2), "%..."))
}
