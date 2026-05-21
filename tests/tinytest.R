if ( requireNamespace("tinytest", quietly=TRUE) ){
  old_omp_thread_limit <- Sys.getenv("OMP_THREAD_LIMIT", unset = NA_character_)
  Sys.setenv("OMP_THREAD_LIMIT" = "1")
  on.exit({
    if (is.na(old_omp_thread_limit)) {
      Sys.unsetenv("OMP_THREAD_LIMIT")
    } else {
      Sys.setenv("OMP_THREAD_LIMIT" = old_omp_thread_limit)
    }
  }, add = TRUE)

  if ( requireNamespace("data.table", quietly=TRUE) ){
    old_dt_threads <- data.table::getDTthreads()
    data.table::setDTthreads(1L)
    on.exit(data.table::setDTthreads(old_dt_threads), add = TRUE)
  }

  tinytest::test_package("automatedRecLin", ncpu = NULL)
}
