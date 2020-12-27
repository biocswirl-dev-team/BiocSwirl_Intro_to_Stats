# Code placed in this file fill be executed every time the
      # lesson is started. Any variables created here will show up in
      # the user's working directory and thus be accessible to them
      # throughout the lesson.


.get_course_path <- function(){
  tryCatch(swirl:::swirl_courses_dir(),
           error = function(c) {file.path(find.package("swirl"),"Courses")}
  )
}

array <- read_csv(file.path(.get_course_path(), "BiocSwirl_Intro_to_Stats", "lab03_simple_linear_regression","pannets_expr_array.csv.gz"))

array_long <- pivot_longer(array, cols = -Gene, 
                           names_to = "Tumour",
                           values_to = "Expr")
array_wide <- pivot_wider(array_long, id_cols = Tumour,
                          names_from = Gene, 
                          values_from = Expr)

actb <- data.frame(
  tumour = rnaseq_wide$Tumour,
  rnaseq = rnaseq_wide$ACTB,
  array = array_wide$ACTB
)