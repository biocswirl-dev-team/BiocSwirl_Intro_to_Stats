# Code placed in this file fill be executed every time the
      # lesson is started. Any variables created here will show up in
      # the user's working directory and thus be accessible to them
      # throughout the lesson.

.get_course_path <- function(){
  tryCatch(swirl:::swirl_courses_dir(),
           error = function(c) {file.path(find.package("swirl"),"Courses")}
  )
}

metadata <- read_csv(file.path(.get_course_path(), "BiocSwirl_Intro_to_Stats", "lab02_hypothesis_testing","pannets_metadata.csv"))

rnaseq <- read_csv(file.path(.get_course_path(), "BiocSwirl_Intro_to_Stats", "lab02_hypothesis_testing","pannets_expr_rnaseq.csv.gz"))

rnaseq_long <- pivot_longer(rnaseq, cols = -Gene, 
                            names_to = "Tumour",
                            values_to = "Expr")

rnaseq_wide <- pivot_wider(rnaseq_long, id_cols = Tumour,
                           names_from = Gene, 
                           values_from = Expr)