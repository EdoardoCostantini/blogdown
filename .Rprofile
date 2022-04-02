# REMEMBER to restart R after you modify and save this file!

# First, execute the global .Rprofile if it exists. You may configure blogdown
# options there, too, so they apply to any blogdown projects. Feel free to
# ignore this part if it sounds too complicated to you.
if (file.exists("~/.Rprofile")) {
  base::sys.source("~/.Rprofile", envir = environment())
}

# Now set options to customize the behavior of blogdown for this project. Below
# are a few sample options; for more options, see
# https://bookdown.org/yihui/blogdown/global-options.html
options(
  # to automatically serve the site on startup, set this option to TRUE
  blogdown.serve_site.startup = FALSE,

  # Knit on save is nice!
  blogdown.knit.on_save = TRUE,

  # build .Rmd to .html (via Pandoc); to build to Markdown, set this option to 'markdown'
  blogdown.method = 'html',

  # Default new blog values
  blogdown.author = "Edoardo Costantini",
  blogdown.ext = ".Rmd",
  blogdown.subdir = "post"
)

# fix Hugo version
options(blogdown.hugo.version = "0.94.2")
