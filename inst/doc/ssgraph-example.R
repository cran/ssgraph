## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set( collapse = TRUE, comment = ">", fig.width = 7, fig.height = 7, fig.align = "center" )

## ----eval = FALSE-------------------------------------------------------------
#  install.packages( "ssgraph" )

## ----loadpkg, message = FALSE, warning = FALSE--------------------------------
library( ssgraph )

## ----fig.align = 'center'-----------------------------------------------------
set.seed( 10 )

data.sim <- bdgraph.sim( n = 100, p = 8, graph = "scale-free", vis = TRUE )

round( head( data.sim $ data, 4 ), 2 )

## ----fig.align = 'center'-----------------------------------------------------
ssgraph.obj <- ssgraph( data = data.sim, method = "ggm", iter = 5000, save = TRUE )

summary( ssgraph.obj )

## ----fig.align = 'center'-----------------------------------------------------
compare( data.sim, ssgraph.obj, main = c( "Target", "ssgraph" ), vis = TRUE )

