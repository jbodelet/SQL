# devtools::install_github("jbodelet/SQL/sql")
library(sql)


#==============
# q= 1 factor:
#==============
  
set.seed(123456)
sim <- sql:::simulate_afm(n = 150, p = 200, q = 1, sde = 1)

sql <- SQL(sim$data, d = 4)
sql
abs( cor(sim$factor_z, sql$factor) )
plot(sql)


