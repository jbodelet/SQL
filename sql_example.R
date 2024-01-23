# devtools::install_github("jbodelet/SQL/sql")
library(sql)


#==============
# q= 1 factor:
#==============
  

sim <- simulate_afm(n = 150, p = 200)

sql <- SQL(sim$data, d = 4)
sql
abs( cor(sim$factor, sql$factor) )
plot(sql)


#==============
# q= 3 factor:
#==============

q <- 3
sim <- simulate_afm(n = 150, p = 200, q = q)

sql <- SQL(sim$data, q = q)
sql
abs( cor(sim$factor, sql$factor) )




