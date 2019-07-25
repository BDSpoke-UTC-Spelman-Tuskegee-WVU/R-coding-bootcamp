# write a R function to calculate the amount of NaCl 
# in grams for x minimolar (mM) y ml solution
# The formula weight of NaCl is 58.443. Hence, for 1.0 liter 1M solution, 
# we need 58.443 gram of NaCl.
# Note: 1mM = 0.001M, 1ml=0.001liter

molar_solution <- function( x, y) {
  58.433 * (x/10^3) * (y/10^3)
}
molar_solution(100, 500)


molar_solution2 <- function( conc, vol, FW) {
  FW * conc * vol /10^6
}

molar_solution3 <- function( conc, vol, formula) {
  FWs = c(58.443, 74.5513, 84.997, 40)
  names(FWs) = c('NaCl', 'KCl', 'NaNO3', 'NaOH');
  FW = FWs[formula]
  FW * conc * vol / 10^6
}

molar_solution( 1, 100)
molar_solution2( 1, 100, 58.433)
molar_solution3( 100, 1000, 'NaNO3' )
molar_solution2( 100, 1000, 74.5513) #for Kcl
molar_solution3( 1000, 1000, 'NaCl' )

