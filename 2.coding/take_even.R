take_even = function( x ) {
  y = c(); # a locale copy inside of loop
  for( i in 1:length(x)) {
    if ( (x[i]%% 2 )== 0 ) { # x mod 2
      y = c(y, x[i] );  # add a new x[i] to y
    }
  }
  y; #the last line, return y to the main program
}
