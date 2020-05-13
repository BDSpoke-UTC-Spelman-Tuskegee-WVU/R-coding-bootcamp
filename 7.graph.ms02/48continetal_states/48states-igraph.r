#continental US demo, igraph usage

require(RCurl)
require(igraph)
URL = "https://docs.google.com/spreadsheet/pub?key=0ArLJZixvTlU7dDI5c3dnTzdRX1dndzBORUk1UC0wYlE&single=true&gid=0&output=csv"
states = read.csv(textConnection(getURL(URL)), colClass = c("character", "character"))

g = graph.data.frame(states, directed=F)
g.degree = degree(g)
g.degree [g.degree == max(g.degree)] #TN and MO have 8 bordering states

g.shortestpath.m = shortest.paths(g)
str(g.shortestpath.m)
sorted.names = sort( rownames(g.shortestpath.m) )
gsm = g.shortestpath.m[, sorted.names]
gsm = gsm[sorted.names, ]

URL2 = "https://docs.google.com/spreadsheet/pub?key=0ArLJZixvTlU7dC1rTlRhTkw5UlY0eWZGT1JOVUVCWGc&single=true&gid=0&output=csv"
state.year.tb = read.csv(textConnection(getURL(URL2)), colClass = c("character", NA))


