m1 = c("Antoine",
       "Hadrien",
       "Hedi",
       "Leyla",
       "Lamia",
       "Louison",
       "Luis",
       "Maël",
       "Mathias",
       "Rayan")
m2 = c("Adam","LeroyMerlin","Pizza","Camille","Victor","Victor","Dieudo","Daddy","Daddy")
for(i in m2){
  print(i)
  res = sample(m1,1)
  m1 = m1[m1 != res]
  print(res)
}