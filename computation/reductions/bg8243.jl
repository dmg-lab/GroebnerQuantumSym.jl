#Reduction of Bg8243

n = 6
G1 = g1_named(n)
u = magic_unitary(n)
bege8243 = bg(8,2,4,3; u=u)
#reduction_string(G1,bege8243; debug=true)

wesc_iterators = [(i,j,k,t) for i = 2:n for j = 2:n for k = 2:n for t in 1:n if j != k];
wesc = [G1["wel$(i)$(j)$(k)"]*u[1,t] for (i,j,k,t) in wesc_iterators];
wesc_names = indexed_name("wesc", [parse(Int, "$i$j$k$t") for (i,j,k,t) in wesc_iterators]);
wesc_ids = [Symbol("wesc$i$j$k$t") for (i,j,k,t) in wesc_iterators];

iesc_iterators = [(i,j,k,t) for i = 2:n for j = 2:n for k = 2:n for t in 1:n if i != k];
iesc = [G1["inj$(i)$(j)$(k)"]*u[t,1] for (i,j,k,t) in iesc_iterators];
iesc_names = indexed_name("iesc", [parse(Int, "$i$j$k$t") for (i,j,k,t) in iesc_iterators]);
iesc_ids = [Symbol("iesc$i$j$k$t") for (i,j,k,t) in iesc_iterators];

uucs_iterators = [(i,j,k,h,f) for i = 1:n for j = 1:n for k = 2:n for h = 2:n for f = 2:n if i != k && j != h && f !=h];
uucs = [u[i,j]*u[k,h]*G1["cs$f"]-u[i,j]*G1["wel$(k)$(h)$(f)"] for (i,j,k,h,f) in uucs_iterators];
uucs_names = indexed_name("uucs", [parse(Int, "$i$j$k$h$f") for (i,j,k,h,f) in uucs_iterators]);
uucs_ids = [Symbol("uucs$i$j$k$h$f") for (i,j,k,h,f) in uucs_iterators];

uurs_iterators = [(i,j,k,h,f) for i = 2:n for j = 2:n for k = 2:n for h = 2:n for f = 2:n if i != k && j != h && f !=k];
uurs = [u[i,j]*u[k,h]*G1["rs$f"]-u[i,j]*G1["inj$(k)$(h)$(f)"] for (i,j,k,h,f) in uurs_iterators];
uurs_names = indexed_name("uurs", [parse(Int, "$i$j$k$h$f") for (i,j,k,h,f) in uurs_iterators]);
uurs_ids = [Symbol("uurs$i$j$k$h$f") for (i,j,k,h,f) in uurs_iterators];

ursu_iterators = [(i,j,k,h,f) for i = 2:n for j = 2:n for k = 2:n for h = 2:n for f = 2:n if i != k && k != h];
ursu = [u[i,j]*G1["rs$k"]*u[h,f] - u[i,j]*G1["inj$(k)$(f)$(h)"] - G1["inj$(i)$(j)$(k)"]*u[h,f] for (i,j,k,h,f) in ursu_iterators];
ursu_names = indexed_name("ursu", [parse(Int, "$i$j$k$h$f") for (i,j,k,h,f) in ursu_iterators]);
ursu_ids = [Symbol("ursu$i$j$k$h$f") for (i,j,k,h,f) in ursu_iterators];

ucsu_iterators = [(i,j,k,h,f) for i = 2:n for j = 2:n for k = 2:n for h = 2:n for f = 2:n if k !=j && f != k];
ucsu = [u[i,j]*G1["cs$k"]*u[h,f] - u[i,j]*G1["wel$(h)$(k)$(f)"] - G1["wel$(i)$(j)$(k)"]*u[h,f] for (i,j,k,h,f) in ucsu_iterators];
ucsu_names = indexed_name("ucsu", [parse(Int, "$i$j$k$h$f") for (i,j,k,h,f) in ucsu_iterators]);
ucsu_ids = [Symbol("ucsu$i$j$k$h$f") for (i,j,k,h,f) in ucsu_iterators];

bg2sr_iterators = [(i,j) for i = 2:n for j = 2:n if i != j && (i,j) != (2,3)];
bg2sr = [sum([G1["bg2_$(i)$(j)$(k)"] for k in 2:n if k != j && (j,k) != (2,3)]) for (i,j) in bg2sr_iterators];
bg2sr_names = indexed_name("bg2sr", [parse(Int, "$i$j") for (i,j) in bg2sr_iterators]);
bg2sr_ids = [Symbol("bg2sr$i$j") for (i,j) in bg2sr_iterators];

function help1(i,j)
  ans = -sum([E1["ucsu$(i)$(j)$(k)$(3)$(3)"] for k in 4:n])
  ans += sum([sum([E1["uucs$(i)$(j)$(k)$(h)$(3)"] for h in 4:n]) for k in 3:n if k != i])
  ans += sum([u[i,j]*E1["rwel$(k)$(3)"] for k in 4:n])
  ans += sum([E1["wesc$(i)$(j)$(k)$(3)"] for k in 4:n])
  ans += sum([E1["bg2sr$(i)$(k)"] for k in 2:n if k != i])
  ans -= sum([E1["bg2_$(i)$(k)3"] for k in 4:n if k != i])
  return ans
end



E1 = named_generators(wesc, wesc_names, wesc_ids)
add!(E1, iesc, iesc_names, iesc_ids)
add!(E1, uucs, uucs_names, uucs_ids)
add!(E1, uurs, uurs_names, uurs_ids)
add!(E1, ursu, ursu_names, ursu_ids)
add!(E1, ucsu, ucsu_names, ucsu_ids)
add!(E1, bg2sr, bg2sr_names, bg2sr_ids)
add!(E1, G1; check=false)

reduction_string(E1, bege8243
-sum([E1["bg2_2$(k)3"] for k in 4:n])+sum([E1["bg8_2$(k)3"] for k in 5:n])
-u[3,2]*E1["cs4"]*u[3,3] +u[3,2]*E1["rwel43"]+E1["wel324"]*u[3,3]
+help1(3,2)-u[3,2]*E1["rwel43"]+E1["ucsu32433"]-E1["wel325"]*u[3,3] - E1["wel326"]*u[3,3] #Scheinbar nicht fine enough
+help1(4,2)+help1(5,2)+help1(6,2)
+sum([sum([E1["ursu2$g$(i)33"] for i in 4:n])-sum([E1["ursu$(h)$g$k$(i)3"] for h in 3:n for i in 2:n for k in 2:n if i != k && k != h && i != 3]) for g in 3:n])

-sum([E1["wesc$h$g$k$i"] for g in 3:n for h in 3:n for k in 2:n for i in 2:n if k != i && k != g && i != 3])
+sum([u[3,g]*E1["wel3$k$i"] for g in 3:n for k in 2:n for i in 2:n if i != 3 && i != k && k != g])

-sum([E1["wel$(g)23"]+4*sum([E1["wel$(g)2$i"] for i in 4:n]) for g in 4:n])
-2*sum([E1["wel$(g)3$i"] for g in 4:n for i in 2:n if i != 3])
+sum([E1["wel$(g)$(h)3"]-2*sum([E1["wel$(g)$h$i"] for i in 2:n if i != 3 && i != h]) for g in 4:n for h in 4:n])


-4*sum([E1["cs2"]*u[3,i] for i in 4:n])

+sum([(4*sum([u[k,2]*E1["cs$i"] for i in 4:n]) -3*u[k,2]*E1["cs3"]) for k in 3:n])

+sum([(6*sum([u[g,h]*E1["cs$i"] for g in 3:n for i in 2:n if i != 3])-3*sum([E1["cs$h"]*u[3,i] for i in 2:n if i != h && i != 3])) for h in 3:n])

+sum([E1["cs$g"]*u[3,3] for g in 4:n])
  -sum([(u[h,g]*E1["cs3"] +2*u[h,g]*E1["cs$g"]) for h in 3:n for g in 4:n])

-12*E1["cs2"]
+3*E1["cs3"]
-13*sum([E1["cs$i"] for i in 4:n])

+4* sum([E1["rs2"]*u[i,3] for i in 2:n if i != 3 && i != 2])
+3* sum([E1["rs3"]*u[i,3] for i in 2:n if i != 3 ])
+sum([(3*sum([E1["rs$g"]*u[i,3] for i in 2:n if i != 3 && i != g]) -E1["rs$g"]*u[3,3]) for g in 4:n])

+3*u[2,3]*E1["rs3"]
-3*sum([u[2,3]*E1["rs$i"] for i in 4:n])
+3*sum([u[2,i]*E1["rs3"] for i in 4:n])

-sum(u[2,3]*E1["rs$i"] for i in 4:n)

-6 * sum(u[h,g]*E1["rs$i"] for h in 3:n for g in 3:n for i in 2:n if h != i && i != 3)

-4*sum([u[g,h]*E1["rs$g"] for h in 3:n for g in 4:n])
-4*sum([u[2,h]*E1["rs$g"] for h in 4:n for g in 4:n])
+sum([u[g,h]*E1["rs3"] for g in 4:n for h in 3:n])

+12*E1["rs2"]
+13*sum([E1["rs$i"] for i in 4:n])
-3*E1["rs3"]

#rinjes
+sum([(-sum([u[2,h]*E1["rinj$(i)3"] for i in 4:n])-sum([u[g,h]*E1["rinj$(i)3"] for g in 3:n for i in 2:n if i != 3 && i != g])+sum([u[g,h]*E1["rinj$k$i"] for g in 3:n for k in 2:n  for i in 2:n if  k != g && i != k])) for h in 3:n])
-4*sum([E1["rinj2$i"] for i in 4:n])
-3*sum([E1["rinj$g$i"] for g in 3:n for i in 2:n if i != g])
+4*sum([E1["rinj$(i)3"] for i in 4:n])

-sum([E1["bg8_$g$h$k"] for g in 3:n for h in 2:n for k in 2:n if g != h && 3 != k && h != k])

-sum([E1["uurs2$g$(k)$(h)3"] for g in 3:n for k in 4:n for h in 3:n if h != g])
+sum([E1["uurs$(h)$j$g$k$i"] for j in 3:n for h in 3:n for k in 3:n for g in 2:n for i in 2:n if i != g && i != 3 && g != h && k != j])
+sum([E1["ucsu$(j)$(g)$(k)3$i"] for g in 3:n for j in 3:n for k in 2:n for i in 2:n if i != k && i != 3 && k != g])
-sum([E1["uucs$h$j$g$k$i"] for h in 3:n for j in 3:n for g in 3:n for k in 2:n for i in 2:n if i != k && k != j && i != 3 && g != h])

+sum([sum([E1["iesc$h$g$k$i"] for h in 3:n for k in 2:n for i in 2:n if i != k & k != h && i != 3])-sum([E1["iesc2$g$(k)3"] for k in 4:n]) for g in 3:n])


+sum([E1["inj23$k"]*u[3,3] for k in 4:n])
-sum([u[h,3]*E1["inj$(g)3$k"] for h in 3:n for g in 2:n for k in 2:n if k != g && k != 3 && g != h])

+sum([(E1["inj2$(h)3"]+2*sum([E1["inj2$(h)$k"] for k in 4:n]) +sum([2*sum([E1["inj$(g)$(h)$k"] for k in 2:n if k != 3 && k != g]) for g in 2:n ]) -sum([E1["inj$(g)$(h)3"] for g in 4:n])) for h in 4:n])

-sum([u[h,j]*E1["rwel$k$i"]  for j in 3:n for h in 3:n for k in 2:n for i in 2:n if k != j && k != i && i != 3]) 

+4*sum([E1["rwel2$i"] for i in 2:n if 2 != i && i != 3])
+3*sum([E1["rwel$k$i"] for k in 3:n for i in 2:n if k != i && i != 3])
- sum([E1["rwel$(h)3"] for h in 4:n])

)
,len=40)
;words=["wesc"])



+E1["bg2sr32"]
;len=30)
,to_file="../data/reduction_strings/n_6_bege8243.txt")


