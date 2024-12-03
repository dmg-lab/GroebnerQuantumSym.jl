using Oscar

n = 8
G1 = g1_named(n)
u = magic_unitary(G1)

#Define helper:

bg2cs(i,j,u) = sum([G1["bg2_$(i)$(j)$(k)"] for k in 2:n if k!= j && k!=3]) + u[i,j]*G1["wel$(j)43"]

bg2css = [ bg2cs(i,j,u) for i in 3:n for j in 4:n if i != j]
bg2css_names = indexed_name("bg2cs",[parse(Int, "$j$i") for i in 3:n for j in 4:n if i!=j]) 
bg2css_ids = [Symbol("bg2cs$j$i") for i in 3:n for j in 4:n if i != j]      

E1 = named_generators(bg2css,bg2css_names,bg2css_ids)

QuantumGB.add!(E1, G1; check=false)

bg243 = bg(8,2,4,3,u=u)

reduction_string(G1, bg243
- sum([G1[Symbol("bg2_2$(i)3")] for i in 4:n])
+ sum([G1[Symbol("bg8_2$(i)3")] for i in 5:n])
- u[3,2] * G1[Symbol("cs4")]*u[3,3] + u[3,2]*G1[Symbol("rwel43")] #cs interaction with rwel
+ sum([G1[Symbol("bg2_32$i")] for i in 4:n])
+ u[3,2] * u[3,4]*G1[Symbol("cs3")] 
- sum([G1["wel324"] * u[k,3] for k in 2:n if k != 3])
+u[3,2] * u[4,4] * G1["cs3"]
+sum([G1["bg2_34$k"] for k in 2:n if k != 4 && k != 3])

)
