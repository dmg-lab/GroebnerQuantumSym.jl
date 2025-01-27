#= Example Case
using QuantumGB

n = 7
t = 6
G1 = g1_named(n)
u = magic_unitary(n)
=#

free_vars::NamedTuple=(ts="""collect(5:n)""",bs="""[1,2]""")::NamedTuple
to_be_reduced=raw"""bg(2,2,3,t,u=u)"""
reduction= raw"""
(G1["rwel24"]*u[t,3]
+sum(G1["bg2_2$(i)$(t)"] for i in (4:n) if i!=t)
-sum(u[j,2]*G1["wel$(t)$(i)3"] for i in (4:n) for j in (2:n))
+sum(G1["rwel2$i"]*u[t,3] for i in (5:n))
+sum(u[j,k]*G1["cs$i"]*u[t,3] for i in (4:n) for j in (3:n) for k in (2:n))
+sum(G1["bg2_$j$i$t"] for i in (2:n) for j in (3:n) if i!=t && i!=j)
-sum(G1["wel$i$j$k"]*u[t,3] for i in (3:n) for j in (2:n) for k in (4:n) if j!=k && (i,j)!=(t,2))
-sum(G1["wel$(i)43"]*u[t,3] for i in (2:n) if (i,4)!=(t,2))
-sum(G1["wel$i$(j)2"]*u[t,3] for i in (3:n) for j in (2:n) if j!= 2 && (i,j)!=(t,2))
-sum(u[j,k]*G1["rs$i"]*u[t,3] for i in (2:n) for j in (2:n) for k in (3:n) if i!=t && i!=j)
+sum(G1["rwel3$i"]*u[t,3] for i in (2:n) if i!=3)
+sum(G1["rwel$k$i"]*u[t,3] for i in (2:n) for k in (5:n) if i!=k)
+sum(G1["inj$k$j$i"]*u[t,3] for i in (2:n) for k in (2:n) for j in (3:n) if i!=k)
-sum(u[k,j]*G1["wel$t$(i)3"] for i in (2:n) for j in (3:n) for k in (2:n) if i!=3 && k !=t)
-sum(u[k,3]*G1["ip$(t)3"] for k in (2:n) if k != t)
-sum(u[k,j]*G1["ip$(t)3"] for k in (2:n) for j in (5:n) )
+sum(u[i,j]*G1["cs2"]*u[t,3] for i in (3:n) for j in (2:n))
+sum(u[i,j]*G1["cs3"]*u[t,3] for i in (3:n) for j in (5:n))
-sum(u[i,k]*G1["cs$k"]*u[t,3] for i in (3:n) for k in (2:n) if k!=4 && k!=3)
-sum(u[i,4]*G1["cs4"]*u[t,3] for i in (3:n))
+G1["rwel42"]*u[t,3]
+sum(G1["rwel4$i"]*u[t,3] for i in (5:n))
+sum(u[i,4]*G1["inj$(j)3$t"] for i in (2:n) for j in (2:n) if j!=t)
+sum(u[t,i]*G1["ip$(t)3"] for i in (4:n))
-sum(G1["wel$k$(j)3"]*u[t,3] for k in (3:n) for j in (5:n))
+(n-2)*sum(G1["rs$i"]*u[t,3] for i in (2:n) if i != t)
-(n-2)*sum(G1["cs$i"]*u[t,3] for i in (2:n))
+2*G1["cs3"]*u[t,3]
-2*sum(G1["inj$(i)3$t"] for i in (2:n) if i!= t)
+(n-2)*sum(G1["wel$t$(i)3"] for i in (2:n) if i !=3)
-G1["wel$(t)23"]
+(n-4)*G1["ip$(t)3"])
"""
