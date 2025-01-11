using Oscar
n = 7
G1 = g1_named(n)
u = magic_unitary(n)

G0 = g0_named(n)

#constructed with (n=8)
s = 5
t = 6
r = 7
w = 5

#reduction_string(G1,bg(8,s,t,2,u=u)*u[r,4]*u[v,3]-u[2,s+u[2,t]*u[4,s]*rinj(3,r,u=u))

reduction_string(G1, bg(8,s,t,2,u=u)*u[r,4]*u[w,3]-u[2,s]*u[4,t]*bg(2,3,r,w,u=u))
+sum(bg(8,s,t,2,u=u)*u[r,i]*u[t,3] for i in (5:n))
-sum(bg(8,s,t,a,u=u)*u[r,i]*u[t,j] for a in (3:n) for i in (2:n) for j in (4:n)) ## conditions.. 


