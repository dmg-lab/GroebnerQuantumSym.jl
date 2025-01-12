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

reduction_string(G1, bg(8,s,t,2,u=u)*u[r,4]*u[w,3]-u[2,s]*u[4,t]*bg(2,3,r,w,u=u)
+sum(bg(8,s,t,2,u=u)*u[r,i]*u[w,3] for i in (5:n))
-sum(bg(8,s,t,3,u=u)*u[r,2]*u[w,j] for j in (4:n)) 
-sum(bg(8,s,t,3,u=u)*u[r,i]*u[w,2] for i in (4:n)) 
-sum(bg(8,s,t,3,u=u)*u[r,i]*u[w,j] for j in (4:n) for i in (4:n) if j!=i) 
-sum(bg(8,s,t,a,u=u)*u[r,2]*u[w,j] for a in (4:n) for j in (4:n) if a!=t)  #j!=3
-sum(bg(8,s,t,a,u=u)*u[r,i]*u[w,2] for a in (4:n) for i in (3:n) if i!=a && a!=t)  #i!=a
-sum(bg(8,s,t,a,u=u)*u[r,i]*u[w,j] for a in (4:n) for i in (3:n) for j in (4:n) if i!=a && a !=t &&i!=j)  
#-sum(bg(8,s,t,a,u=u)*u[r,i]*u[w,j] for a in (4:n) for i in (3:n) for j in (3:n) if i!=a && a !=t && i!=j&& j!=t && a!=t) 
#-sum(bg(8,s,t,a,u=u)*u[r,i]*u[w,j] for a in (4:n) for j in (3:n) for i in (j+1:n) if j!=a && a!=t) 
)

