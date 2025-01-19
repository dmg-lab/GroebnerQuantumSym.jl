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


reduction_string(G1,bg(8,s,t,2,u=u)*u[r,4]*u[v,3]-u[2,s+u[2,t]*u[4,s]*rinj(3,r,u=u))

print(
reduction_string(G1, bg(8,s,t,2,u=u)*u[r,4]*u[w,3]-u[2,s]*u[4,t]*bg(2,3,r,w,u=u)
+sum(bg(8,s,t,2,u=u)*u[r,i]*u[w,3] for i in (5:n))
-sum(bg(8,s,t,3,u=u)*u[r,2]*u[w,j] for j in (4:n)) 
-sum(bg(8,s,t,3,u=u)*u[r,i]*u[w,2] for i in (4:n)) 
-sum(bg(8,s,t,3,u=u)*u[r,i]*u[w,j] for j in (4:n) for i in (4:n) if j!=i) 
-sum(bg(8,s,t,a,u=u)*u[r,2]*u[w,j] for a in (4:n) for j in (4:n) if a!=t)  #j!=3
-sum(bg(8,s,t,a,u=u)*u[r,i]*u[w,2] for a in (4:n) for i in (3:n) if i!=a && a!=t)  #i!=a
-sum(bg(8,s,t,a,u=u)*u[r,i]*u[w,j] for a in (4:n) for i in (3:n) for j in (4:n) if i!=a && a !=t &&i!=j)  
-sum(u[2,s]*u[a,t]*u[3,c]*rinj(r,w,u=u) for a in (4:n) for c in (3:n) if c!=t)
+sum(u[2,s]*u[a,t]*u[3,c]*row_sum(r,u)*u[w,3] for a in (4:n) for c in (3:n))
-sum(u[2,s]*u[4,t]*u[3,t]*u[r,d]*row_sum(w,u) for d in (3:n))
-sum(u[2,s]*u[a,t]*u[3,c]*u[r,d]*row_sum(w,u) for a in (4:n) for c in (3:n) for d in (3:n) if c!=t )
+sum(u[2,s]*u[a,t]*u[3,c]*inj(r,j,w,u=u) for a in (4:n) for c in (3:n) for j in (4:n))
+sum(u[2,s]*u[a,t]*inj(3,j,r,u=u)*u[w,2] for a in (4:n) for j in (3:n) if j != t)
+sum(u[2,s]*u[a,t]*inj(3,3,r,u=u)*u[w,e] for a in (4:n) for e in (4:n))
+sum(u[2,s]*u[a,t]*inj(3,b,r,u=u)*u[w,e] for a in (4:n) for b in (4:n) for e in (4:n) if b!=e && b !=t)
-sum(u[2,s]*inj(a,t,3,u=u)*u[r,c]*u[w,3] for a in (4:n) for c in (2:n))
-sum(u[2,s]*inj(a,t,3,u=u)*u[r,c]*u[w,c] for a in (4:n) for c in (4:n))  ## THIS WOULD BE NEW 
+sum(u[2,s]*u[4,t]*u[3,t]*u[r,d]*row_sum(w,u) for d in (3:n))
-sum(u[2,s]*u[j,t]*bg(2,3,r,w,u=u) for j in (5:n))
-sum(u[3,s]*col_sum(t,u)*u[3,2]*u[r,j]*u[w,3] for j in (4:n))
+sum(u[3,s]*col_sum(t,u)*u[3,i]*u[r,j]*u[w,k] for i in (3:n) for j in (2:n) for k in (2:n) if i!=j && k!=3 && j!=k && i!=t)
;len=300
))








