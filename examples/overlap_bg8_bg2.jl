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


print(
reduction_string(G1, bg(8,s,t,2,u=u)*u[r,4]*u[w,3]-u[2,s]*u[4,t]*bg(2,3,r,w,u=u)
+sum(bg(8,s,t,2,u=u)*u[r,i]*u[w,3] for i in (5:n))
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
+sum(u[3,s]*u[a,t]*bg(2,c,r,w,u=u) for a in (2:n) for c in (2:n) if a!=3 && c!=r && c!=3 && c!=a)
-sum(u[x,s]*col_sum(t,u)*u[3,2]*u[r,j]*u[w,3] for x in (3:n) for j in (4:n))
+sum(u[x,s]*col_sum(t,u)*u[3,i]*u[r,j]*u[w,k] for x in (3:n) for i in (3:n) for j in (2:n) for k in (2:n) if i!=j && k!=3 && j!=k && i!=t)
+sum(u[x,s]*u[b,t]*col_sum(2,u)*u[r,j]*u[w,3] for x in (3:n) for b in (3:n) for j in (4:n))
-sum(u[x,s]*u[b,t]*col_sum(i,u)*u[r,j]*u[w,d] for x in (3:n) for b in (3:n) for i in (3:n) for j in (2:n) for d in (2:n) if j!=i && d!=3 && d!=j && i!=t)
+sum(u[3,s]*rwel(t,2,u=u)*u[r,b]*u[w,3] for b in (4:n))
-sum(u[3,s]*rwel(t,a,u=u)*u[r,b]*u[w,c] for a in (3:n) for b in (2:n) for c in (2:n) if a!=t && b!=c && c!=3 && a!=b)
-sum(u[3,s]*u[a,t]*u[c,b]*row_sum(r,u)*u[w,3] for a in (2:n) for b in (3:n) for c in (2:n) if c != r && a!=3 && (a,c)!=(2,2) &&c!=3 && a!=c)
+sum(u[3,s]*u[a,t]*u[c,j]*u[r,d]*row_sum(w,u) for a in (2:n) for c in (2:n) for d in (3:n) for j in (3:n) if a!=3 && c != r && j!=t && c!=3 && (a,c)!=(2,2))
+sum(u[3,s]*u[a,t]*u[c,j]*rinj(r,w,u=u) for a in (2:n) for c in (2:n) for j in (3:n) if a!=3 && c!=3 && (a,c)!=(2,2) && j!=t && c != r) # it could be c!=n
+sum(u[3,s]*u[a,t]*ip(r,j,u=u)*u[w,j] for a in (2:n) for j in (4:n) if j!=t && a!=3)
-sum(u[3,s]*u[a,t]*ip(r,t,u=u)*u[w,3] for a in (2:n) if a!=3)
+sum(u[3,s]*ip(a,t,u=u)*u[r,j]*u[w,j] for a in (4:n) for j in (4:n))
+sum(u[3,s]*ip(r,t,u=u)*u[r,j]*u[w,3] for j in (2:n))
-sum(bg(8,s,t,a,u=u)*u[r,b]*u[w,c] for a in (3:n) for b in (2:n) for c in (2:n) if a != t && c!=b && a!=b && c!=3)
-sum(u[3,s]*u[a,t]*inj(b,j,r,u=u)*u[w,c] for a in (2:n) for b in (2:n) for j in (3:n) for c in (2:n) if a!=3 && c != 3 && c!=j && b!=r && j!=t && b!=3 && (a,b)!=(2,2))
-sum(u[3,s]*u[a,t]*u[b,j]*inj(r,c,w,u=u) for a in (2:n) for j in (3:n) for b in (2:n) for c in (4:n) if a!=3 && b!=3 && (a,b)!=(2,2))
+sum(u[3,s]*inj(b,t,a,u=u)*u[r,c]*u[w,3] for b in (2:n) for a in (2:n) for c in (2:n) if b!=3 && b!=a && a!=3)
+sum(u[3,s]*inj(b,t,a,u=u)*u[r,c]*u[w,c] for b in (2:n) for a in (2:n) for c in (4:n) if b!=3 && b!=a && a!=3)
-sum(wel(3,s,t,u=u)*u[b,2]*u[r,c]*u[w,3] for b in (2:n) for c in (4:n))
+sum(wel(3,s,t,u=u)*u[3,2]*u[r,c]*u[w,3] for c in (4:n))
+sum(wel(3,s,t,u=u)*u[b,i]*u[r,c]*u[w,d] for b in (2:n) for i in (3:n) for c in (2:n) for d in (2:n) if  d!=3 && c!=d && i!=t && b!=3 && i!=c)
-sum(u[3,s]*wel(b,t,j,u=u)*u[r,c]*u[w,3] for j in (2:n) for b in (4:n) for c in (2:n) if j!=t && (j,c)!=(2,2) && b!=r && (j,c)!=(2,3))
-sum(u[3,s]*u[b,t]*wel(r,2,c,u=u)*u[w,3] for c in (4:n) for b in (2:n) if b!=3)
+sum(u[3,s]*u[b,t]*wel(r,a,c,u=u)*u[w,d] for b in (2:n) for a in (3:n) for c in (2:n) for d in (2:n) if b!=3 && a!=c && d!=3 && (c,d)!=(2,2) && a!=t)
-sum(u[3,s]*u[b,t]*wel(r,t,c,u=u)*u[w,3] for b in (2:n) for c in (2:n) if b!=3 && t!=c)
;
words=["cs"],
len=500
))





