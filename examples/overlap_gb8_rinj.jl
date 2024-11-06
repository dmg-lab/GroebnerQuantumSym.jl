s = 5
t = 6
r = 7



reduction_string(G1,bg(8,t,s,2,u=u)*u[r,3]+u[2,t]*u[4,s]*rinj(3,r,u=u))




# u[i,a]*u[j,b]*u[k,c]*u[m,d]

# graded-4 component finished?
reduction_string(G1,u[s,2]*bg(8,3,t,r,u=u)+rinj(s,2,u=u)*u[4,t]*u[3,r]
+sum(rinj(s,2,u=u)*u[i,t]*u[3,r] for i in (5:n))
-sum(rinj(s,3,u=u)*u[i,t]*u[2,r] for i in (4:n))
-sum(rinj(s,j,u=u)*u[i,t]*u[m,r] for m in (4:n) for i in (2:n) for j in (3:n) if i != 3 && j != s)
-sum(rinj(s,j,u=u)*u[3,t]*u[m,r] for m in (2:n) for j in (4:n) if  m != 3 && j!= s)
-sum(rinj(s,j,u=u)*u[i,t]*u[2,r] for i in (4:n) for j in (4:n) if s!= j)
-sum(inj(s,a,b,u=u)*u[c,t]*u[3,r] for a in (3:n) for b in (2:n) for c in (2:n) if s!=b)
-sum(u[s,2]*inj(b,t,c,u=u)*u[3,r] for b in (2:n) for c in (2:n) if c != 3 && b != c)
-sum(u[s,a]*inj(b,t,c,u=u)*u[3,r] for a in (3:n) for b in (2:n) for c in (2:n) if c != 3 && b != c && t!=a)
-sum(u[s,a]*u[j,2]*col_sum(t,u)*u[3,r] for a in (3:n) for j in (3:n))
-sum(u[s,a]*u[j,b]*col_sum(t,u)*u[3,r] for a in (2:n) for j in (3:n) for b in (3:n) if b != t && a!=b)
+sum(u[s,a]*u[i,b]*u[j,t]*col_sum(r,u) for a in (2:n) for i in (3:n) for b in (2:n) for j in (3:n) if i!=s && (a,b)!=(2,2))
+sum(wel(s,a,b,u=u)*u[i,t]*u[3,r] for a in (2:n) for i in (2:n) for b in (2:n) if b != a)
-sum(u[s,2]*wel(3,b,t,u=u)*u[m,r] for b in (3:n) for m in (2:n) if b != t && m != 3)
-sum(u[s,a]*wel(3,b,t,u=u)*u[m,r] for a in (3:n) for b in (2:n) for m in (2:n) if b != t && m != 3)
+sum(u[s,a]*wel(2,a,t,u=u)*u[3,r] for a in (3:n))
-sum(u[s,a]*u[b,t]*wel(3,t,r,u=u) for a in (3:n) for b in (4:n))
+sum(u[s,a]*u[2,a]*wel(3,t,r,u=u) for a in (3:n))
+sum(u[s,t]*u[b,t]*wel(3,t,r,u=u) for b in (4:n) if b!=s)
-sum(u[s,2]*u[b,t]*wel(3,t,r,u=u) for b in (4:n))
+sum(u[s,2]*u[j,b]*rwel(t,r,u=u) for j in (3:n) for b in (3:n) if j!= s)
+sum(u[s,a]*u[j,b]*rwel(t,r,u=u) for a in (3:n) for j in (3:n) for b in (2:n) if j!= s)
+sum(u[s,2]*bg(8,j,t,r,u=u) for j in (4:n) if j!= t)
+sum(u[s,a]*bg(8,j,t,r,u=u) for a in (3:n) for j in (2:n) if j!=a && j!=t)
+ip(s,t,u=u)*u[3,t]*u[3,r]
-sum(u[s,a]*ip(3,t,u=u)*u[d,r] for a in (2:n) for d in (2:n))
+sum(u[s,t]*ip(b,t,u=u)*u[3,r] for b in (3:n))
-sum(u[s,a]*ip(b,t,u=u)*u[3,r] for a in (2:n) for b in (4:n))
+sum(u[s,a]*row_sum(2,u)*u[k,t]*u[3,r] for a in (3:n) for k in (4:n))
-sum(u[s,a]*row_sum(j,u)*u[k,t]*u[2,r] for a in (3:n) for j in (3:n) for k in (3:n) if j !=s && (j,k)!=(3,3))
-sum(u[s,a]*row_sum(j,u)*u[k,t]*u[m,r] for a in (3:n) for j in (3:n) for k in (2:n) for m in (4:n) if j != s && (j,k)!=(3,3)) 
)








