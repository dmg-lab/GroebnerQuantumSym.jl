using Oscar
n = 8
G1 = g1_named(n)
u = magic_unitary(n)

G0 = g0_named(n)

#constructed with (n=8)
s = 5
t = 6
r = 7
w = 8


r_all = bg(8,s,t,2,u=u)*u[r,4]*u[w,3]-u[2,s]*u[4,t]*bg(2,3,r,w,u=u);

r_4 = r_all + (sum(bg(8,s,t,2,u=u)*u[r,i]*u[w,3] for i in (5:n))
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
+sum(u[x,s]*u[a,t]*bg(2,c,r,w,u=u) for x in (3:n) for a in (2:n) for c in (2:n) if (a,x)!=(3,3) && c!=r && c!=3 && c!=a) 
-sum(u[x,s]*col_sum(t,u)*u[3,2]*u[r,j]*u[w,3] for x in (3:n) for j in (4:n))
+sum(u[x,s]*col_sum(t,u)*u[3,i]*u[r,j]*u[w,k] for x in (3:n) for i in (3:n) for j in (2:n) for k in (2:n) if i!=j && k!=3 && j!=k && i!=t)
+sum(u[x,s]*u[b,t]*col_sum(2,u)*u[r,j]*u[w,3] for x in (3:n) for b in (3:n) for j in (4:n))
-sum(u[x,s]*u[b,t]*col_sum(i,u)*u[r,j]*u[w,d] for x in (3:n) for b in (3:n) for i in (3:n) for j in (2:n) for d in (2:n) if j!=i && d!=3 && d!=j && i!=t)
+sum(u[x,s]*rwel(t,2,u=u)*u[r,b]*u[w,3] for x in (3:n) for b in (4:n))
-sum(u[x,s]*rwel(t,a,u=u)*u[r,b]*u[w,c] for x in (3:n) for a in (3:n) for b in (2:n) for c in (2:n) if a!=t && b!=c && c!=3 && a!=b)
-sum(u[x,s]*u[a,t]*u[c,b]*row_sum(r,u)*u[w,3] for x in (3:n) for a in (2:n) for b in (3:n) for c in (2:n) if c != r && (x,a)!=(3,3) && (a,c)!=(2,2) &&c!=3 && a!=c)
+sum(u[x,s]*u[a,t]*u[c,j]*u[r,d]*row_sum(w,u) for x in (3:n) for a in (2:n) for c in (2:n) for d in (3:n) for j in (3:n) if (a,x)!=(3,3) && c != r && j!=t && c!=3 && (a,c)!=(2,2))
+sum(u[x,s]*u[a,t]*u[c,j]*rinj(r,w,u=u) for x in (3:n) for a in (2:n) for c in (2:n) for j in (3:n) if (x,a)!=(3,3) && c!=3 && (a,c)!=(2,2) && j!=t && c != r) # it could be c!=n
+sum(u[x,s]*u[a,t]*ip(r,j,u=u)*u[w,j] for x in (3:n) for a in (2:n) for j in (4:n) if j!=t && (x,a)!=(3,3))
-sum(u[x,s]*u[a,t]*ip(r,t,u=u)*u[w,3] for x in (3:n) for a in (2:n) if (x,a)!=(3,3))
+sum(u[x,s]*ip(a,t,u=u)*u[r,j]*u[w,j] for x in (3:n) for a in (3:n) for j in (4:n) if (x,a)!=(3,3) && a!=3)
+sum(u[x,s]*ip(r,t,u=u)*u[r,j]*u[w,3] for x in (3:n) for j in (2:n))
-sum(bg(8,s,t,a,u=u)*u[r,b]*u[w,c] for a in (3:n) for b in (2:n) for c in (2:n) if a != t && c!=b && a!=b && c!=3)
-sum(u[x,s]*u[a,t]*inj(b,j,r,u=u)*u[w,c] for x in (3:n) for a in (2:n) for b in (2:n) for j in (3:n) for c in (2:n) if (x,a)!=(3,3) && c != 3 && c!=j && b!=r && j!=t && b!=3 && (a,b)!=(2,2))
-sum(u[x,s]*u[a,t]*u[b,j]*inj(r,c,w,u=u) for x in (3:n) for a in (2:n) for j in (3:n) for b in (2:n) for c in (4:n) if (x,a)!=(3,3) && b!=3 && (a,b)!=(2,2))
+sum(u[x,s]*inj(b,t,a,u=u)*u[r,c]*u[w,3] for x in (3:n) for b in (2:n) for a in (2:n) for c in (2:n) if (x,b)!=(3,3) && b!=a && a!=3)
+sum(u[x,s]*inj(b,t,a,u=u)*u[r,c]*u[w,c] for x in (3:n) for b in (2:n) for a in (2:n) for c in (4:n) if (x,b)!=(3,3) && b!=a && a!=3)
-sum(wel(3,s,t,u=u)*u[b,2]*u[r,c]*u[w,3] for b in (2:n) for c in (4:n))
+sum(wel(3,s,t,u=u)*u[3,2]*u[r,c]*u[w,3] for c in (4:n))
+sum(wel(3,s,t,u=u)*u[b,i]*u[r,c]*u[w,d] for b in (2:n) for i in (3:n) for c in (2:n) for d in (2:n) if  d!=3 && c!=d && i!=t && b!=3 && i!=c)
-sum(u[x,s]*wel(b,t,j,u=u)*u[r,c]*u[w,3] for x in (3:n) for j in (2:n) for b in (4:n) for c in (2:n) if j!=t && (j,c)!=(2,2) && b!=r && (j,c)!=(2,3))
-sum(u[x,s]*u[b,t]*wel(r,2,c,u=u)*u[w,3] for x in (3:n) for c in (4:n) for b in (2:n) if (b,x)!=(3,3))
+sum(u[x,s]*u[b,t]*wel(r,a,c,u=u)*u[w,d] for x in (3:n) for b in (2:n) for a in (3:n) for c in (2:n) for d in (2:n) if (b,x)!=(3,3) && a!=c && d!=3 && (c,d)!=(2,2) && a!=t)
-sum(u[x,s]*u[b,t]*wel(r,t,c,u=u)*u[w,3] for x in (3:n) for b in (2:n) for c in (2:n) if (x,b)!=(3,3) && t!=c));

# to graded 3 component: 

r_4_rinj = r_4 +
(-sum(u[a,s]*u[b,j]*rinj(r,w,u=u) for a in (2:n) for b in (2:n) for j in (3:n) if b!=r && b!=3 && j!=t && (a,b)!=(2,2))
-sum(u[a,s]*u[b,j]*rinj(r,w,u=u) for a in (3:n) for b in (2:n) for j in (3:n) if b!=r && b!=3 && j!=t && (a,b)!=(2,2))
+sum(u[2,s]*u[b,t]*rinj(3,w,u=u) for b in (4:n))
-sum(u[a,s]*u[b,t]*rinj(c,w,u=u) for a in (3:n) for b in (2:n) for c in (2:n) if b!=a && c!=b && c!=3 && c!=w)
-(n-3)*sum(u[a,s]*u[b,t]*rinj(r,w,u=u) for a in (3:n) for b in (2:n) if b!=3)
-sum(u[a,s]*u[3,b]*rinj(r,w,u=u) for a in (4:n) for b in (3:n) if b!=s)
-(n-4)*sum(u[a,s]*u[3,t]*rinj(r,w,u=u) for a in (4:n))
-sum(u[a,t]*u[b,j]*rinj(r,w,u=u) for a in (2:n) for b in (2:n) for j in (3:n) if b!=r && b!=3 && j!=t && (a,b)!=(2,2)));

r_4_bg = r_4_rinj +
(-sum(u[a,s]*bg(2,i,r,w,u=u) for a in (2:n) for i in (2:n) if i != r && i!=3 && (a,i)!=(2,2))
-sum(u[a,s]*bg(2,i,r,w,u=u) for a in (3:n) for i in (2:n) if i != r)
-sum(u[a,t]*bg(2,i,r,w,u=u) for a in (2:n) for i in (2:n) if i != r && i!=3 && (a,i)!=(2,2))
+sum(bg(8,s,t,i,u=u)*u[w,d] for i in (2:n) for d in (2:n) if d!=3 && (i,d)!=(2,2) && i!=t)
+sum(bg(8,s,t,i,u=u)*u[w,d] for i in (3:n) for d in (2:n) if i!=t)
+sum(bg(8,s,t,i,u=u)*u[r,d] for i in (3:n) for d in (2:n) if i!=t));

r_4_rwel = r_4_bg +
(-sum(rwel(s,2,u=u)*u[r,a]*u[w,3] for a in (4:n))
+sum(rwel(s,a,u=u)*u[r,j]*u[w,d] for a in (3:n) for j in (2:n) for d in (2:n) if j!=d && (j,d)!=(2,2) &&  d!=3 && a!=t && s!=a) 
+sum(u[a,s]*rwel(t,c,u=u)*u[w,d] for a in (3:n) for c in (2:n) for d in (2:n) if d!=3 && (d,c)!=(2,2) && t!=c)
+sum(u[a,s]*rwel(t,c,u=u)*u[w,d] for a in (3:n) for c in (3:n) for d in (2:n) if t!=c)
+sum(u[a,s]*rwel(t,c,u=u)*u[r,d] for a in (3:n) for c in (3:n) for d in (2:n) if t!=c)
-sum(rwel(t,2,u=u)*u[r,a]*u[w,3] for a in (4:n))
+sum(rwel(t,a,u=u)*u[r,b]*u[w,d] for a in (3:n) for b in (2:n) for d in (2:n) if a!=t && d!=3 && (b,d)!=(2,2)));

print(
reduction_string(G1, r_4_rwel + 
(sum(u[2,s]*wel(r,2,j,u=u)*u[w,3] for j in (4:n))
-sum(u[2,s]*wel(r,b,j,u=u)*u[w,d] for b in (3:n) for j in (2:n) for d in (2:n) if j!=b&& b!=s && d!=3 && j!=d)
+sum(u[2,s]*wel(r,t,j,u=u)*u[w,d] for j in (2:n) for d in (2:n) if j!=d && j!=t && (j,d)!=(2,2))
+2*sum(u[3,s]*wel(r,2,j,u=u)*u[w,3] for j in (4:n))
-2*sum(u[3,s]*wel(r,3,j,u=u)*u[w,d] for j in (4:n) for d in (2:n) if d!=3 && d!=j)
-sum(wel(3,s,3,u=u)*u[r,c]*u[w,d] for c in (2:n) for d in (2:n) if c!=d)
);
words=["wel"],
len=1000
)
)






















