#question: distances from database or a SIF?
#question: Can (should) we create a new distances file for new data
#check pmatch vs. match in extracting the wks. MAke it all "match" and test. 
#number of rows returned in match for cbind.
# i <- inhibiting PTM.
# a <- activating PTM
#c <- total protein measurement (no specific PTM, includes everything)
#ts: integral_dose(fs*(xi+sigma_j(2^p*xj*product_k(wk))))

nodes: 
total protein level (c)
phospho level (a or i)

a -> b

wk=-1,0,1

Signed PC output
a upreg b (a(c,i,a) & b(c))
a downreg b 
a phos b
a dephosp b 

a upreg b (from Augustin)
wk=1
a(c) -> b(c)
a(a) -> b(c)

a downreg b
wk = -1
a(c) -> b(c)
a(a) -> b(c)

a phosp b
AKT phosph RB
AKTp308 (a) + oncogenic
RBp806 (i) - oncogenic
wk = 1

SRCp416 (a) oncogenic
SRCp517 (i) suppress

AKT phosp SRC

ATPp473 phosp SRCp416 wk=+1
ATPp473 phosp SRCp517 wk=0 unless there is explicit information

tumor supressor phosp oncogene
TS phosp O
TS(a)
O(a)

active form phosp active wk=1
active phosp inh wk=0



O(a) phosp O2
wk=1 O2(a)
wk=0 O2(i)









wk=1
a(a) -> b(a)
a(a) -> b(i) ? 
a(c) -> b(a) ?
a(i) -> XX 






get the data from WQ. melanoma with synergy/antogonism to MEKi/RAFi.

