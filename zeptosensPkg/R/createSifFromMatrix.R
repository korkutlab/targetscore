#' Create sif from network matrix.
#' 
#' @param t.net p*p Gaussian gene Network matrix estimated.
#' @param genelist gene name list corresponding to the Gaussian gene Network estimated.
#' @return edgelist decription dataframe contain G(V,E) including Vertex and Edges(1 positive/-1 negative) edgevalue and vertex number of Gaussian Graphical Model.
#' @concept zeptosensPkg
#' @export
createSifFromMatrix=function(t.net,genelist){
    index=genelist
    edge.number=sum(t.net!=0)
    node1=array(0,dim=c(edge.number,1))
    node1_n=array(0,dim=c(edge.number,1))
    node2=array(0,dim=c(edge.number,1))
    node2_n=array(0,dim=c(edge.number,1))
    edges=array(0,dim=c(edge.number,1))
    edges.value=array(0,dim=c(edge.number,1))
    a=1
    for(i in 1:nrow(t.net)){
        for(j in 1:ncol(t.net)){
            if(t.net[i,j]!=0){
                
                node1[a]=index[i]
                node1_n[a]=i
                node2[a]=index[j]
                node2_n[a]=j
                edges.value[a]=t.net[i,j]
                edges[a]=ifelse(t.net[i,j]>0,1,-1)
                a=a+1
            }
        }
    }
    edgelist=as.data.frame(cbind(node1,edges,node2,edges.value,node1_n,node2_n))
    colnames(edgelist)=c("node1","edges","node2","edge.value","node1_n","node2_n")
    return(edgelist)
}

