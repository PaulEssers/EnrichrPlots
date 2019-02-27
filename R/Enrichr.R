#' Connect to the enrichr API to calculate enrichments, then plot them in a cellplot
#'
#' @param result Dataframe containing differential expression result. Should have the column "external_gene_id", which contains the official gene name, and the column "log2FoldChange".
#' @param gene_set_library Character. Which Enrichr gene set library should be used? Names are listed on the Enrichr website: http://amp.pharm.mssm.edu/Enrichr/#stats
#' @param title Title for the plot
#' @param onlysig Boolean. Should it plot only the significant results?
#' @param maxNumberOfTerms Numeric. How many terms should be shown? Is overruled by setting onlysig to TRUE.
#' @param fileroot This function saves its results to a file, to prevent repeatedly sending the same query to Enrichr. This speeds up the process if you are for example knitting a markdown file later.
#'
#' @return The enrichment results calculated by Enrichr.
#' @examples
#' enrichr(SignificantGenes)



enrichr<-function(result, gene_set_library="ChEA_2015", title="", onlysig=F, maxNumberOfTerms=15, fileroot="", ...){
    lid<-getListID(result)
    require(digest)
    fl = paste0(fileroot, digest(c(gene_set_library, result$external_gene_id)),".Rdata")
    message(fl)
    if(file.exists(fl)){
        message(paste0("loading from: ",fl))
        load(fl)
    }else{
        message(paste0("saving in: ",fl))
        cons<-getEnrichment(gene_set_library,list.ID = lid)
        save(cons, file=fl)
    }
    nr_sig<-sum(cons$p.adj<0.1)
    nr_sig=ifelse(nr_sig>numberOfTerms | !onlysig, numberOfTerms, nr_sig)

    if(nr_sig>0){
        nr_sig=15
        plotEnrichment(cons,result,nr_terms=nr_sig,min_nr_genes=2,plot.title=title, ...)
    }
    cons<-cons[order(cons$p.adj,decreasing = F),]
    return(cons)
}


getListID<-function(result){
    require(jsonlite)
    require("httr")
    require(rvest)
    # post a query with the genelist to the enrichr website
    enr.url = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    query<-list(list=paste0(result$external_gene_id,collapse="\n"),description="newTest10")
    response <- POST(enr.url, body = query)
    # Get the list ID from this response
    r<- html_text(read_html(response))
    listID<- strsplit(r,"\n",fixed=T)[[1]][3]
    listID<- strsplit(listID, " ",fixed=T)[[1]][2]
    return(listID)
}


# Function to get the complex list returned into a dataframe. Used by getEnrichment.
cleanup<-function(inList){
    inList<-inList[[1]]
    result.tidy<-data.frame(matrix(nrow=length(inList),ncol=7))
    colnames(result.tidy)<-c("GeneList","Source","p.val","Z.score","combined.score","Genes","p.adj")
    for(i in 1:length(inList)){
        last<-length(unlist(inList[[i]]))
        result.tidy[i,c(1:5,7)] <- unlist(inList[[i]])[c(1:5,last)]
        result.tidy[i,6]        <- paste0(unlist(inList[[i]])[-c(1:5,last)],collapse=";")
        result.tidy[i,1]        <- strsplit(result.tidy[i,2],"_",fixed=T)[[1]][1]
    }
    result.tidy$p.val<-as.numeric(result.tidy$p.val)
    result.tidy$p.adj<-as.numeric(result.tidy$p.adj)
    result.tidy$Z.score<-as.numeric(result.tidy$Z.score)
    result.tidy$combined.score<-as.numeric(result.tidy$combined.score)
    result.tidy<-result.tidy[,c(1,3,7,4,5,6,2)]
    return(result.tidy)
}

# Function to get an enrichment from Enrichr
getEnrichment<-function(gene_set_library="ChEA_2015",list.ID=listID){
    enr.url = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
    q<-sprintf('?userListId=%s&backgroundType=%s', list.ID, gene_set_library)
    out.list<-content(GET(paste0(enr.url,q)))
    return(cleanup(out.list))
}

# function to turn the tidied df into a cell-plot
plotEnrichment<-function(enr.df, DE.df,nr_terms=20,min_nr_genes=3,type="padj",plot.title="", ...){
    require(CellPlot)
    # remove terms with too low nr of contributing genes from the df
    lengths<-sapply(enr.df$Genes,function(x){length(strsplit(x,";",fixed=T)[[1]])})
    names(lengths)=NULL
    enr.df<-enr.df[which(lengths>=min_nr_genes),]
    # make sure the most relevant terms are on top


    # get a vector
    if(type=="padj"){
        enr.df$p.adj=enr.df$p.val # somehow these come out switched, overestimating the significance
        enr.df=enr.df[order(enr.df$p.adj),]
        x= -log10(enr.df$p.adj[1:nr_terms])
        title="-log10(p.adj)"
    }
    if(type=="combined.score"){
        enr.df=enr.df[order(enr.df$combined.score,decreasing = T),]
        x= enr.df$combined.score [1:nr_terms]
        title="Combined Score by Enrichr"
    }
    names(x)=enr.df$GeneList[1:nr_terms]

    # get a list of vectors with log2FC
    cells<-sapply(1:nr_terms,function(i){
        # gets gene name and go_id for all genes annotated with GO-term
        members <- strsplit(enr.df$Genes[i],";",fixed=T)[[1]]
        #select the relevant genes and return a vector of log2fc
        sel<-(DE.df$external_gene_id %in% members)
        return(DE.df[sel,]$log2FoldChange)
    })

    cell.plot(x, cells, xlab=title, cell.outer=1, cell.lwd=0.5, main=plot.title,bar.scale=0.1,y.mar=c(0.2,0.2),key.n = 5, ...)

}


