library(dplyr)
library(ape)
library(phytools)

results <- data.frame(
  seqfile = character(),
  m0_kappa = numeric(),
  m0_beta_0 = numeric(),
  m0_beta_1 = numeric(),
  m0_treescaling = numeric(),
  m0_ml = numeric(),
  m1_kappa = numeric(),
  m1_beta_0 = numeric(),
  m1_beta_1 = numeric(),
  m1_treescaling = numeric(),
  m1_ml = numeric(),
  stringsAsFactors = FALSE
)



#This code extracts BODYWEIGHT values from the continuous data


#Reading in the FASTA files of interest, set path as location where FASTA files are
seqfiles = list.files(path = "/Users/emilydunican/Desktop/seqfiles_growth", full.names = TRUE)


for (seqfile in seqfiles) {
  fasta_lines = readLines(seqfile)
  
#Creating variables to allow for the storage of species name and their corresponding sequences 
  species_name = ""
  sequence = ""
  
#Lists to store species names and sequences
  species_names = character()
  sequences = character()
  
#Loop through the lines of the FASTA file, extracting species names - e.g. >emily_dunican will now be species, emily_dunican. 
#Identifying the species name by looking at the > in front of the name 
  for (i in 1:length(fasta_lines)) {
    line = fasta_lines[i]
    if (startsWith(line, ">")) { 
      # Extract species name 
      species_name = gsub("^>(.+)", "\\1", line)  
      
      if (nchar(sequence) > 0) {
        species_names <- c(species_names, species_name)
        sequences <- c(sequences, sequence)
        
        sequence <- ""
      }
    } else {  
      sequence <- paste0(sequence, line)
    }
  }
  
  species_names = c(species_names, species_name)
  sequences = c(sequences, sequence)
  fasta_df = data.frame(Species = species_names, Sequence = sequences, stringsAsFactors = FALSE)
  
#Read in the continuous data and extract the species names from this list
  continuous_data = read.table("/Users/emilydunican/Desktop/species_data.txt", header = TRUE)
  continuous_species = continuous_data$Species
  
#Extract species from fasta data
  gene_species = fasta_df$Species
  
#Find common species between fasta and continuous data and filter the continuous data based on common species 
  common_species = intersect(continuous_species, gene_species)
  common_continuous_data = continuous_data[continuous_data$Species %in% common_species,]
  
#Filter fasta data based on common species
  common_fasta_df = fasta_df[fasta_df$Species %in% common_species, ]
  
#Merge filtered continuous data with fasta data based on species and remove the adult weight data 
  merged_data = merge(common_continuous_data, common_fasta_df, by = "Species")
  merged_data$max_longevity.months. = NULL
  
  merged_data = merged_data[complete.cases(merged_data), ]
  
#editing the original tree by removing all of the non common species - needed to make a vector 
#called not_in_merged_species and remove that from the tree
  not_in_merged_species = setdiff(continuous_species, merged_data$Species)
  
  
#Filter the tree to the subset of species - using the not_in_merged_species vector to drop those species from the tree, 
#result is the subset_tree, a subset of the orginal tree that contains the overlapping species we have 
  
  phy = read.tree("/Users/emilydunican/Desktop/Mammals_tree_all_10.nwk")
  subset_tree = drop.tip(phy, not_in_merged_species, trim.internal = TRUE, subtree = FALSE,
                         root.edge = 0, rooted = is.rooted(phy), collapse.singles = TRUE,
                         interactive = FALSE)
  

#To attach data to tree, make matrix and then make the data numeric 
  bodyweight_data = as.matrix(merged_data$adult_weight.g.)
  bodyweight_data = as.numeric(bodyweight_data)
  
#add Species and add labels back
  bodyweight_data = setNames(bodyweight_data,merged_data$Species)
  
  
#Reconstruct ancestral states using fastAnc function 
  anc_states = fastAnc(subset_tree, bodyweight_data, vars=FALSE, CI=FALSE)
  anc_states_df = data.frame(node = as.numeric(names(anc_states)), reconstructed = as.numeric(anc_states))
  contMap(subset_tree,bodyweight_data)
  
#Reconstruct ancestral states (maximum likelihood model)
  anc = fastAnc(subset_tree, bodyweight_data)
  anc_df = data.frame(node = as.numeric(names(anc)), reconstructed = as.numeric(anc))
  
  
##Create dataframe of branches, this will be used for the BODYSIZE vector 
  branches = subset_tree$edge
  branches = data.frame(branches)
  colnames(branches) = c("node1","node2")
  
  
#Add reconstructed values to the dataframe
  branches = left_join(branches, anc_states_df, by = c("node1"="node"))
  colnames(branches) = c("node1","node2","reconstructed_node1")
  branches = left_join(branches, anc_df, by = c("node2"="node"))
  colnames(branches) = c("node1","node2","reconstructed_node1","reconstructed_node2")
  
#Reconstructed_node2 column contains some NAs because these are the terminal branches
#map tip labels to node numbers to longevity
  tip_labels = subset_tree$tip.label
  tip_label_to_node = data.frame(subset_tree$tip.label, 1:length(subset_tree$tip.label))
  colnames(tip_label_to_node) = c("Species","node")
  tip_label_to_node = left_join(tip_label_to_node, merged_data[,c("Species","adult_weight.g.")])
  
#Fill in terminal node values
  branches = left_join(branches, tip_label_to_node, by = c("node2"="node")) %>% dplyr::select(-Species)
  branches = branches %>%
    mutate(reconstructed_node2 = ifelse(is.na(reconstructed_node2), adult_weight.g., reconstructed_node2)) %>%
    dplyr::select(-adult_weight.g.)
  branches$mean = (branches$reconstructed_node1 + branches$reconstructed_node2) / 2
  
  
#Plot phylogenetic tree with ancestral states
  plot.phylo(subset_tree, show.tip.label = TRUE)
  nodelabels(anc_states, frame = "rect", cex = 0.7)
  
 
  treefile = subset_tree  
  
#Creating a .csv file that contains the average values of the bodyweight, which will be read in 
  write.csv(branches, file = "branches.csv")
  
  
  #!/usr/bin/env Rscript
  
  library('ape')
  
  
  
##Calculate the vector bodyweight, each entry of which is the average bodyweight for the two nodes at either end of each branch
  branches = read.csv("branches.csv") 
  bodyweight = branches$mean
  
  
  ## Usage: Rscript codon_ml seqfilesdataset - to run type Rscript codon_ml seqfilesdataset in the command line. 
  ## This dataset should contain the FASTA files of genes of interest.
  ## The tree file should include branch lengths. Note that relative branch lengths are treated as fixed.
  ## The sequence file should contain a codon-aware alignment in fasta format and must include the 
  ## stop codon as the last position of the alignment. Sequences that include a gap or a codon other than
  ## a stop codon (i.e. sequences for which the stop codon is not positionally homologous with the last 
  ## position in the alignment) are excluded.
  ## Output: a file called stopcodon.rscript.out, which contains the following: The maximum log likelihood;
  ## ML estimate of kappa; ML estimate of omega; ML estimate of the treescaling parameter; 
  ## convergence of optimizer (0 = success,1 = failure)
  
  args = commandArgs(trailingOnly=TRUE)
  
  model = 'MG'
  codon_frequencies_model = 'f1x4'
  
  #####Requires ape and expm
  require(ape)
  require(expm)
  
  #tree_file = args[1]
  #seq_file = args[2]
  out_file = paste(seqfile,".stopout",sep="")
  
  ######Setup########
  setup = function(beta_0,beta_1,kappa,treefile, seqfile) {
    amino_acids = c("F", "L", "S", "Y", "C", "W", "P", "H", "Q", "R", "I", "M", "T", "K", "N", "V", "A", "D", "E", "G","X")
    codons = c("ttt", "ttc", "tta", "ttg", "tct", "tcc", "tca", "tcg", "tat", "tac", "tgt", "tgc", "tgg", "ctt", "ctc", "cta", "ctg", "cct", "ccc", "cca", "ccg", "cat", "cac","caa", "cag", "cgt", "cgc", "cga", "cgg", "att", "atc", "ata", "atg", "act", "acc", "aca", "acg", "aat", "aac","aaa", "aag", "agt", "agc", "aga", "agg", "gtt", "gtc","gta", "gtg", "gct", "gcc", "gca", "gcg", "gat", "gac","gaa", "gag", "ggt", "ggc", "gga", "ggg", "tag", "tga","taa")
    aa = c("F", "F", "L", "L", "S", "S", "S", "S", "Y", "Y","C", "C", "W", "L", "L", "L", "L", "P", "P", "P", "P","H", "H", "Q", "Q", "R", "R", "R", "R", "I", "I", "I","M", "T", "T", "T", "T", "N", "N", "K", "K", "S", "S","R", "R", "V", "V", "V", "V", "A", "A", "A", "A", "D","D", "E", "E", "G", "G", "G", "G")
    codon_numbers = list(ttt = 1, ttc = 2, tta = 3, ttg = 4,tct = 5, tcc = 6, tca = 7, tcg = 8, tat = 9, tac = 10,tgt = 11, tgc = 12, tgg = 13, ctt = 14, ctc = 15, cta = 16,ctg = 17, cct = 18, ccc = 19, cca = 20, ccg = 21, cat = 22,cac = 23, caa = 24, cag = 25, cgt = 26, cgc = 27, cga = 28,cgg = 29, att = 30, atc = 31, ata = 32, atg = 33, act = 34,acc = 35, aca = 36, acg = 37, aat = 38, aac = 39, aaa = 40,aag = 41, agt = 42, agc = 43, aga = 44, agg = 45, gtt = 46,gtc = 47, gta = 48, gtg = 49, gct = 50, gcc = 51, gca = 52,gcg = 53, gat = 54, gac = 55, gaa = 56, gag = 57, ggt = 58,ggc = 59, gga = 60, ggg = 61)
    nucleotides = list(a = 1, c = 2, g = 3, t = 4)
    purine = c(1, 0, 1, 0)
    
    tree = treefile
    tree=reorder(tree,"postorder")
    edges = tree$edge
    
    seqs = read.dna(seqfile,fo="fasta")
    chars = as.character(seqs)
    chars = chars[tree$tip.label,]
    
    #set up PI (codon stationary frequencies)
    f1 = c(0, 0, 0, 0)
    f2 = c(0, 0, 0, 0)
    f3 = c(0, 0, 0, 0)
    ff = c(0,0,0,0)
    PI = array(rep(0, 61*61), dim = c(61, 61))
    PI3x4 = array(rep(0, 61*61), dim = c(61, 61))
    PI1x4 = array(rep(0, 61*61), dim = c(61, 61))
    
    #convert to codons - produces cdAlign which is a matrix representing the coding sequence alignment. Each
    #row is a sequence, each column is a numeric representation (1-61) of the codon
    cdAlign = array(rep(0,nrow(chars)*ncol(chars)/3),c(nrow(chars),ncol(chars)/3))
    for(i in 1:nrow(chars)) {
      s = paste(chars[i,],collapse="")
      for(j in 1:(ncol(chars)/3)) {
        ss = substr(s,3*(j-1)+1,3*j)
        if(!is.null(codon_numbers[[ss]])) {
          cdAlign[i,j] = codon_numbers[[ss]]
          f1[nucleotides[[substr(ss,1,1)]]] = f1[nucleotides[[substr(ss,1,1)]]] + 1
          f2[nucleotides[[substr(ss,2,2)]]] = f2[nucleotides[[substr(ss,2,2)]]] + 1
          f3[nucleotides[[substr(ss,3,3)]]] = f3[nucleotides[[substr(ss,3,3)]]] + 1
          ff[nucleotides[[substr(ss,1,1)]]] = ff[nucleotides[[substr(ss,1,1)]]] + 1
          ff[nucleotides[[substr(ss,2,2)]]] = ff[nucleotides[[substr(ss,2,2)]]] + 1
          ff[nucleotides[[substr(ss,3,3)]]] = ff[nucleotides[[substr(ss,3,3)]]] + 1
        }
        else {
          cdAlign[i,j] = NA
        }
      }
    }
    dimnames(cdAlign)[[1]]=tree$tip.label
    f1 = f1/sum(f1)
    f2 = f2/sum(f2)
    f3 = f3/sum(f3)
    ff = ff/sum(ff)
    
    for (i in 1:61) {
      PI3x4[i, i] = f1[nucleotides[[substr(codons[i], 1, 1)]]] * f2[nucleotides[[substr(codons[i],2, 2)]]]* f3[nucleotides[[substr(codons[i],3, 3)]]]
      PI1x4[i, i] = ff[nucleotides[[substr(codons[i], 1, 1)]]] * ff[nucleotides[[substr(codons[i],2, 2)]]]* ff[nucleotides[[substr(codons[i],3, 3)]]]
    }
    
    # Toggle f3x4 f1x4
    #PI = PI3x4/sum(PI3x4)
    PI = PI1x4/sum(PI1x4)
    
    Mg_mult = array(c(rep(0, 61*61)), dim = c(61, 61))
    
    ###############Set up the generator matrix############
    
    
    ### Changed to allow Omega to depend on model 
    Omega = rep(0,nrow(edges))
    for(i in 1:nrow(edges)) {
      Omega[i] = beta_0 + beta_1 * bodyweight[i]
    }
    
    Nonsyn_matrix = array(c(rep(0,61*61)),dim=c(61,61)) # Matrix indicating nonsynonymous codon pairs
    R = array(rep(0, nrow(edges)*61*61), dim = c(nrow(edges),61, 61))
    
    for (i in 1:61) {
      for (j in 1:61) {
        diffs = 0
        for (k in 1:3) {
          nuc1 = nucleotides[[substr(codons[i], k, k)]]
          nuc2 = nucleotides[[substr(codons[j], k, k)]]
          if (nuc1 != nuc2) {
            diffs = diffs + 1
            if(codon_frequencies_model == 'f1x4') {
              Mg_mult[i,j] = ff[nuc2]
            } else if(codon_frequencies_model == 'f3x4') {
              if(k==1) {Mg_mult[kk,i,j] = f1[nuc2]}
              if(k==2) {Mg_mult[kk,i,j] = f2[nuc2]}
              if(k==3) {Mg_mult[kk,i,j] = f3[nuc2]}
            } else {print("Codon frequencies model undefined")
              return(NA)
            }
            if (purine[nuc1] == purine[nuc2]) {
              R[,i, j] = kappa
            }
            else {
              R[,i,j] = 1
            }
            if (aa[i] != aa[j]) {
              Nonsyn_matrix[i,j] = R[1,i,j]
              
              R[,i,j] = 0
            }
          }
          if (diffs > 1) {
            R[,i, j] = 0
            Nonsyn_matrix[i,j] = 0
          }
        }
      }
    }
    
    for(k in 1:nrow(edges)) { 
      R[k,,] = R[k,,] + Nonsyn_matrix * Omega[k]
    }
    
    
    if(model == 'MG') {
      for(k in 1:nrow(edges)) { 
        R[k,,] = R[k,,] * Mg_mult
      }
    } else if(model == 'GY') {
      for(k in 1:nrow(edges)) { 
        R[k,,] = R[k,,] %*% PI
      }
    } else {print("Model undefined")
      return(NA)
    }
    
    #Scale the generator matrix
    #scale_fac = sum(PI%*%R)
    #Scale the generator matrix (using only the sense codon part, since this is the part that influences the branch length estimates
    #scale_fac = sum(PI[1:61,1:61]%*%R[1:61,1:61])
    for(k in 1:nrow(edges)) { 
      scale_fac = sum(PI%*%R[k,,])
      R[k,,] = (1/scale_fac)*R[k,,]
    }
    
    
    for(i in 1:61) {
      for(k in 1:nrow(edges)) { 
        R[k,i,i] = -sum(R[k,i,-i])
      }
    }
    
    return(list(edges=edges,tree=tree,R=R,cdAlign=cdAlign,PI=PI))
  }
  
  #####Felsenstein pruning algorithm
  felsen=function(pos, tree, node, cdAlign, TPM) {
    edges = tree$edge
    prod = t(t(rep(1,61)))
    s = which(edges[,1] == node)
    for(k in s) {
      if(edges[k,2] <= length(tree$tip.label)){
        if(!is.null(dim(cdAlign))) {
          if(!is.na(cdAlign[tree$tip.label[edges[k,2]],pos])) {
            prod = prod * t(t(TPM[k,,cdAlign[tree$tip.label[edges[k,2]],pos]]))
          }
        }
        else {
          if(!is.na(cdAlign[tree$tip.label[edges[k,2]]])) {
            prod = prod * t(t(TPM[k,,cdAlign[tree$tip.label[edges[k,2]]]]))
          }
        }
      }
      else {
        prod = prod * TPM[k,,]%*%felsen(pos,tree,edges[k,2],cdAlign, TPM)
      }
    }
    return(prod)
  }
  
  
  ########End of setup########
  
  #Takes parameters: kappa, omega, treescale
  lik_fun = function(pars,treefile,seqfile, model=0) {
    kappa = pars[1]
    #omega = pars[2]
    beta_0 = pars[2]
    beta_1 = pars[3]
    treescale = pars[4]
    if(model==0) {
      beta_1 = 0
    } 
    
    if(kappa < 0 | treescale < 0) {
      return(NA)
    }
    out = setup(beta_0,beta_1,kappa,treefile,seqfile)
    edges = out$edges
    tree = out$tree
    cdAlign = out$cdAlign
    R = out$R
    PI = out$PI
    TPM = array(rep(0,61*61*nrow(edges)),dim=c(nrow(edges),61,61))
    for(i in 1:nrow(edges)){
      TPM[i,,] = expm(treescale*tree$edge.length[i]*R[i,,])
    }
    likelihood = 0
    for(i in 1:ncol(cdAlign)) {
      likelihood = likelihood + log(sum(diag(PI)*felsen(i, tree, edges[nrow(edges),1], cdAlign,TPM)))
    }
    print(c(pars,likelihood))
    return(likelihood)
  }
  
  
  
  ## The next lines run the code - you can modify them to run the code on whatever you want
  
  m0.out = optim(c(2,0.1,0,1),lik_fun,model=0,treefile=treefile,seqfile=seqfile,control=list(fnscale=-1))
  m1.out = optim(c(m0.out$par[1:2],0.01,m0.out$par[4]), lik_fun, model = 1, treefile = treefile, seqfile = seqfile, control = list(fnscale = -1))
  print(m0.out)
  print(m1.out)
  
  
  results <- rbind(results, data.frame(
    seqfile = basename(seqfile),
    m0_kappa = m0.out$par[1],
    m0_beta_0 = m0.out$par[2],
    m0_beta_1 = m0.out$par[3],
    m0_treescaling = m0.out$par[4],
    m0_ml = m0.out$value,
    m1_kappa = m1.out$par[1],
    m1_beta_0 = m1.out$par[2],
    m1_beta_1 = m1.out$par[3],
    m1_treescaling = m1.out$par[4],
    m1_ml = m1.out$value,
    stringsAsFactors = FALSE
  ))
}

write.table(results, file = "/Users/emilydunican/Desktop/results_growth.csv", sep = ",", row.names = FALSE)
