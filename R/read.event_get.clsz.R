#*********************************************
#*********************************************
#' Extracts clusters of segmentation data.
#'
#' @param out  is a list containing the segmentation data, typically out$sgsc when run from read.event().
#' @param beamsfiles  is a vector of the file names of the beams-files.
#' @param beamsfilesind  is a vector of the indexes of the beams-files in the list of files.
#' @param TIME  is the list returned from UNIX_time().
#' @param pamkpar  is a list of parameters ('krange', 'criterion', 'alpha' and 'mindist') used pamk() when clustering the segmented voxels 'sgsc'. The voxels 'sgsc' are also ordered so that the voxels belonging to the largest cluster lead and the second to largest cluster follows. A suggested set of parameters are pamkpar=list(krange=1:4,criterion="asw",alpha=0.05,N=1e2,mindist=100).
#' @param esnm  is the name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom fpc pamk
#' @importFrom stats dist
#'
#' @export
#' @rdname read.event_get.clsz
#'
read.event_get.clsz<-function(out, beamsfiles, beamsfilesind, tlist, pamkpar, esnm){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2013-05-12 - Clean version.
	

	##################################################
	##################################################
	if(!"sgsc" %in% names(out)){
		warning("Segmentation cluster size 'clsz' requested but not returned. Segmentation data 'sgsc' must be present for 'clsz' to be returned")
		}
	else{
		# Read beam configuration data required to extract the sonar volume dimensions and to regenerate uniformly distributed points in each segmented voxel (by the function read.event_sgsc2pos.TSD()), used in the clustering:
		#thisbeams=read.event_read_other_beams(c("lenb","freq","dira","dire","lenb","esnm","numb","eqba","sint","asps"), beamsfiles, beamsfilesind, onlyonestep=TRUE, TIME)
		thisbeams=read.event_read_files(var=c("lenb","freq","dira","dire","lenb","esnm","numb","eqba","sint","asps"), files=beamsfiles, filesind=beamsfilesind, tlist=tlist)
		
		J=max(thisbeams$lenb)
		I2=length(unique(thisbeams$freq))
		I1=length(thisbeams$freq)/I2
		# Set defaults of the paramterers used in the pamk(). The number of points in each segmented voxel, N, and the minimum required distance between cluster medoids are included in 'pamkpar', but removed when the function pamk() is run:
		if(!is.list(pamkpar)){
			pamkpar=as.list(pamkpar)
			}
		if(length(pamkpar$krange)==0){
			pamkpar$krange=1:4
			}
		if(length(pamkpar$criterion)==0){
			pamkpar$criterion="asw"
			}
		if(length(pamkpar$alpha)==0){
			pamkpar$alpha=0.01
			}
		if(length(pamkpar$N)==0){
			pamkpar$N=1e2
			}
		if(length(pamkpar$mindist)==0){
			pamkpar$mindist=100
			}
		N=pamkpar$N
		pamkpar$N=NULL
		mindist=pamkpar$mindist
		pamkpar$mindist=NULL
		# Assume large data set, and allways use clara():
		pamkpar$usepam=FALSE
		
		# Define the cluster size variable to return:
		out$clsz=vector("list",length(out$sgsc))
		for(i in seq_along(out$sgsc)){
			if(length(out$sgsc[[i]])>1){
				##################################################
				####### Old attempts, saved for reference: #######
				##################################################
				# Cluster the segmentation data by position in of the voxels:
				#thispsxx=data$psxx[sgsc]
				#thispsyx=data$psyx[sgsc]
				#thispszx=data$pszx[sgsc]
				#cl=pamk(ind2arr.ind(out$sgsc[[i]],c(J,I1,I2)),krange=intersect(krange,seq_len(length(out$sgsc[[i]])-1)),criterion=criterion,usepam=FALSE,alpha=alpha)$pam$clustering
				# Cluster the segmentation data by voxels index:
				#cl=pamk(ind2arr.ind(out$sgsc[[i]],c(J,I1,I2)),krange=intersect(krange,seq_len(length(out$sgsc[[i]])-1)),criterion=criterion,usepam=FALSE,alpha=alpha)$pam$clustering
				#N=min(floor(1e7/out$sgsc[[i]]),N)
				##################################################
				##################################################
				
				# Generate uniformly distributed points in each segmented cluster, using the simplified function read.event_sgsc2pos.TSD(), which assumes equal number of points in each voxel and uses the vessel coordinate system:
				sgscN=read.event_sgsc2pos.TSD(c(thisbeams,list(sgsc=out$sgsc[[i]])),N=N,esnm=esnm)
				
				# Apply clustering with automatic selection of the number of clusters:
				cl=do.call("pamk",c(list(data=sgscN),pamkpar))
				# Extract the distance between the medoids, and merge the clusters that are closer together than a minumum distance 'mindist'. If cluster 2 is too close to clustser 4 and cluster 4 is too close to cluster 5 but cluster 2 and 5 are suffucuently distant, then all three clusters are merged:
				distcl=as.matrix(dist(cl$pam$medoids))
				cl=cl$pam$clustering
				# Get the pairs of clusters to merge:
				mergecl=unique(t(apply(which(distcl<mindist & distcl>0,arr.ind=TRUE),1,sort)))
				if(length(mergecl)>0){
					newlabels=seq_len(nrow(distcl))
					# Merge each pair sucsessively:
					for(i_merge in seq_len(nrow(mergecl))){
						newlabels[mergecl[i_merge,2]]=newlabels[mergecl[i_merge,1]]
						}
					newlabels=match(newlabels,sort(unique(newlabels)))
					# Re-label the clusters, merging the close ones:
					cl=newlabels[cl]
					}
				# Remove row names of 'cl' and set the dimension to average the cluster allocation of each voxel:
				rownames(cl)=NULL
				dim(cl)=c(N,length(out$sgsc[[i]]))
				cl=round(colMeans(cl))
				# Shift the cluster indices to that the largest cluster has label 1, the second to largest cluster has label 2, and so on:
				out$clsz[[i]]=table(cl)
				ranksize=length(out$clsz[[i]])+1-rank(out$clsz[[i]])
				cl=ranksize[cl]
				# Add sizes of the clusters to 'data':
				out$clsz[[i]]=sort(out$clsz[[i]],decreasing=TRUE)
				# Reorganize 'sgsc' so that the voxels in the largest cluster comes first, the voxels in the second to largest cluster follows, and so on:
				out$sgsc[[i]]=out$sgsc[[i]][order(cl)]
				}
			}
		}
	out
	##################################################
	##################################################
	}
