#*********************************************
#*********************************************
#' Returning a categorical array (0 and 1) for the condition defined by 'a', 'b' and 'cond'.
#'
#' @param x  is the input array.
#' @param a  is the lower border given as a single value or as a vector c(a,b).
#' @param b  is the upper border. If 'b' is a character, it is interpreted as 'cond'..
#' @param cond  specifies the border condition type:
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname contains
#'
contains<-function(x,a,b=NULL,cond="less"){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2008-10-08 - Finished.
	# Update: 2009-02-23 - Added support for vectors of border points.
	# Update: 2009-02-25 - Expanded to support all condition types {"x<a", "x<=a", x==a", "x>a", "x>=a", "a<x<b", "a<=x<b", "a<x<=b", "a<=x<=b", "x<a,x>b", "x<=a,x>b", "x<a,x>=b" and "x<=a,x>=b"}
	# Update: 2009-05-06 - Cleaned up by introducing the condition variable 'condnum' and transforming to this at the preparation-part. Added support for entering 'cond' as the 3rd input. For onesided tests (cond in 1:5) 'b' is simply ignored, and for thosided tests (cond in 6:13) b=a[-1] and a=a[-length(a)].
	# Last: 2009-05-18 - Support for more than one condition removed, as it was regarded unnecessary and greatly complicated the function. Old version saved in the "-unused"-directory.
	########### DESCRIPTION: ###########
	# Returning a categorical array (0 and 1) for the condition defined by 'a', 'b' and 'cond'.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---x--- is the input array.
	# ---a--- is the lower border given as a single value or as a vector c(a,b).
	# ---b--- is the upper border. If 'b' is a character, it is interpreted as 'cond'..
	# ---cond--- specifies the border condition type:
	#		1 - "l" - "a>x" - "x<a" - "less"
	#		2 - "le" - "a>=x" - "x<=a" - "lessequal"
	#		3 - "e" - "a=x" - "x=a" - "a==x" - "x==a" - "equal"
	#		4 - "g" - "a<x" - "x>a" - "greater"
	#		5 - "ge" - "a<=x" - "x>=a" - "greaterequal"
	#		6 - "g,l" - "a<x<b" - "b>x>a" - "greater,less" - "b" - "between"
	#		7 - "ge,l" - "a<=x<b" - "b>x>=a" - "greaterequal,less"
	#		8 - "g,le" - "a<x<=b" - "b>=x>a" - "greater,lessequal"
	#		9 - "ge,le" - "a<=x<=b" - "b>=x>=a" - "greaterequal,lessequal" - "o" - "onbetween"
	#		10 - "l,g" - "a>x,b<x" - "a>x,x>b" - "x<a,b<x" - "x<a,x>b" - "less,greater"
	#		11 - "le,g" - "a>=x,b<x" - "a>=x,x>b" - "x<=a,b<x" - "x<=a,x>b" - "lessequal,greater"
	#		12 - "l,ge" - "a>x,b<=x" - "a>x,x>=b" - "x<a,b<=x" - "x<a,x>=b" - "less,greaterequal"
	#		13 - "le,ge" - "a>=x,b<=x" - "a>=x,x>=b" - "x<=a,b<=x" - "x<=a,x>=b" - "lessequal,greaterequal"
	
	
	##################################################
	##################################################
	##### Preparation #####
	# Dimension of the input 'x':
	dimx=dim(x)
	# If 'b' is missing or is a character, onesided matching is done:
	oneborder=FALSE
	if(is.character(b)){
		oneborder=TRUE
		cond=b
		}
	else if(is.null(b)){
		oneborder=TRUE
		}
		
	# Stripping 'cond' of whitespaces:
	cond=gsub(" ","",cond)
	# Transforming 'cond' into an integer in the range [1,13]:
	condnum=numeric(0)
	if(cond[1]=="l" || cond[1]=="a>x" || cond[1]=="x<a" || cond[1]=="less"){
		condnum=1
		}
	else if(cond[1]=="le" || cond[1]=="a>=x" || cond[1]=="x<=a" || cond[1]=="lessequal"){
		condnum=2
		}
	else if(cond[1]=="e" || cond[1]=="a=x" || cond[1]=="x=a" || cond[1]=="a==x" || cond[1]=="x==a" || cond[1]=="equal"){
		condnum=3
		}
	else if(cond[1]=="g" || cond[1]=="a<x" || cond[1]=="x>a" || cond[1]=="greater"){
		condnum=4
		}
	else if(cond[1]=="ge" || cond[1]=="a<=x" || cond[1]=="x>=a" || cond[1]=="greaterequal"){
		condnum=5
		}
	else if(cond[1]=="g,l" || cond[1]=="a<x<b" || cond[1]=="b>x>a" || cond[1]=="greater,less" || cond[1]=="b" || cond[1]=="between"){
		condnum=6
		}
	else if(cond[1]=="ge,l" || cond[1]=="a<=x<b" || cond[1]=="b>x>=a" || cond[1]=="greaterequal,less"){
		condnum=7
		}
	else if(cond[1]=="g,le" || cond[1]=="a<x<=b" || cond[1]=="b>=x>a" || cond[1]=="greater,lessequal"){
		condnum=8
		}
	else if(cond[1]=="ge,le" || cond[1]=="a<=x<=b" || cond[1]=="b>=x>=a" || cond[1]=="greaterequal,lessequal" || cond[1]=="o" || cond[1]=="onbetween"){
		condnum=9
		}
	else if(cond[1]=="l,g" || cond[1]=="a>x,b<x" || cond[1]=="a>x,x>b" || cond[1]=="x<a,b<x" || cond[1]=="x<a,x>b" || cond[1]=="less,greater"){
		condnum=10
		}
	else if(cond[1]=="le,g" || cond[1]=="a>=x,b<x" || cond[1]=="a>=x,x>b" || cond[1]=="x<=a,b<x" || cond[1]=="x<=a,x>b" || cond[1]=="lessequal,greater"){
		condnum=11
		}
	else if(cond[1]=="l,ge" || cond[1]=="a>x,b<=x" || cond[1]=="a>x,x>=b" || cond[1]=="x<a,b<=x" || cond[1]=="x<a,x>=b" || cond[1]=="less,greaterequal"){
		condnum=12
		}
	else if(cond[1]=="le,ge" || cond[1]=="a>=x,b<=x" || cond[1]=="a>=x,x>=b" || cond[1]=="x<=a,b<=x" || cond[1]=="x<=a,x>=b" || cond[1]=="lessequal,greaterequal"){
		condnum=13
		}
	else if(is.numeric(cond)){
		if(0<cond[1] && 14>cond[1]){
			condnum=cond[1]
			}
		else{
			stop("'cond' does not match the predefined values")
			condnum=1
			}
		}
	else{
		stop("'cond' does not match the predefined values")
		condnum=1
		}
		
	# If tests are to be done two-sided, 'b' need to be given: 
	if(condnum %in% 6:13 && oneborder){
		if(length(a)<2){
			stop("Two border conditions must be given (see info(\"contains\"))")
			}
		else{
			b=a[2]
			a=a[1]
			}
		}
	# Ensuring clean performance:
	#a=a[1]
	#b=b[1]
	 
	
	##### Execution and output #####
	# If 'a' and 'b' have length==1, the simplest version is executed:
	if(condnum==1){
		out=a>x
		dim(out)=dimx
		return(out)
		}
	if(condnum==2){
		out=a>=x
		dim(out)=dimx
		return(out)
		}
	if(condnum==3){
		out=a==x
		dim(out)=dimx
		return(out)
		}
	if(condnum==4){
		out=a<x
		dim(out)=dimx
		return(out)
		}
	if(condnum==5){
		out=a<=x
		dim(out)=dimx
		return(out)
		}
	if(condnum==6){
		out=a<x & x<b
		dim(out)=dimx
		return(out)
		}
	if(condnum==7){
		out=a<=x & x<b
		dim(out)=dimx
		return(out)
		}
	if(condnum==8){
		out=a<x & x<=b
		dim(out)=dimx
		return(out)
		}
	if(condnum==9){
		out=a<=x & x<=b
		dim(out)=dimx
		return(out)
		}
	if(condnum==10){
		out=a>x | x>b
		dim(out)=dimx
		return(out)
		}
	if(condnum==11){
		out=a>=x | x>b
		dim(out)=dimx
		return(out)
		}
	if(condnum==12){
		out=a>x | x>=b
		dim(out)=dimx
		return(out)
		}
	if(condnum==13){
		out=a>=x | x>=b
		dim(out)=dimx
		return(out)
		}
	##################################################
	##################################################
	}
