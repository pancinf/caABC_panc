#Generates an input bed file for ABC
###
###
###Libraries
library(optparse)

###
###
###Arguments
option_list = list(
  make_option("--act", type="character", default=NULL,
              help="Path to input act file"),

  make_option("--actOut", type="character", default=NULL,
              help="Path to input list"),

  make_option("--actCol", type="numeric", default=NULL,
              help="Path to output act file")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$act) | is.null(opt$actOut) | is.null(opt$actCol)) {
  print("You have to specify all inputs. Here is the help:")
  print_help(opt_parser)
  quit()
}
###
###
###
EnhancerList <- read.delim(opt$act, header=TRUE)
EnhancerList <- EnhancerList[,c(1,2,3,opt$actCol)]
write.table(x = EnhancerList,file = opt$actOut,sep = "\t",append = FALSE,quote = FALSE,row.names = FALSE,col.names = FALSE)
