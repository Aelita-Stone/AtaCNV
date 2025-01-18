#' Generate input matrix from fragment file for AtaCNV pipeline
#'
#' @param fragment_file path of fragment file.
#' @param cell_barcodes a string vector giving barcodes that represent cells in
#' scATAC-seq expriment.
#' @param genome reference genome, "hg19", "hg38" or "mm10".
#' @param output_dir output directory
#'
#' @return a cell-by-bin read count matrix
#' @export
#'
#' @examples
#'
generate_input <- function(fragment_file,
                           cell_barcodes,
                           genome="hg38",
                           output_dir="./"){
  if(genome=="hg19"){
    data("bin_info_hg19")
    bin <- bin_info_hg19
  }else if(genome=="hg38"){
    data("bin_info_hg38")
    bin <- bin_info_hg38
  }else if(genome=="mm10"){
    data("bin_info_mm10")
    bin <- bin_info_mm10
  }else{
    print("Genome should be one of {hg19, hg38, mm10}.")
  }
  print(paste0("Reference genome: ", genome))
  chr_list <- unique(bin$chr)
  chr_length <- rep(0, length(chr_list))
  for(i in 1:length(chr_list)){
    chr_length[i] <- max(bin$end[bin$chr==chr_list[i]])
  }
  barcode <- cell_barcodes
  K <- 1000000

  count <- matrix(data=0, nrow=length(barcode), ncol=nrow(bin))
  rownames(count) <- barcode
  colnames(count) <- bin2str(bin)

  t1 <- proc.time()
  # con <- file(fragment_file, "r")
  con <- gzcon(file(fragment_file, open="r"))
  flag <- 0
  nchr <- 1
  nchr_start <- 0
  while(TRUE){
    flag <- flag+1
    line <- readLines(con, n=1)
    # print(line)
    if(line %% 1000 == 0){
      print(line)
    }
    if(length(line) == 0){
      break
      print("all finished")
    }
    temp <- strsplit(line, "\t")[[1]]
    cell_barcode <- temp[4]
    if(!(cell_barcode %in% barcode) | !(temp[1]) %in% chr_list){
      next
    }
    if(temp[1] != chr_list[nchr]){
      print(paste(chr_list[nchr], "finished"))
      nchr_start <- nchr_start + floor(chr_length[nchr]/K)+1
      nchr <- nchr + 1
      if(nchr > length(chr_list)){
        break
      }
    }
    pos <- as.numeric(temp[2])
    count[cell_barcode, nchr_start+floor(pos/K)+1] <-
      count[cell_barcode, nchr_start+floor(pos/K)+1] + as.numeric(temp[5])
  }
  close(con)
  t2 <- proc.time()
  print(t2-t1)

  print(paste0("Saving to ", output_dir))
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  saveRDS(count,file=paste0(output_dir, "count.rds"))
}
