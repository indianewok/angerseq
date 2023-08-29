# Introduction

To whom it may concern:

Please find here the records of my trials, tribulations, and
stutter-steps in trying to figure out how to demultiplex this data.

I’ve attempted to embed it into this R Notebook so that whoever’s
looking at this can at least have an idea as to whatever semblance of
logic possessed me while I did this.

I’ve uploaded the commands as a package to be used, but they’re
kinda…wrapped in stuff so I figured I’d show how I used them down here.
Enough talk!

## Disclaimer

This is how to demultiplex **one** run.

This isn’t the smartest way to do this, and it certainly isn’t the best.

You need a **really big-boy computer** to run this code, specifically.

# Libraries

      libraries<-c("stringr","stringi", "seqRFLP","seqinr","phylotools","parallel","future",
                   "R.utils","data.table","pbapply","pbmcapply","stringdist","tidyverse",
                   "ggplot2","ggpubr", "tgutil", "progress", "DescTools","processx","trqwe",
                   "pafr","tools")
      invisible(lapply(libraries, library, character.only = TRUE))
      rm(libraries)

These are all the libraries you’ll need (maybe a few extras, I’m still
trimming dependencies). Additionally, you’ll need to install
[seqkit](https://anaconda.org/bioconda/seqkit) and
[canu](https://anaconda.org/bioconda/canu) for this pipeline to work. I
avoid using the package reticulate by just copy-pasting the path to the
seqkit and canu executables, which you can snag via the following steps:

1.  Download seqkit and canu to the same python environment (conda
    install -c bioconda seqkit & conda install -c canu).

2.  After it’s done downloading, just type “which seqkit” and “which
    canu”. The system should return the path of the working executable,
    which you can just add like this: seqkit\_path&lt;-“path/to/seqkit”.

# Processing & Intro Stats

      prepared_anger<-prepare2_anger(file_path, path_to_seqkit, path_to_output)

The first command that I do is run the prepare2\_anger command, which
does the following:

1.  Given a path to a directory, it first checks whether the file is one
    big one or a bunch of small ones.

2.  It’ll calculate the relevant QC/interesting factoids of the
    sequencing with seqkit stats, returning stuff like GC content,
    average length of read, number of reads, etc.

3.  If it’s a big file, which it will be if you get your data straight
    from the Nanopore sequencing, it’ll use seqkit split2 to split the
    file into a thousand smaller sub-files, which makes processing way
    easier.

# Importing data into R

    fastq_files<-list.files(path = path_to_output, full.names = TRUE, pattern = ".gz")
    batch_size<-ceiling(length(fastq_files)/10)
    anger_frame<-list()
    for(i in seq(1, length(fastq_files), by = batch_size)){
      j<-min(i+batch_size-1, length(fastq_files))
      anger_frame[[match(i, seq(1, length(fastq_files), by = batch_size))]]<-pblapply(
        fastq_files[i:j], function(x){
          read_fastq(fn = x)
      })
      anger_frame[[match(i, seq(1, length(fastq_files), by = batch_size))]]<-rbindlist(anger_frame[[match(i, seq(1, length(fastq_files), by = batch_size))]])
    }
    anger_frame<-pblapply(seq_along(anger_frame), function(y){
      colnames(anger_frame[[y]])<-c("id","fastq_files","qc")
      anger_frame[[y]][,c(1,3,2)]
    })

This is code that does two things: it imports *all* of the reads into R
in ten batches, makes them into a list of data.frames, and reformats
them to have three columns: id, fastq\_files, and qc.

# Doing the demux and concatenate error correction

    print("Now doing adapter distances and barcodes!")
    anger_frame<-anger_distance(anger_frame = anger_frame,
                                  primer_1 = primer_1,
                                  primer_2 = primer_2)
      anger_frame<-anger_pbarc(whitelist_path,anger_frame = anger_frame,
                               adapter_distance = 5,
                               primer_1 = primer_1)
      print("Now cleaning and doing error correction!")
      anger_frame<-lapply(anger_frame, function(x){
        stripping_in_anger(x) %>%
          erred_in_anger(.)
      })

This is the demultiplexing code, with a cute little progress bar for you
to see just how long this’ll take. Again, I had a really nice computer
to do this on, so runtimes will be very variable.

# Barcode Error Correction

      true_bcs<-list()
      lazy_bcs<-list()
      for(j in 1:length(anger_frame)){
        true_bcs[[j]]<-anger_frame[[j]]$true_bcs
        lazy_bcs[[j]]<-anger_frame[[j]]$lazy_bcs
      }
      true_bcs<-unlist(true_bcs) %>% table(.) %>% data.frame(.)
      lazy_bcs<-unlist(lazy_bcs) %>% table(.) %>% data.frame(.) %>% .[order(.$Freq, decreasing = TRUE),]
      lazy_bcs<-lazy_bcs[order(lazy_bcs$Freq, decreasing = TRUE),]
      lazy_bcs<-lazy_bcs[1:100000,]
      lazy_bcs$matched_true<-NA
      print("Now finding matching true barcodes for lazy barcodes that are 2 edits away, at most!")
      lazy_bcs$matched_true<-amatch(lazy_bcs$., table = true_bcs$., method = "lv", maxDist = 2)
      lazy_bcs<-lazy_bcs[which(!is.na(lazy_bcs$matched_true)),]
      lazy_bcs$matched_true<-true_bcs$.[lazy_bcs$matched_true]
      for(i in 1:length(anger_frame)){
        match_index<-match(anger_frame[[i]]$putative_bcs[which(anger_frame[[i]]$cr_barcode == FALSE)],
                           lazy_bcs$.) %>% .[!is.na(.)] %>% lazy_bcs$matched_true[.] %>% as.character(.)
        anger_frame[[i]][which(!is.na(match_index))]$putative_bcs<-match_index
        anger_frame[[i]][which(!is.na(match_index))]$cr_barcode<-TRUE
        print(paste0("Done matching and cleaning part number ",i,"!"))
      }

This is an attempt at barcode error correction: I take the top 10K
mismatched barcodes that show up the most frequently, and attempt to see
whether they are two or less edits away from any barcode in the “true”
list. It’s a cheap solve, as it doesn’t scan against the entire
whitelist every time, but it’s a good enough fix.

# QC Score Cleanup

      true_calls<-list()
      for(j in 1:length(anger_frame)){
        true_calls[[i]]<-anger_frame[[i]]$anger_frame
        true_calls[[i]]$qc_updated<-mapply(x = true_calls[[i]]$cleaned_seq, y = true_calls[[i]]$fastq_files, z = true_calls[[i]]$qc,
                                           function(x,y,z){
                                             str_locate(pattern = x, string = y) %>% str_sub(string = z, .)
                                           })
        true_calls[[i]]<-true_calls[[i]] %>% .[which(.$cr_barcode == TRUE),c("id","cleaned_seq","qc_updated","putative_bcs")]
      }

The sequencing data has already been stripped of the adapter sequences,
and anything left over; here, I trim the QC to match that.

# Saving individual fastq files

      true_calls<-rbindlist(true_calls) %>% split.data.frame(., f = .$putative_bcs)
      progress<-progress::progress_bar$new(total = length(true_calls),
                                           format = "  progress [:bar] :percent eta: :eta",)
      dir.create(path = paste0(path_to_output,"/fastq_output/"))
      save_folder<-paste0(path_to_output,"/fastq_output/")
      print(paste0("Writing the fastqs now...made a new folder at ",save_folder))
      for(k in 1:length(true_calls)){
        true_calls[[k]]<-true_calls[[k]][,c("id","cleaned_seq","qc_updated")]
        colnames(true_calls[[k]])<-c("id","seq","qual")
        write_fastq(df = true_calls[[k]], fn = paste0(save_folder,names(true_calls)[[k]],".fastq"))
        progress$tick()
      }

This chunk basically creates a new file directory where all of the new
.fastqs will be saved, creates the new fastq files, and saves them. The
format will be something like “AATCGTTTACTAGA.fastq”. From here on out,
I’m pasting the code that I used to process *all* of the data, not just
one file at a time, so you’ll have to edit it accordingly to just do one
if that’s all you’ve got.

# Filtering Cells by paired scRNA-seq metadata

      mcv_metadata<-read_csv("/path_to_metadata.csv")
      barcodes<-list()
      for(i in unique(mcv_metadata$dataset)){
        barcodes[[i]]<-mcv_metadata[grep(pattern = i, x = mcv_metadata$dataset),"...1"] %>%
          unlist(.) %>%
          gsub(pattern = "^LA._", replacement = "", x = .) %>%
          gsub(pattern = "-1$", replacement = "", x = .)
      }
      LA1_files<-list.files(path = save_folder_1, full.names = FALSE) %>% basename(.) %>% gsub(pattern = ".fastq", replacement = "", x = ., fixed = TRUE)
      LA2_files<-list.files(path = save_folder_2, full.names = FALSE) %>% basename(.) %>% gsub(pattern = ".fastq", replacement = "", x = ., fixed = TRUE)
      LA3_files<-list.files(path = save_folder_3, full.names = FALSE) %>% basename(.) %>% gsub(pattern = ".fastq", replacement = "", x = ., fixed = TRUE)
      barcode_list<-list(LA1_files, LA2_files, LA3_files)
      LA1_files<-list.files(path = save_folder_1, full.names = TRUE)
      LA2_files<-list.files(path = save_folder_2, full.names = TRUE)
      LA3_files<-list.files(path = save_folder_3, full.names = TRUE)
      filepath_list<-list(LA1_files, LA2_files, LA3_files)
      
      valid_files<-mapply(x = barcodes, y = barcode_list, z = filepath_list, function(x,y,z){
        match(x, y) %>% na.omit(.) %>% as.numeric(.) %>% z[.]
      })

Here’s a brief overview of what’s going on:

1.  Import the metadata from a Seurat object saved as a .csv.
2.  Look at a column in the metadata called “dataset” to see which
    barcodes come from which data sources.
3.  Store the barcodes, in a list called “barcodes”, matching whichever
    dataset they came from.
4.  Generate two lists: one is barcode list and one is filepath list.
    +barcode\_list is a list of *all* barcodes that the demux returned
    as being associated with cells in a certain dataset. +filepath\_list
    is a list of all of the paths to those files.
5.  Generate a list of lists called “valid\_files”, which basically
    overlaps the barcodes that passed all of your QC for your Seurat
    object with any that are also in the associated demux output. These
    will be the files that are used fo downstream processing.

There’s 100% smarter ways to do this, I just didn’t know any of them
when I did it.

# Canu Processing & Scratched-Out Minimap2 code

      for(i in 1:length(valid_files)){
        #save_path_minimap<-paste0(unique(dirname(valid_files[[i]])),"/minimapped/")
         save_path_canu<-paste0(unique(dirname(valid_files[[i]])),"/canu/")
         save_path_canu_commands<-paste0(save_path_canu,"canu_commands/")
        #dir.create(save_path_minimap)
         dir.create(save_path_canu)
         dir.create(save_path_canu_commands)
        #minimap2<-"/dartfs/rc/lab/A/AckermanLab/CMV/PostRot/cmv/miniconda3/envs/medaka/bin/minimap2"
        #pbsapply(valid_files[[i]], function(x){
        #minimap2_args<-c("/dartfs-hpc/rc/home/3/f005f43/trust4/TRUST4/mouse/GRCm38_bcrtcr.fa","-x","splice","-c",x,"-o",
        #                 paste0(save_path_minimap, file_path_sans_ext(basename(x)),".paf"))
        #processx::run(minimap2, minimap2_args, echo_cmd =FALSE, echo = FALSE, spinner = FALSE)
        #})
        batch_size<-ceiling(length(valid_files[[i]])/10)
          for(j in seq(1, length(valid_files[[i]]), by = batch_size)){
            k<-min(j+batch_size-1, length(valid_files[[i]]))
            print(paste0("Now writing canu commands for ", basename(dirname(unique(dirname(valid_files[[i]]))))," samples ",j," through ", k))
            canu<-path_to_canu
            cat(
            sapply(valid_files[[i]][j:k], function(x){
            paste(canu,"-d", paste0(save_path_canu, file_path_sans_ext(basename(x))),"-correct","stopOnLowCoverage=0.1",
                  "minInputCoverage=0.1","genomesize=15k","-minReadLength=100","-minOverlapLength=50","-correctedErrorRate=0.4","-p",
                  file_path_sans_ext(basename(x)),
                           "-nanopore-raw", x, sep = " ")
            }),
            file = paste0(save_path_canu_commands,"bulk_",j,"_",k,".txt"), sep = "\n")
        }
      }

This is basically two commands, though one is commented out. The most
important one, the canu one, is one I’ll describe:

1.  Create a new directory for Canu output
2.  Create a new directory to save the canu commands
3.  Define a batch size of 10, though it can be bigger or smaller
4.  For each valid barcode in “valid\_files”, create a canu command that
    can be run through the shell.
5.  Tag that on to a literal .txt document which has all of these
    commands.

*Listen*. I need to be real with you. Like, really real. Sometimes, when
all you have is a hammer (a really juvenile understanding of R, a very
nice computer, and a lack of patience), every problem (trying to go back
and forth between R and the command line, trying to run Canu on like ~8K
things at once) looks like a nail. So what do you do? Sometimes, instead
of figuring out how to write a for-loop in the command line–which you
can do now, it’s not a big deal–instead, you create a command that
literally copy-pastes the same words with different paths to .fastq
files, that you can then call from the command line using nohup sh. &,
and then just let it run until it’s done.

They don’t ask how, they ask how many. Or something.

# Canu Post-Processing

      files<-list.files(path="/dartfs/rc/lab/A/AckermanLab/CMV/PostRot/RAGE-Seq_RDS/TCR_outs/canu_outs", full.names = TRUE)
      files<-lapply(files, function(x){list.files(path = x, full.names = TRUE, recursive = FALSE)})
      files<-lapply(files, function(x){list.files(path = x, full.names = TRUE, recursive = FALSE)})
      assembled_files<-lapply(seq_along(files), function(x){files[[x]][which(grepl(pattern = ".contigs.fasta", x = files[[x]]))]})
      assembled_files<-unlist(assembled_files)
      test<-sapply(assembled_files, function(x){file.size(x)})
      files<-assembled_files[which(test > 0)]
      files<-rbindlist(files)

# Per-Dataset IGBLAST

Then, after doing canu error correction, I make a list of files that
come out of canu error correction, and use this command to make a large
fastq file that has everything error corrected. This generates basically
three extra-extra large files that I then run igblast on individually,
with a command in bash kinda like this:

    for i in $(ls path_to_TCR_outs/demux_corrected_fa/*.fa);do ~/igdata/ncbi-igblast-1.20.0/bin/igblastn -query $i -ig_seqtype TCR -germline_db_V mouse/imgt_mouse_tr_v -germline_db_D mouse/imgt_mouse_tr_d -germline_db_J mouse/imgt_mouse_tr_j -outfmt 19 -organism mouse -out path_to_TCR_outs/demux_corrected_fa/$(basename -s .fa $i).csv -auxiliary_data ~/igdata/ncbi-igblast-1.20.0/optional_file/mouse_gl.aux -clonotype_out path_to_TCR_outs/demux_corrected_fa/$(basename -s .fa $i)_clon.csv;echo $(basename -s .fa $i);done

# Final Step

Finally, after doing all of this, I can use the igblast\_importer to
import the large igblast files, and using the metadata from a Seurat
object, use the igmeta\_converter to convert the igBLAST data and append
it to the metadata at hand. Then, you can just reassign the metadata to
the Seurat object, and all of your data should be loaded in and good to
go!
