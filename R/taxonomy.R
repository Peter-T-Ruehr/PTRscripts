#' Get taxize Gbif taxonomy data for a list of taxonomic genera.
#' @param genera A vector of `character` values for taxonomic genera.
#' @param taxon A `character` string for taxonomic context to search. Can be
#' `Animalia`, `Arthropoda`, or `Insecta`.
#' @return A data frame in the `tibble` format with the colums `genus`, `class`,
#' `order`, and `family`.
#' @export
#' @examples
#' # add leading zeros to numbers
#'
#' taxonomy_gbif(c("Caenorhabditis", "Drosophila", "Xenopus"), taxon = "Animalia")
#'
#' taxonomy_gbif(c("Argiope", "Limulus", "Folsomia"), taxon = "Arthropoda")
#'
#' taxonomy_gbif(c("Calopteryx", "Forficula", "Drosophila"), taxon = "Insecta")
#'
taxonomy_gbif <- function(genera,
                          taxon = "Insecta"){
  require(tibble)
  require(taxize)
  require(rgbif)
  require(dplyr)

  # get taxonomy data from ncbi: slow, but many taxonomic levels
  # curr.taxonomy.all.genera <- classification(genera, db = 'gbif', output = "classification")
  if(is.null(genera)){
    genera <- readClipboard()[-1]
  }
  full.taxonomy <- empty_tibble(names = c("genus"), nrow = 0, type = "factor") # length(curr.taxonomy.all.genera)
  # full.taxonomy$genus <- names(curr.taxonomy.all.genera)

  if(taxon == "Animalia"){
    taxon_key <- 1
  } else if(taxon == "Arthropoda"){
    taxon_key <- 54
  } else {
    taxon_key <- get_gbifid(taxon)
  }

  for(g in 1:length(genera)){ # length(genera)
    print(paste0(g, ": ", genera[g]))
    # for(g in 1:length(curr.taxonomy.all.genera)){
    # print(g)
    curr.genus <- genera[g]
    # curr.genus <- names(curr.taxonomy.all.genera)[g]
    if(g==1){last.genus <- ""}

    if(curr.genus == ""){curr.genus <- "NA"}

    if(curr.genus != last.genus){

      classification(taxon_key, db = 'gbif', output = "classification", rows = 1)[[1]]

      curr.lookup <- name_lookup(curr.genus,
                                 rank = "genus",
                                 status = "ACCEPTED",
                                 higherTaxonKey = taxon_key)

      # get gbif ID
      # curr.gbif.ID <- get_gbifid(curr.genus,
      #                            class = "Insecta",
      #                            rank = "genus")


      if(length(curr.lookup[2]$data > 0)){
        curr.gbif.ID <- curr.lookup[2]$data$key
        # get taxize data for curr. genus
        curr.taxonomy_ <- classification(curr.gbif.ID, db = 'gbif', output = "classification", rows = 1)[[1]]
        # curr.taxonomy_ <- curr.taxonomy.all.genera[[g]]


        # check if "class" is part of rank (not the case for basal BaHe)
        if("class" %in% curr.taxonomy_$rank){
          # get class line
          class.line <- which(curr.taxonomy_$rank == "class")

          # print(class.line)
          # only save below & incl. subclass and delete clade lines
          curr.taxonomy <- curr.taxonomy_[(class.line:nrow(curr.taxonomy_)),2:1] %>%
            filter(rank != "clade" & rank != "no rank")

          # convert to transposed data frame
          curr.taxonomy.t <- as.data.frame(t(curr.taxonomy))

          # make rank names colnames
          colnames(curr.taxonomy.t) <- as.character(unlist(curr.taxonomy.t[1,]))

          # delete rank line
          curr.taxonomy.finished <- as_tibble(curr.taxonomy.t[-1,])
        } else {
          curr.taxonomy.finished <- empty_tibble(names = "genus", nrow = 1)
          curr.taxonomy.finished$genus <- curr.genus
        }

      }else {
        curr.taxonomy.finished <- empty_tibble(names = "genus", nrow = 1)
        curr.taxonomy.finished$genus <- curr.genus
      }
    } else {
      # do nothing <- use last curr.taxonomy.finished
    }
    # bind rows of cumulative taxonomy and curr. taxonomy
    full.taxonomy <- bind_rows(full.taxonomy, curr.taxonomy.finished)

    last.genus <- curr.genus

    print_progress(g, length(genera))
    # print_progress(g, length(curr.taxonomy.all.genera))
  }


  # delete class column (they're all Insecta Pterygota)
  full.taxonomy.save <- full.taxonomy %>%
    select(-class)

  # write.xlsx2(as.data.frame(full.taxonomy.save), paste0("N:/PAPERS/PTR_Bite force across clades and scales/phylogeny/", today(), "_Chesters_full_taxonomy_gbif.xlsx"), row.names = F)


  return(full.taxonomy)
}
