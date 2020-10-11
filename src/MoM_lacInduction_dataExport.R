# ZENODO ##########
(function(export_dir = "/scicore/home/nimwegen/GROUP/MM_Data/Thomas/_lacInductionArticle/Julou_2020_lacInduction_RawImages",
          .force=FALSE) {
# this anonymous function exports raw data (one dataset per compressed file) after the following scheme:
# - traverse R data to list all datasets used and corresponding files
# - create a directory per experiments and symlinks to corresponding files
# - create a compressed archive per dataset (following symlinks)
  
  # browser()
  dir.create(export_dir, showWarnings=FALSE)
  
  # delete existing symlinks and directories to avoid creating undesired archive files
  if (length( fs::dir_ls(export_dir, type = c("directory", "symlink")) )) {
    if (!.force)
      stop('existing directories and/or symlinks in the export directory: please delete them or run with .force=FALSE (note that this option will delete them but force the compression of all files).')
  } else {
    system(sprintf('find %s -maxdepth 1 -type l -exec rm {} +', export_dir))
    system(sprintf('find %s -maxdepth 1 -mindepth 1 -type d -exec rm -rf {} +', export_dir)) 
  }
  
  
  myconditions %>% 
    map(~.$paths) %>% 
    unlist() %>% 
    tibble(cpath=.) %>% 
    filter(str_detect(cpath, "data_thomas")) %>% 
    extract(cpath, 'date', 'data_thomas/(\\d{8})/', remove = FALSE) %>% 
    mutate(path = map(cpath, ~{#browser();
      fs::dir_ls(fs::path(., '..')) %>% as.character()
    })) %>% 
    unnest(path) %>% 
    mutate(path=fs::path_norm(path)) %>% 
    
    # filter(str_detect(path, 'Pos\\d+$')) %>% 
    # mutate(dsize = map(path, ~fs::dir_map(., fs::file_size, recurse = TRUE) %>% unlist() %>% sum() %>% fs::as_fs_bytes() )) %>% 
    # unnest(dsize) %>% 
    # pull(dsize) %>% sum() %>% fs::as_fs_bytes()
    
    filter(!str_detect(path, 'Pos\\d'), !str_detect(path, 'curated')) %>%
    filter(!str_detect(path, '20151127_flatfield'),
           !str_detect(path, '20190605_ASC662'),
           !str_detect(path, '20170108_gluLac_lac_snapAfter'),
           !path %in% c('data_thomas/20180514/20180515_gfpFlatfield_pos',
                        'data_thomas/20180606/20180606_gluLac_lac_switch16h_1')) %>% 
    mutate(filename = fs::path_file(path), 
           epath = fs::path(export_dir, date, filename),
           is_dir = fs::is_dir(path)) %>% 
    # print(n=Inf)
    # copy files and dirs
    # with( pwalk(list(path, epath), ~{ #browser();
    #   fs::dir_create(fs::path_dir(..2), recurse = TRUE)
    #   if (fs::is_file(..1)) fs::file_copy(..1, ..2)
    #   if (fs::is_dir(..1)) fs::dir_copy(..1, ..2)
    # }) ) %>% 
    # create symlinks
    with( pwalk(list(path, epath), ~{ #browser();
      fs::dir_create(fs::path_dir(..2), recurse = TRUE)
      fs::link_create(fs::path_real(..1), ..2)
    }) ) %>%
    identity()
  
  # create files list before compressing
  # use bash command as fs doesnt allow to follow symlinks when traversing directories
  # # deprecated
  # system2("tree", c("-l", export_dir), stdout = fs::path(export_dir, paste0(fs::path_file(export_dir), '_fileList.txt')))
  # system2("find", c(export_dir, "-maxdepth 1 -mindepth 1 -type d -exec tree {} +"), 
  #         stdout = fs::path(export_dir, paste0(fs::path_file(export_dir), '_fileList.txt')))
  # same but alphabetically sorted 
  system2("find", c("-L", export_dir, "-mindepth 1 -type d -print0 | sort -z | xargs -r0 tree"), 
          stdout = fs::path(export_dir, paste0(fs::path_file(export_dir), '_fileList.txt')))
  message("Done listing files...")
  
  # compress experiment by experiment (and resolve symlinks)
  fs::dir_walk(export_dir, function(.p) {# browser();
    if (fs::is_dir(.p) | (fs::is_link(.p))) {
      # tar(paste0(.p, '.tar.gz'), .p, compression = 'gzip', extra_flags = '--dereference')
      # fs::dir_delete(.p)
      # message("Done compressing ", .p)
      
      .f <- fs::path_file(.p)
      if (fs::file_exists(paste0(.p, '.tar.gz')) && !.force) {
        message ('skipping the archiving of ', .f)
        fs::dir_delete(.p)
      } else {
#       sprintf("sbatch <<EOF
# #!/bin/bash
# #SBATCH --job-name=tar_raw%s
# #SBATCH --mem=32G
# #SBATCH --time=1-0:00:00
# #SBATCH --qos=1day
# #SBATCH -o slogs/$JOB_NAME.o$JOB_ID
# #SBATCH -e slogs/$JOB_NAME.e$JOB_ID
# cd %s
# tar --dereference -czf %s.tar.gz %s
# rm -rf %s
# EOF", .f, export_dir, .f, .f, .f) %>%
#         system()
      # message("Job submitted for archiving ", .p)
    } }
  })
  
})()


(function(export_dir = "/scicore/home/nimwegen/GROUP/MM_Data/Thomas/_lacInductionArticle/Julou_2020_lacInduction_GL_Images",
          .force=FALSE) {
  # this anonymous function exports preprocessed image data (one dataset per compressed file) after the following scheme:
  # - traverse R data to list all datasets used
  # - create symlink to each dataset (renamed in a systematic manner)
  # - create a compressed archive per dataset (following symlinks)
  
  # browser()
  # delete existing symlinks to avoid creating undesired archive files
  if (length( fs::dir_ls(export_dir, type = "symlink") )) {
    if (!.force)
      stop('existing symlinks in the export directory: please delete them or run with .force=FALSE (note that this option will delete them but force the compression of all files).')
  } else {
    system(sprintf('find %s -maxdepth 1 -type l -exec rm {} +', export_dir))
  }

  fs::dir_create(export_dir, recurse = TRUE)
  myconditions %>% 
    map(~.$paths) %>% 
    unlist() %>% 
    tibble(path=.) %>% 
    filter(str_detect(path, "data_thomas")) %>% 

    # copy files and dirs
    mutate(filename=fs::path_file(path)) %>% 
    extract(filename, 'date', '^(\\d{8})_', remove = TRUE) %>% 
    mutate(epath = fs::path(export_dir, date)) %>% 
    # print(n=Inf)
    with( pwalk(list(path, epath), ~{ # browser();
      fs::dir_create(fs::path_dir(..2), recurse = TRUE)
      fs::link_create(fs::path_real(..1), ..2)
    }) ) %>%
    identity()
  
  # create files list before compressing
  # use bash command as fs doesnt allow to follow symlinks when traversing directories
  # # deprecated
  # system2("tree", c("-l", export_dir), stdout = fs::path(export_dir, paste0(fs::path_file(export_dir), '_fileList.txt')))
  system2("find", c(export_dir, "-maxdepth 1 -type l -print0 | sort -z | xargs -r0 tree"), 
          stdout = fs::path(export_dir, paste0(fs::path_file(export_dir), '_fileList.txt')))
  message("Done listing files...")
  
  # compress experiment by experiment (and resolve symlinks)
  fs::dir_walk(export_dir, function(.p) {# browser();
    if (fs::is_dir(.p) | (fs::is_link(.p))) {
      # tar(paste0(.p, '.tar.gz'), .p, compression = 'gzip', extra_flags = '--dereference')
      # fs::dir_delete(.p)
      # message("Done compressing ", .p)
      .f <- fs::path_file(.p)
      if (fs::file_exists(paste0(.p, '.tar.gz')) && !.force) {
        message ('skipping the archiving of ', .f)
        fs::link_delete(.p)
      } else {
        sprintf("sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=tar_GL%s
#SBATCH --mem=32G
#SBATCH --time=1-0:00:00
#SBATCH --qos=1day
#SBATCH -o slogs/$JOB_NAME.o$JOB_ID
#SBATCH -e slogs/$JOB_NAME.e$JOB_ID
cd %s
tar --dereference -czf %s.tar.gz %s
rm %s
EOF", .f, export_dir, .f, .f, .f) %>%
        system()
      # message("Job submitted for archiving ", .p)
    } }
  })

})()


(function(export_dir = "/scicore/home/nimwegen/GROUP/MM_Data/Thomas/_lacInductionArticle/Julou_2020_lacInduction_GL_Preproc") {
  # this anonymous function exports parsed output (tabular data) after the following scheme:
  # - traverse R data to list all datasets used
  # - create symlink to each dataset (renamed in a systematic manner)
  # - create a compressed archive (following symlinks)
  
  # browser()
  fs::dir_create(export_dir, recurse = TRUE)

  # delete existing symlinks to avoid creating undesired archive files
  if (length( fs::dir_ls(export_dir, type = "symlink") )) {
    if (!.force)
      stop('existing symlinks in the export directory: please delete them or run with .force=FALSE (note that this option will delete them but force the compression of all files).')
  } else {
    system(sprintf('find %s -maxdepth 1 -type l -exec rm {} +', export_dir))
  }
  
  # start from myconditions rather than listing ./preproc in order to export only used data
  myconditions %>%
    map(~.$paths) %>% 
    unlist() %>% 
    tibble(cpath=.) %>% 
    filter(str_detect(cpath, "data_thomas")) %>% 
    extract(cpath, 'date', 'data_thomas/(\\d{8})/', remove = FALSE) %>% 
    mutate(path = fs::path("preproc", date),
           epath = fs::path(export_dir, date) ) %>% 
    # print(n=Inf)
    # create symlinks
    with( pwalk(list(path, epath), ~fs::link_create(fs::path_real(..1), ..2) ) ) %>%
    identity()
  
  # create files list before compressing
  system2("tree", c("-l", export_dir), stdout = fs::path(export_dir, '..', paste0(fs::path_file(export_dir), '_fileList.txt')))
  
  # compress all together (and resolve symlinks)
  sprintf("cd %s; tar --dereference -czf %s.tar.gz %s", 
          fs::path_dir(export_dir), fs::path_file(export_dir), fs::path_file(export_dir)) %>% 
    system()

  fs::dir_delete(export_dir)
})()


# IDR metadata #######
(function() {
  # this anonymous function exports IDR metadata after the following scheme:
  # - traverse R data to list all datasets used and corresponding files
  # - create a directory per experiments and symlinks to corresponding files
  # - create a compressed archive per dataset (following symlinks)
  
  bind_rows(
    myconditions %>% 
      map(~{tibble(condition=.$condition, cpath=.$paths)}) %>% 
      bind_rows() %>% 
      filter(str_detect(cpath, "data_thomas")) %>% 
      extract(cpath, 'date', 'data_thomas/(\\d{8})/', remove = FALSE) %>% 
      mutate(path = map(cpath, ~{#browser();
        fs::dir_ls(fs::path(., '..')) %>% as.character()
      })) %>% 
      unnest(path) %>% 
      mutate(path=fs::path_norm(path)) %>% 
      
      filter(!str_detect(path, 'Pos\\d'), !str_detect(path, 'curated')) %>%
      filter(!str_detect(path, '20151127_flatfield'),
             !str_detect(path, '20190605_ASC662'),
             !str_detect(path, '20170108_gluLac_lac_snapAfter'),
             !path %in% c('data_thomas/20180514/20180515_gfpFlatfield_pos',
                          'data_thomas/20180606/20180606_gluLac_lac_switch16h_1')) %>% 
      filter(fs::is_dir(path)) %>% 
      mutate(path=map_chr(path, ~fs::dir_ls(., glob="*.ome.tif")[1] ),
             file=fs::path_file(as.character(path))
      ) %>% 
      # with(interaction(date, file)) %>% unique %>% length
      mutate(type=ifelse(str_detect(path, fixed('flatfield', ignore_case=TRUE)), 'flatfield', 'raw')) %>% 
      select(condition, dataset=date, file, path, type) %>% 
      mutate(path=str_replace(path, 'data_thomas', 'Julou_2020_lacInduction_RawImages'),
             preproc='') %>% 
      identity(),
    
    myconditions %>% 
      map(~{tibble(condition=.$condition, cpath=.$paths)}) %>% 
      bind_rows() %>% 
      filter(str_detect(cpath, "data_thomas")) %>% 
      mutate(filename=fs::path_file(cpath)) %>% 
      extract(filename, 'date', '^(\\d{8})_', remove = TRUE) %>%
      
      # filter(date %in% c('20151204', '20180313')) %>%
      
      mutate(path=map(cpath, ~fs::dir_ls(., type = "directory", recurse=T) %>% as.character()) ) %>% 
      unnest(path) %>% 
      mutate(path=map_chr(path, ~fs::dir_ls(., glob="*.tif")[1]),
             file=fs::path_file(as.character(path)) ) %>% 
      mutate(type='derived') %>% 
      extract(path, c("pos", "gl"), ".*\\d{8}_.*[Pp]os(\\d+).*_GL(\\d+).*", remove=FALSE) %>%
      left_join(
        myfiles %>% 
          mutate(preproc=data2preproc(path)) %>% 
          extract(path, c("date", "pos", "gl"), ".*(\\d{8})_.*[Pp]os(\\d+).*_GL(\\d+).*", remove=FALSE) %>%
          ungroup() %>% select(date, pos, gl, preproc) %>% 
          identity(),
        by=c('date', 'pos', 'gl')
      ) %>% 
      select(condition, dataset=date, file, path, type, preproc) %>% 
      mutate(path=str_replace(path, '\\./data_thomas', 'Julou_2020_lacInduction_GL_Images')) %>% 
      mutate(preproc=str_replace(preproc, '\\./preproc', 'GL_Preproc')) %>% 
      identity(),
  ) %>% 
    
    left_join(readxl::read_xlsx('data_thomas/MM_Experiments_Table.xlsx') %>% 
                transmute(dataset=date, strain=str_replace_all(strain, '\r?\n', ' | '), 
                          media = str_replace_all(media, '\r?\n', ' | '), 
                          steps = str_replace_all(steps, '\r?\n', ' | '), interval_min, 
                          flow_control = str_replace_all(flow_control, '\r?\n', ' | ')),
              by='dataset') %>% 
              
    transmute(
      'Source Name' = 'Julou_2020-lacInduction',  # an identifier for the sample
      'Characteristics [Organism]' = 'Escherichia coli',  # the species of the sample e.g.Homo sapiens
      'Term Source 1 REF' = 'NCBITaxon',
      'Term Source 1 Accession' = '', # leave blank
      'Characteristics [Cell Line]' = 'ASC662 (MG1655 lacZ-GFPmut2)', # the name of the cell line used e.g HeLa (if appropriate)  Other Characteristics columns can be added as necessary e.g. Characteristics [Organism Part]
      'Term Source 2 REF' = 'EFO',
      'Term Source 2 Accession' = '', # leave blank
      'Protocol REF' = 'growth protocol',
      'Protocol REF ' = 'treatment protocol', 
      'Protocol REF  ' = 'image acquisition and feature extraction protocol',
      'Protocol REF   ' = 'image analysis protocol', # data analysis protocol
      'Assay Name' = condition,  # an name for the imaging assay. This column can be used to group several images together If several images (e.g. fields, or replicate images) have been taken of the same sample then repeat the assay name on each row corresponding to the imaging file.  It can also be used to group raw and processed images from the same assay together
      'Experimental Condition [Media]' = media, # add an experimental condition if there is one e.g. 'targeted protein', or 'antibody'
      'Experimental Condition [Steps]' = steps, # add an experimental condition if there is one e.g. 'targeted protein', or 'antibody'
      'Experimental Condition [Flow Control]' = flow_control, # add an experimental condition if there is one e.g. 'targeted protein', or 'antibody'
      'Comment [Preculture]' = strain, 
      # 'Comment [Acquisition Interval]' = interval_min, 
      'Comment [Gene Identifier]' = '', # enter an identifier for any associated genes e.g. ensembl identifier for targeted protein
      'Comment [Gene Symbol]' = '', # enter gene symbol for any associated genes e.g. BRCA1
      'Comment [Gene Annotation Comments]' = '', # any comments about the gene annotation e.g. what gene annotation build the gene identifiers come from
      'Dataset Name' = dataset, # this column can be used to group images into datasets. These datasets will be used to group the images in the Image Data Resource. The name can be the same as the assay name.
      'Image File' = file,  # the name of the image file
      'Comment [Image File Path]' = path, # the path to the file
      'Comment [Image File Comments]' = '', # any comments about image files (can say its missing if you no longer have it)
      'Comment [Image File Type]' = type, # is it a raw image or derived?
      'Channels' = "Phase Contrast, GFP [epifluorescence]",  # the names of the channels and what was labeled in each channel
      'Processed Data File' = preproc, # name of the file with the results in
    ) %>% 
    
  write_delim('share/Julou_2020-lacInduction-assays.txt', delim='\t')

})()



# PLoS (lags export) #######

mycells_switching %>% ungroup %>% 
  filter_article_ds() %>% 
  filter(! date %in% discarded_dates) %>% 
  filter(!discard_arrested | condition=='switch_gly_lac') %>%
  filter(switch_idx == 1 | (str_detect(condition, '^switch_[0-9]+h$') & switch_idx < 3) ) %>% 
  left_join(
    mycells_switching_memory %>% 
      ungroup() %>% 
      mutate(gfp_inherit=parent_gfp/2^divs_since_par) %>% 
      select(ugen=ugen.x, gfp_inherit, parent_gfp, divs_since_par)
  ) %>% 
  left_join(select(mycells_switchingLow_memory, ugen, gfp_inherit, fluo_focus)) %>% 
  select(condition, ugen, date, pos, gl, switch_idx, lag_gfp_200=lag_200, growth_lag, 
         gfp_ini, gr_before=logl_time_slope_before, gfp_inherit, parent_gfp, n_divs_since_parent=divs_since_par, fluo_focus, time_birth, time_div) %>% 
  write.csv("share/Julou_2020-lacInduction-lags.csv", row.names=FALSE)
  # identity()


