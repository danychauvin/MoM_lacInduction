(function(export_dir = "/scicore/home/nimwegen/GROUP/MM_Data/Thomas/_lacInductionArticle/Julou_2020_lacInduction_RawImages",
          .force=FALSE) {
# this anonymous function exports raw data (one dataset per compressed file) after the following scheme:
# - traverse R data to list all datasets used and corresponding files
# - create a directory per experiments and symlinks to corresponding files
# - create a compressed archive per dataset (following symlinks)
  
  # browser()
  
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
      fs::dir_ls(fs::path(., '..'))
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
  system2("find", c(export_dir, "-maxdepth 1 -mindepth 1 -type d -print0 | sort -z | xargs -r0 tree"), 
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
      sprintf("sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=tar_raw%s
#SBATCH --mem=32G
#SBATCH --time=1-0:00:00
#SBATCH --qos=1day
#SBATCH -o slogs/$JOB_NAME.o$JOB_ID
#SBATCH -e slogs/$JOB_NAME.e$JOB_ID
cd %s
tar --dereference -czf %s.tar.gz %s
rm -rf %s
EOF", .f, export_dir, .f, .f, .f) %>%
        system()
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
