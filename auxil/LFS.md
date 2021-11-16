# Uploading large files using lfs

## Instructions on installing lfs
- https://github.com/git-lfs/git-lfs?utm_source=gitlfs_site&utm_medium=source_link&utm_campaign=gitlfs

## Cloning a Repository

Before you can edit your repository you will need to clone/download it.    

    git clone <reponame>.git

## Adding new files

You start by tracking file type (this is the extention), if it is already being track you do not have to track it again. Your list of tracked types is in .gitattributes.    

Example .gitattributes tracking all csv files:    

    *.csv filter=lfs diff=lfs merge=lfs -text

Once the file type is being tracked you just have to add it:    

    git add <filename>

You can verify it is being tracked by typing:    

    git lfs ls-files
       
Once it has been added you need to commit it then push it to your repo:    

    git commit -m "Notes you want included along with you commit"
    git push

## If the filetype is not added then you will need to add it    

    git lfs track "*.psd"
