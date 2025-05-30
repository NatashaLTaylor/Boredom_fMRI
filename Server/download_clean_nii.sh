#!/bin/bash

#Script copies all relevant clean data from the HCP database for a list of subjects to University server

# curl function to download relevant HCP data
DBHOST=https://db.humanconnectome.org
PROJECT=HCP_1200
REST_URL_PREFIX=$DBHOST/data/archive/projects/$PROJECT/subjects
Folder_rfMRI=rfMRI_REST1_LR #change to rfMRI_REST1_RL
##################################################################################################################
# GET a jsession
# MY_CURL_OPTIONS: Options that you may want to pass to curl such as --connect-timeout, --retry etc
db_user=### #Username for the database
db_password=### #Password for the database
subject_list_file=subject1.txt
# Convert the input file to Unix-style line endings
dos2unix "$subject_list_file" 2>/dev/null
#curl -u user # https://db.humanconnectome.org/data/archive/projects/HCP_1200/subjects/100307/experiments/100307_CREST/resources/100307_CREST/files/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR.nii.gz -O 
jsession=`curl -u $db_user:$db_password $DBHOST`

while IFS= read -r line; do
      subject=$(echo "$line" | xargs)
      mkdir "$subject"
      cd "$subject"
      subject_url_prefix=$REST_URL_PREFIX/${subject}/experiments/${subject}_CREST/resources/${subject}_CREST/files/
      file_relative_path1=MNINonLinear/Results/$Folder_rfMRI/${Folder_rfMRI}_hp2000_clean.nii.gz
      file_relative_path2=MNINonLinear/Results/$Folder_rfMRI/${Folder_rfMRI}.nii.gz

      #curl -u $db_user:$db_password -O $subject_url_prefix$file_relative_path1 
      curl -u $db_user:$db_password -O "$subject_url_prefix$file_relative_path1"
      curl -u $db_user:$db_password -O "$subject_url_prefix$file_relative_path2"
      echo $subject 
      cd ..
done < $subject_list_file