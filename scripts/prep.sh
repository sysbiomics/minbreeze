#!/bin/bash


ASKFILE=$1
CACHE_DIR="$HOME/.sysmiome/minflow/cache/qiime2/"
QIIME2_REPO="https://qiime2-data.s3.amazonaws.com/2021.8/common/"
mkdir -p $CACHE_DIR

declare -A md5checksums
md5checksums=(
  ["gg-13-8-99-515-806-nb-classifier.qza"]="9e82e8969303b3a86ac941ceafeeac86"
  ["gg-13-8-99-515-806-nb-weighted-classifier.qza"]="8fb808c4af1c7526a2bdfaafa764e21f"
  ["gg-13-8-99-nb-classifier.qza"]="6bbc9b3f2f9b51d663063a7979dd95f1"
  ["gg-13-8-99-nb-weighted-classifier.qza"]="2baf87fce174c5f6c22a4c4086b1f1fe"
  ["sepp-refs-gg-13-8.qza"]="9ed215415b52c362e25cb0a8a46e1076"
  ["sepp-refs-silva-128.qza"]="7879792a6f42c5325531de9866f5c4de"
  ["silva-138-99-515-806-nb-classifier.qza"]="e05afad0fe87542704be96ff483824d4"
  ["silva-138-99-nb-classifier.qza"]="b8609f23e9b17bd4a1321a8971303310"
  ["silva-138-99-nb-weighted-classifier.qza"]="48965bb0a9e63c411452a460d92cfc04"
  ["silva-138-99-seqs-515-806.qza"]="a914837bc3f8964b156a9653e2420d22"
  ["silva-138-99-seqs.qza"]="de8886bb2c059b1e8752255d271f3010"
  ["silva-138-99-tax-515-806.qza"]="e2c40ae4c60cbf75e24312bb24652f2c"
  ["silva-138-99-tax.qza"]="f12d5b78bf4b1519721fe52803581c3d")

# Prep folder and things
CACHE_DIR="$HOME/.sysmiome/minflow/cache/qiime2/"
QIIME2_REPO="https://qiime2-data.s3.amazonaws.com/2021.8/common/"
mkdir -p $CACHE_DIR


CheckCache () {
  
  file=$1
  # Not on list, just exit
  if [ -z ${md5checksums[${file}]} ]
  then
	  echo "Not on the list"
  	  exit 1
  fi

  # File exists, check md5
  if [ -f ${CACHE_DIR}/${file} ]
  then
  	  expected_checksum=${md5checksums[$file]}
  	  computed_checksum=$(md5sum $CACHE_DIR/$file | awk '{print $1}')

      if [ "$expected_checksum" == "$computed_checksum" ]; then
        echo "$file: download successful: checksum matches."
        exit 0
      else
        echo "$file: error: checksum does not match."
        exit 1
      fi

  fi

  # File not exists download and check

  wget --no-check-certificate --no-proxy '${QIIME2_REPO}/file' -P $CACHE_DIR
  expected_checksum=${md5checksums[$file]}
  computed_checksum=$(md5sum $CACHE_DIR/$file | awk '{print $1}')

  if [ "$expected_checksum" == "$computed_checksum" ]; then
	echo "$file: download successful: checksum matches."
	exit 0
	else
	echo "$file: error: checksum does not match."
	exit 1
  fi
}

CheckCache $ASKFILE


# Loop through the files and checksums
#for file in "${!md5checksums[@]}"; do

    # Compute the checksum of the downloaded file
#    expected_checksum=${md5checksums[$file]}
#    computed_checksum=$(md5sum $CACHE_DIR/$file | awk '{print $1}')

    # Compare the checksums and print a message
#    if [ "$expected_checksum" == "$computed_checksum" ]; then
#        echo "$file: download successful: checksum matches."
#    else
#        echo "$file: error: checksum does not match."
#    fi
#done
