#!/bin/bash


extract_fl () {
	local zipfile=$1
	local lookfor=$2
	local location_z=$(unzip -l $zipfile  | gawk '{print $4}' | grep ${lookfor})  # Really fraggile, but work

	unzip -j $zipfile $location_z
}

extract_fl rooted-tree.qza tree.nwk
mv tree.nwk root-tree.nwk
extract_fl taxonomy_gg.qza taxonomy.tsv
mv taxonomy.tsv taxonomy_gg.tsv

