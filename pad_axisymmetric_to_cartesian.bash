#! /bin/bash


# Loop over input arguments
for orig_file in "$@"
do
    new_file="padded_"$orig_file
    echo "$orig_file $new_file"
    awk '{if ($1=="ZONE"){}else{print $1" 0.0 "$2" "$3" "$4" "$5}}' $orig_file > $new_file
done
