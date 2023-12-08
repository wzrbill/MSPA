#!/bin/bash
# example usage
# bash BatchRename NY SH
# And all key word will be replaced from NY to SH
fromName=$1
toName=$2

for file in ./*; do
    if [[ "$file" == *"$fromName"* ]]; then
        sed -i "s/${fromName}/${toName}/g" "$file"
        new_name="${file//$fromName/$toName}"
        
        mv "$file" "$new_name"
        
        echo "Renamed: $file -> $new_name"
    fi
done
