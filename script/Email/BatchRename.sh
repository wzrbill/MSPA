#!/bin/bash
# Example usage
# bash BatchRename.sh Email Twitch
# All key word will be replace from Email to Twitch
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
