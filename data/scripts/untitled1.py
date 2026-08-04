#!/bin/bash

folders=(
    "../Aves_Amphibia_trees/R_clade_parts"
    "../Aves_Amphibia_trees/P_Y_clade_parts"
    "../Aves_Amphibia_trees/B_clade_parts"
    "../Aves_Amphibia_trees/G_clade_parts"
)

for folder in "${folders[@]}"
do
    echo "Обработка папки: $folder"

    for file in "$folder"/*.fasta
    do
        if [ -f "$file" ]; then
            echo "Запуск IQ-TREE для: $file"

            iqtree \
                -s "$file" \
                -st DNA \
                -m MFP \
                -bb 1000 \
                -nt AUTO

        fi
    done
done

echo "Все выравнивания обработаны"